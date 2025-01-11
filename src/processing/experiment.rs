use anyhow::{Context, Result, bail};
use getset::{Getters, MutGetters, Setters};
use itertools::izip;
use itertools_num;
use nlopt::*;
use rustfft::{FftPlanner, num_complex::Complex};
use std::{
    iter,
    f64::consts::PI,
    path::PathBuf,
    collections::HashMap
};

use crate::importing;
use crate::importing::ImportingError;
use crate::processing::{
    UNBOUNDED, Baseline1D, Bound, Bounds, Constraint, Phase1D, Species1D
};


//=============================================================================
// Experiment1D 

#[derive(Clone, Debug, Getters, MutGetters, Setters)]
#[getset(get = "pub")]
pub struct Experiment1D {

    path: PathBuf,

    #[getset(set = "pub")]
    name: String,
    #[getset(set = "pub")]
    frequency: f64,
    #[getset(set = "pub")]
    sampling_frequency: f64,

    chemical_shift: Vec<f64>,
    intensity: Vec<Complex<f64>>,
    species: HashMap<String, Species1D>,

    #[getset(get = "pub", get_mut = "pub")]
    baseline: Baseline1D,
    #[getset(get = "pub", get_mut = "pub")]
    phase: Phase1D,
}


impl PartialEq for Experiment1D {
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
    }
}

impl Eq for Experiment1D {}


//=============================================================================
// Constructors

impl Experiment1D {

    pub fn new(
        path: PathBuf,
        frequency: f64,
        sampling_frequency: f64,
        chemical_shift: Vec<f64>, 
        intensity: Vec<Complex<f64>>,
    ) -> Result<Experiment1D> {

        // Get name from path
        let name = path.file_name()
            .expect("Error extracting experiment name.")
            .to_os_string()
            .into_string()
            .expect("Error parsing experiment name.");

        // Basic validity check
        if chemical_shift.len() != intensity.len() {
            bail!("Chemical shift and intensity must be of equal length.")
        }

        // TODO: generalize default parameters for baseline/phase

        // Initializing default baseline
        let baseline = Baseline1D::from_span(
            &chemical_shift, f64::NAN, 3
        )?;

        // Initializing default phase
        let phase = Phase1D::new(1)?;

        let mut experiment = Experiment1D {
            path: path,
            name: name,
            frequency: frequency,
            sampling_frequency: sampling_frequency,
            chemical_shift: chemical_shift,
            intensity: intensity,
            species: HashMap::new(),
            baseline: baseline,
            phase: phase,
        };

        experiment.scale();

        Ok( experiment )
    }

    //-------------------------------------------------------------------------
    pub fn from_path(
            path: &PathBuf, 
            processing: &Option<String>
        ) -> Result<Experiment1D> {

        type Loader = 
            fn(&PathBuf, &Option<String>) -> 
                Result<Experiment1D, ImportingError>;

        // List of import functions to try
        let loaders: Vec<(&str, Loader)> = vec![
            ("RS2D", importing::import_rs2d_processed_1d),
            ("Bruker", importing::import_bruker_1d)
        ];

        // Try import functions one at a time
        for (name, f) in loaders {

            match f(path, processing) {

                // If data parsed, scale and return
                Ok(experiment) => return Ok(experiment),

                // If mismatch, ignore. Otherwise error out
                Err(ImportingError::Mismatch) => {},

                Err(ImportingError::MalformedData(error)) => {
                    let out = Err(error).with_context(|| {
                        format!(
                            "Data in \"{}\" matched {} file format, \
                            but encountered error during import.",
                            path.to_string_lossy(),
                            name,
                        )
                    });

                    return out
                }
            }
        }

        bail!(
            "Data in \"{}\" could not be imported using available methods.", 
            path.to_string_lossy()
        )
    }


}


//=============================================================================
// Setters

impl Experiment1D {


    //-------------------------------------------------------------------------
    pub fn zero_fill(&mut self, factor: f64) {

        if factor == 0.0 {
            return
        }

        let n = self.intensity.len();
        let n_zeros = ((n as f64)*factor).ceil() as usize;
        let n_new = n + n_zeros;

        // Take FFT
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        fft.process(&mut self.intensity);

        // Add zeros
        self.intensity = self.intensity.clone()
            .into_iter()
            .chain(iter::repeat(Complex::new(0.0, 0.0)).take(n_zeros))
            .collect();

        // Invert
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_inverse(n_new);
        fft.process(&mut self.intensity);

        // Modify chemical shift to match
        let lower = self.chemical_shift[0];
        let upper = self.chemical_shift[n-1];
        self.chemical_shift = itertools_num::linspace(lower, upper, n_new)
            .collect();
    }


    //-------------------------------------------------------------------------
    pub fn line_broaden(&mut self, width: f64) {

        if width == 0.0 {
            return
        }

        let n = self.intensity.len();

        // Take FFT
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        fft.process(&mut self.intensity);

        // Apodize
        let scalar: f64 = -width*PI/self.sampling_frequency*2.0;

        for (i, value) in self.intensity.iter_mut().enumerate() {
            *value *= (scalar * (i as f64)).exp();
        }

        // Invert
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_inverse(n);
        fft.process(&mut self.intensity);

        // Re-scale
        self.scale();
    }


    //-------------------------------------------------------------------------
    pub fn set_baseline(
            &mut self, span: f64, degree: usize, bounds: Bound
        ) -> Result<()> {

        let x = &self.chemical_shift;

        let mut baseline = Baseline1D::from_span(x, span, degree)?;
        baseline.set_bounds(bounds);

        self.baseline = baseline;

        Ok( () )
    }


    //-------------------------------------------------------------------------
    pub fn set_phase(
            &mut self, order: usize, bounds: f64
        ) -> Result<()> {

        let mut phase = Phase1D::new(order)?;
        phase.set_bounds(bounds);

        self.phase = phase;

        Ok( () )
    }
}


//=============================================================================
// Basic processing

impl Experiment1D {

    //-------------------------------------------------------------------------
    pub fn scale(&mut self) {

        // Scale y data to a max real height of 1
        let max: f64 = self.intensity.iter()
            .map(|y| y.re)
            .fold(std::f64::NEG_INFINITY, |max, y| max.max(y));

        for y in self.intensity.iter_mut() {
            *y /= max;
        }
    }


    //-------------------------------------------------------------------------
    pub fn autophase(&mut self) {

        // Define objective function
        fn f_obj(
                x: &[f64], 
                grad: Option<&mut [f64]>, 
                y: &mut Vec<Complex<f64>>
            ) -> f64 {

            let factor = Complex::new(0.0, x[0]).exp();

            let threshold = 0.05 * y.iter()
                .fold(f64::NEG_INFINITY, |a, b| {
                    a.max(b.re)
                }); 

            y.iter()
                .filter(|y| y.re.abs() >= threshold)
                .fold(0.0, |a, b| {
                    let Complex { re, .. } = b * factor;
                    a + re
                })
        }

        // Initial phase guess is 0
        let mut x = vec![0.0];
        let xl = vec![-PI];
        let xu = vec![PI];

        let y = self.intensity.clone();

        let mut opt = Nlopt::new(
            Algorithm::Bobyqa, 1, f_obj, Target::Maximize, y
        );

        opt.set_lower_bounds(&xl).unwrap();
        opt.set_upper_bounds(&xu).unwrap();
        opt.set_xtol_rel(1e-6).unwrap();

        let res = opt.optimize(&mut x);

        if res.is_ok() {
            let factor = Complex::new(0.0, x[0]).exp();
            for y in self.intensity.iter_mut() {
                *y *= factor;
            }

            // Re-scale
            self.scale();
        }
    }


    //-------------------------------------------------------------------------
    pub fn reference_align(&mut self, ppm: f64, threshold: f64) {

        if ppm.is_nan() {
            return
        }

        // Smooth using SG
        let norm: Vec<_> = self.intensity.iter()
            .map(|y| y.norm())
            .collect();

        let values: Vec<_> = norm.windows(7)
            .map(|y| {
                (   5.0*(y[0] + y[6]) -  
                    30.0*(y[1] + y[5]) +
                    75.0*(y[2] + y[4]) + 
                    131.0*y[3]        
                )/231.0
            })
            .collect();

        let value_threshold = threshold * values.iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();        

        let derivatives: Vec<_> = norm.windows(7)
            .map(|y| {
                (   22.0*(y[0] - y[6]) +  
                    67.0*(-y[1] + y[5]) +
                    58.0*(-y[2] + y[4])
                )/252.0
            })
            .collect();

        let derivative_threshold = threshold * derivatives.iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

        let iterator = izip!(
            values.into_iter(),
            derivatives.into_iter(),
            self.chemical_shift.iter().skip(3)
        );

        let mut offset = iterator.rev()
            .skip_while(|(v, d, p)| {
                (*v <= value_threshold) || (d.abs() >= derivative_threshold)
            })
            .map(|(_, _, p)| *p)
            .next().
            unwrap_or(0.0);

        offset -= ppm;

        for x in self.chemical_shift.iter_mut() {
            *x -= offset;
        }
    }
}


//=============================================================================
// Data (and phasing)

impl Experiment1D {

    //-------------------------------------------------------------------------
    pub fn domain(&self) -> [f64; 2] {

        // Figure out min/max values of data
        let mut xmin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;

        for x in self.chemical_shift.iter() {
            xmin = xmin.min(*x);
            xmax = xmax.max(*x);
        }

        return [xmin, xmax]
    }

    //-------------------------------------------------------------------------
    pub fn phase_values(&self, x: &[f64]) -> Vec<f64> {

        let terms = self.phase.terms();
        let b = terms[0];
        
        let mut a: f64 = 0.0;
        if terms.len() > 1 {
            a = terms[1];
        }

        // Apply to given values of x
        x.iter()
            .map(|x| a*x + b)
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn data(
        &self, 
        bounds: Option<[f64; 2]>,
        phase: bool,
    ) -> (Vec<f64>, Vec<Complex<f64>>) {

        // Get x,y data
        let mut x = self.chemical_shift.clone();
        let mut y = self.intensity.clone();

        // Filter
        if let Some(bounds) = bounds {
            (x, y) = x 
                .into_iter()
                .zip(y)
                .filter(|(x, _)| (x > &bounds[0]) && (x < &bounds[1]))
                .unzip();
        }

        // Apply phase
        if phase {
            let mut theta = self.phase_values(&x[..]);

            for (y, theta) in y.iter_mut().zip(theta.iter()) {

                let y_r = y.re;
                let y_i = y.im;

                let sin = theta.sin();
                let cos = theta.cos();

                y.re =  y_r*cos + y_i*sin;
                y.im = -y_r*sin + y_i*cos;
            }
        }

        // And output
        (x, y)
    }
}


//=============================================================================
// Species

impl Experiment1D {

    //-------------------------------------------------------------------------
    pub fn add_species(&mut self, species: &Vec<Species1D>) {
        
        for species in species.iter() {

            // Initialize heights based on experimental data
            let mut species = species.clone();
            let mut terms = species.terms(self.frequency, &self.domain());

            // Width around peak is 2 Hz to either side
            let width = 2.0/self.frequency;

            let positions = terms.clone().into_iter()
                .step_by(4);

            let heights = terms.iter_mut()
                .skip(2)
                .step_by(4);

            for (position, height) in positions.zip(heights) {

                if ! height.is_nan() {
                    continue
                }
                
                let x_low = position - width;
                let x_high = position + width;
                
                // Assuming that chemical shift is decreasing
                *height = self.chemical_shift.iter()
                    .zip(self.intensity.iter())
                    .skip_while(|(x, _)| **x > x_high)
                    .take_while(|(x, _)| **x > x_low)     
                    .map(|(_, y)| {
                        y.re
                    })
                    .max_by(|a, b| a.partial_cmp(b).expect("Sort error."))
                    .expect("Error estimating height.");
            }

            species.replace(&terms, self.frequency, &UNBOUNDED); 
            self.species.insert(species.name().clone(), species);
        }
    }


    //-------------------------------------------------------------------------
    pub fn n_peaks(&self, species: &Vec<String>, ppm: &[f64; 2]) -> usize {
        
        species.iter()
            .map(|x| {
                self.species()
                    .get(x)
                    .expect("Missing species.")
                    .n_peaks(self.frequency, ppm)
            })
            .sum()
    }


    //-------------------------------------------------------------------------
    pub fn terms(&self, species: &Vec<String>, ppm: &[f64; 2]) -> Vec<f64> {

        let terms: Vec<f64> = species.iter()
            .map(|x| {
                self.species()
                    .get(x)
                    .expect("Missing species.")
                    .terms(self.frequency, ppm) 
            })
            .flatten()
            .collect();

        terms
    }

    //-------------------------------------------------------------------------
    pub fn replace(
        &mut self, 
        species: &Vec<String>, 
        terms: &[f64],
        ppm: &[f64; 2]
    ) {

        if terms.len()/4 != self.n_peaks(species, ppm) {
            panic!("Invalid number of terms.")
        }

        let mut i = 0;

        for s in species {
            let s = self.species.get_mut(s)
                .expect("Missing species.");

            let j = i + s.n_peaks(self.frequency, ppm)*4;
            s.replace(&terms[i..j], self.frequency, ppm);
            i = j;
        }
    }


    //-------------------------------------------------------------------------
    pub fn bounds(
        &self, 
        species: &Vec<String>, 
        ppm: &[f64; 2], 
        offset_bounds: &Bounds,
        general_bounds: &Bounds
    ) -> Vec<[f64; 2]> {

        let mut offset_bounds = offset_bounds.bounds();
        offset_bounds[1][0] /= self.frequency;
        offset_bounds[1][1] /= self.frequency;
        let offset_bounds = IntoIterator::into_iter(offset_bounds);

        let mut general_bounds = general_bounds.bounds();
        general_bounds[1][0] /= self.frequency;
        general_bounds[1][1] /= self.frequency; 
        let general_bounds = IntoIterator::into_iter(general_bounds);

        // Generate offset bounds then update with general ones
        let bounds: Vec<[f64; 2]> = self.terms(species, ppm)
            .into_iter()
            .zip(offset_bounds.cycle())
            .map(|(x, [l, u])| [x+l, x+u])
            .zip(general_bounds.cycle())
            .map(|([l1, u1], [l2, u2])| {
                [f64::max(l1, l2), f64::min(u1, u2)]
            })
            .collect();

        // Update species bounds with defaults
        species.iter()
            .map(|x| {
                self.species()
                    .get(x)
                    .expect("Missing species.")
                    .bounds(self.frequency, ppm) 
            })
            .flatten()
            .zip(bounds.into_iter())
            .map(|([l1, u1], [l2, u2])| {
                [f64::max(l1, l2), f64::min(u1, u2)]
            })
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn constraints(
        &self, 
        species: &Vec<String>, 
        ppm: &[f64; 2]
    ) -> Vec<Constraint> {

        let mut offset = 0;

        let mut constraints: Vec<Constraint> = Vec::new();

        // Pull out species constraints
        for s in species.iter() {
            
            let species = self.species.get(s)
                .expect("Missing species.");

            constraints.extend(species.constraints(self.frequency, ppm, offset));
            offset += species.n_peaks(self.frequency, ppm);
        }

        constraints
    }


    //-------------------------------------------------------------------------
    pub fn species_lineshape(
        &self, 
        species: &String, 
        x: &[f64]
    ) -> Vec<Complex<f64>> {
        
        //---
        // If species does not exist, output zeros
        //---
        let species = self.species().get(species);

        if species.is_none() {
            let zero = Complex::new(0.0, 0.0);
            let n = x.len();
            return iter::repeat(zero).take(n).collect()
        }

        //---
        // If species does exist, ensure that all NaN terms are replaced
        // with defaults.
        //---
        
        let species = species.unwrap();

        return species.species_lineshape(x, self.frequency)
    }
}

/*
//=============================================================================
// Baselines

impl Experiment1D {

    //-------------------------------------------------------------------------
    pub fn add_baseline(&mut self, baseline: Baseline1D) {
        self.baselines.push(baseline);
    }


    //-------------------------------------------------------------------------
    pub fn baseline(&self, x: &[f64]) -> Vec<Complex<f64>> {

        let zero = Complex::<f64>::new(0.0, 0.0);
        let n = x.len();

        let mut xmin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;
        for x in x.iter() {
            xmin = xmin.min(*x);
            xmax = xmax.max(*x);
        }

        // If there are no baselines, output zero
        if self.baselines.len() == 0 {
            return iter::repeat(zero).take(n).collect()
        }

        // Get all baselines and then merge them
        // Ignoring non-applicable baselines
        let mut baselines: Vec<_> = self.baselines
            .iter()
            .filter(|y| {
                (xmin < y.ppm()[1]) && (xmax > y.ppm()[0])
            })
            .map(|y| y.baseline(x).into_iter())
            .collect();

        // If there are no baseline within bounds, output zero
        if baselines.len() == 0 {
            return iter::repeat(zero).take(n).collect()
        }

        // Transposing baseline vector to have outer vector correspond to points
        let baselines: Vec<Vec<Complex<f64>>> = (0 .. n)
            .map(|_| {
                baselines 
                    .iter_mut()
                    .map(|b| b.next().unwrap())
                    .collect::<Vec<Complex<f64>>>()
            })
            .collect();

        // Weight is degree + 1
        let weights: Vec<f64> = self.baselines
            .iter()
            .map(|y| (2 + y.np() - y.nk()) as f64)
            .collect();

        let mut baseline: Vec<Complex<f64>> = Vec::new();

        for i in 0 .. n {
            
            let (values, weights): (Vec<Complex<f64>>, Vec<f64>) = baselines[i]
                .iter()
                .zip(weights.iter())
                .filter(|(y, _)| y != &&zero)
                .map(|(y, z)| (*y, *z))
                .unzip();

            if values.len() == 0 {
                baseline.push(zero.clone());
            } else {
                let denominator: f64 = weights.iter().sum();
                let numerator: Complex<f64> = values
                    .into_iter()
                    .zip(weights.into_iter())
                    .map(|(x, w)| x*w)
                    .sum();

                baseline.push(numerator/denominator);
            }
        }

        baseline
    }
}

*/
