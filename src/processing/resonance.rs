use anyhow::{Context, Result, bail};
use getset::{Getters, Setters};
use itertools::Itertools;
use cascade::cascade;
use rustfft::num_complex::Complex;
use serde::Deserialize;

use crate::processing;
use crate::processing::{
    Constraint, ConstraintLeeways, Bounds, Parameters, 
    PeaksTemplate, Peaks
};
use crate::fitting::{f_peak, f_peak_area};


//==============================================================================
// Resonance

#[derive(Clone, Debug, Default, Getters, Setters)]
#[getset(get = "pub")]
pub struct Resonance {

    #[getset(set = "pub")]
    name: String,

    peaks: Peaks,
    
    #[getset(set = "pub")]
    peak_names: Vec<String>,

    #[getset(set = "pub")]
    parameters: Parameters,
    #[getset(set = "pub")]
    constraint_leeways: ConstraintLeeways,
    #[getset(set = "pub")]
    general_bounds: Bounds,
    #[getset(set = "pub")]
    offset_bounds: Bounds,
}


//=============================================================================
// YAML parsing

impl Resonance {

    //-------------------------------------------------------------------------
    pub fn from_template(
        template: ResonanceTemplate
    ) -> Result<Self> {

        // Ensure that there are no references
        let peaks_template = template.peaks.clone();

        match peaks_template {
            PeaksTemplate::Label(_) => {
                bail!("Invalid reference to another resonance")
            },
            _ => {}
        };

        // Peaks are converted first
        let peaks = Peaks::from_template(peaks_template)
            .context("Invalid peak definition")?;

        let width = match template.width {
            Some(value) => value,
            None => f64::NAN,
        };

        let height = match template.height {
            Some(value) => value,
            None => f64::NAN,
        };

        let fraction = match template.fraction {
            Some(value) => value,
            None => f64::NAN,
        };

        let n_peaks = peaks.n_peaks();

        let mut resonance = cascade!{
            Resonance::from_peaks(peaks);
                ..set_name(template.resonance);
                ..set_width(width);
                ..set_height(height);
                ..set_fraction(fraction);
                ..set_constraint_leeways(template.constraint_leeways);
                ..set_general_bounds(template.general_bounds);
                ..set_offset_bounds(template.offset_bounds);
        };

        // Rename peaks if applicable
        let n_names = template.peak_names.len();

        if n_names > 0 {
            if n_names != n_peaks {
                bail!("Invalid number of peak names")
            }

            resonance.set_peak_names(template.peak_names);
        }

        Ok( resonance )
    }
}


//=============================================================================
// Constructors

impl Resonance {

    //-------------------------------------------------------------------------
    pub fn from_pattern(pattern: String) -> Result<Resonance> {

        let name = pattern.clone();
        let peaks = Peaks::from_pattern(pattern.clone())
            .with_context(|| {
                format!("Could not generate Resonance from \"{}\".", &pattern)
            })?;

        let resonance = Resonance {
            name: name,
            peaks: peaks,
            ..Default::default()
        };

        Ok(resonance)
    }


    //-------------------------------------------------------------------------
    pub fn from_positions(
        positions: Vec<f64>, 
        areas: Vec<f64>
    ) -> Result<Resonance> {

        let n = positions.len();
        let name = format!("{}..{} m", positions[0], positions[n]);
        let peaks = Peaks::from_positions(positions, areas)
            .context("Could not generate Resonance.")?;

        let mut resonance = Resonance {
            name: name,
            peaks: peaks,
            ..Default::default()
        };

        let peak_names: Vec<String> = (0 .. resonance.n_peaks())
            .map(|x| format!("{}", x + 1))
            .collect();

        resonance.set_peak_names(peak_names);

        Ok(resonance)
    }


    //-------------------------------------------------------------------------
    pub fn from_peaks(peaks: Peaks) -> Resonance {

        let mut resonance = Resonance {
            peaks: peaks,
            ..Default::default()
        };

        let peak_names: Vec<String> = (0 .. resonance.n_peaks())
            .map(|x| format!("{}", x + 1))
            .collect();

        resonance.set_peak_names(peak_names);

        resonance
    }
}


//=============================================================================
// Simple getters/setters

impl Resonance {

    //-------------------------------------------------------------------------
    pub fn set_width(&mut self, width: f64) {
        self.parameters.width = width;
    }

    pub fn set_height(&mut self, height: f64) {
        self.parameters.height = height;
    }

    pub fn set_fraction(&mut self, fraction: f64) {
        self.parameters.fraction = fraction;
    }
}


//=============================================================================
// Methods

impl Resonance {

    //-------------------------------------------------------------------------
    pub fn is_in(&self, frequency: f64, ppm: &[f64; 2]) -> bool {

        self.terms(frequency)
            .iter()
            .step_by(4)
            .any(|x| (x > &ppm[0]) && (x < &ppm[1]))
    }

    //-------------------------------------------------------------------------
    pub fn n_peaks(&self) -> usize {

        self.peaks.n_peaks()
    }

    //-------------------------------------------------------------------------
    pub fn terms(&self, frequency: f64) -> Vec<f64> {

        self.peaks.terms(frequency)
    }


    //-------------------------------------------------------------------------
    pub fn update(&mut self, frequency: f64, defaults: &Parameters) {

        let mut parameters = self.parameters.clone();
        parameters.update(defaults);

        self.peaks.update(frequency, &parameters);
    }


    //-------------------------------------------------------------------------
    pub fn replace(&mut self, terms: &[f64], frequency: f64) {

        self.peaks.replace(terms, frequency);
    }


    //-------------------------------------------------------------------------
    pub fn bounds(&self, frequency: f64) -> Vec<[f64; 2]> {

        let mut offset_bounds = self.offset_bounds().bounds();
        offset_bounds[1][0] /= frequency;
        offset_bounds[1][1] /= frequency;
        let offset_bounds = IntoIterator::into_iter(offset_bounds);

        let mut general_bounds = self.general_bounds().bounds();
        general_bounds[1][0] /= frequency;
        general_bounds[1][1] /= frequency; 
        let general_bounds = IntoIterator::into_iter(general_bounds);

        // Generate offset bounds then update with general ones
        self.terms(frequency)
            .into_iter()
            .zip(offset_bounds.cycle())
            .map(|(x, [l, u])| [x+l, x+u])
            .zip(general_bounds.cycle())
            .map(|([l1, u1], [l2, u2])| {
                [f64::max(l1, l2), f64::min(u1, u2)]
            })
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn constraints(
        &self, 
        frequency: f64,
        offset: usize,
    ) -> Vec<Constraint> {

        let (parameters, areas) = self.peaks.extract(frequency)
            .expect("Coupling could not be converted into peak parameters.");

        let mut constraints: Vec<Constraint> = Vec::new();

        // Each pair of peaks must have constant position offset, constant
        // height ratio, equal width, and equal fraction.
        // Every peak is referenced to the first.
        let leeways = &self.constraint_leeways;

        let p0 = parameters[0].position;
        let a0 = areas[0];
        let j0 = 0 + offset;

        for i in 1..areas.len() {
            let p = parameters[i].position - p0;
            let a = areas[i]/a0;
            let j1 = i + offset;

            constraints.extend(
                vec![
                    Constraint::position((j0, j1), p, leeways.position),
                    Constraint::width((j0, j1), leeways.width),
                    Constraint::height((j0, j1), a, leeways.height),
                    Constraint::fraction((j0, j1), leeways.fraction) 
                ]
            );
        }

        constraints
    }


    //-------------------------------------------------------------------------
    pub fn peak_lineshapes(
        &self, 
        x: &[f64], 
        frequency: f64
    ) -> Vec<Vec<Complex<f64>>> {

        self.terms(frequency)
            .chunks(4)
            .map(|p| {
                let mut y = vec![0.0; x.len()*2];
                f_peak(&p[..], &x[..], &mut y[..], 1e-10);

                y.into_iter()
                    .chunks(2)
                    .into_iter()
                    .map(|mut c| {
                        Complex::new(c.next().unwrap(), c.next().unwrap())
                    })
                    .collect()
            })
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn peak_areas(&self, frequency: f64) -> Vec<f64> {

        self.terms(frequency)
            .chunks(4)
            .map(|y| f_peak_area(&y[..], 1e-10) )
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn resonance_lineshape(
        &self, 
        x: &[f64], 
        frequency: f64
    ) -> Vec<Complex<f64>> {

        let lineshapes = self.peak_lineshapes(x, frequency);
        processing::sum_across(lineshapes)
    }


    //-------------------------------------------------------------------------
    pub fn resonance_area(&self,  frequency: f64) -> f64 {

        self.peak_areas(frequency).iter().sum()
    }


}


//=============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct ResonanceTemplate {

    resonance: String,

    peaks: PeaksTemplate,
    #[serde(default)] 
    peak_names: Vec<String>,

    #[serde(default)] 
    constraint_leeways: ConstraintLeeways,

    #[serde(default)] 
    width: Option<f64>,
    #[serde(default)] 
    height: Option<f64>,
    #[serde(default)] 
    fraction: Option<f64>,

    #[serde(default)] 
    general_bounds: Bounds,
    #[serde(default)] 
    offset_bounds: Bounds,
}
