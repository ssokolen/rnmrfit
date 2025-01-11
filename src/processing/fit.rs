use anyhow::{bail, Result};
use cascade::cascade;
use getset::{Getters, MutGetters, Setters};
use itertools::{Itertools, izip};
use ndarray::prelude::*;
use nlopt::{Algorithm, Target, Nlopt};
use rustfft::num_complex::Complex;
use serde::Deserialize;
use std::{
    f64::consts::PI,
    iter,
    ops::Range,
};

use crate::processing::*;
use crate::fitting;



//=============================================================================
// Misc structs

#[derive(Clone, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Output {
    Plot,
    Area,
    Lineshape,
}


//=============================================================================
// Fit1D 

#[derive(Clone, Debug, Default, Getters, MutGetters, Setters)]
#[getset(get = "pub")]
pub struct Fit1D {

    #[getset(set = "pub")]
    name: String,

    #[getset(set = "pub")]
    species: Vec<String>,
    #[getset(set = "pub")]
    experiments: Vec<String>,

    #[getset(set = "pub")]
    ppm: [f64; 2],

    #[getset(set = "pub")]
    parameters: Parameters,
    #[getset(set = "pub")]
    constraint_leeways: ConstraintLeeways,
    #[getset(set = "pub")]
    general_bounds: Bounds,
    #[getset(set = "pub")]
    offset_bounds: Bounds,

    #[getset(set = "pub")]
    outputs: Vec<Output>,
}


//=============================================================================
// YAML parsing

impl Fit1D {

    //-------------------------------------------------------------------------
    pub fn from_template(
        template: Fit1DTemplate,
        job: &Job1D
    ) -> Result<Self> {

        // Starting with name
        let name = template.fit;

        // Ensuring that species and experiments exist
        let species = template.species.patterns()?;
        job.check_species_globs(&species)?;

        let species = job.get_species_from_globs(&species);

        let experiments = template.experiments.patterns()?;
        job.check_experiment_globs(&experiments)?;

        let experiments = job.get_experiments_from_globs(&experiments);

        let [ppm1, ppm2] = template.ppm;
        let lower = ppm1.min(ppm2);
        let upper = ppm1.max(ppm2);

        // Use a not to catch NAN values
        if ! ((upper - lower) > 1e-4) {
            bail!("Invalid ppm range or too small \"{} - {}\"", ppm1, ppm2)
        }

        let ppm = [lower, upper];

        // Initializing
        let mut fit = cascade!{
            Fit1D::new(species, experiments);
                ..set_name(name);
                ..set_ppm(ppm);
                ..set_width(template.width);
                ..set_height(template.height);
                ..set_fraction(template.fraction);
                ..set_constraint_leeways(template.constraint_leeways);
                ..set_general_bounds(template.general_bounds);
                ..set_offset_bounds(template.offset_bounds);
                ..set_outputs(template.outputs);
        };

        Ok( fit )
    }
}


//=============================================================================
// Constructors

impl Fit1D {

    //-------------------------------------------------------------------------
    pub fn new(
        species: Vec<String>,
        experiments: Vec<String> 
    ) -> Fit1D {

        Fit1D {
            species: species,
            experiments: experiments,
            ..Default::default()
        }
    }
}


//=============================================================================
// Simple getters/setters

impl Fit1D {

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
// Fitting

impl Fit1D {

    //-------------------------------------------------------------------------
    pub fn fit(&mut self, job: &mut Job1D, tol: f64) -> Result<()> {

        // Get species defaults
        let species: Vec<Species1D> = self.species.iter()
            .map(|x| {
                let mut species = job.species()
                    .get(x)
                    .expect("Missing species.")
                    .clone();

                // Update bounds
                let mut general_bounds = species.general_bounds().clone();
                general_bounds.update(&self.general_bounds);
                species.set_general_bounds(general_bounds);

                let mut offset_bounds = species.offset_bounds().clone();
                offset_bounds.update(&self.offset_bounds);
                species.set_offset_bounds(offset_bounds);

                species
            })
            .collect();

        // Extracting ranges
        let ranges = job.ranges().clone();

        // Loop through experiments
        for experiment in self.experiments.iter() {

            // Get experimental data from name
            let experiment = job.experiments_mut()
                .get_mut(experiment)
                .expect("Missing experiment.");

            // Initialize and update species with defaults
            experiment.add_species(&species);

            // Extracting data
            let (mut x, y) = experiment.data(None, false);
            let n = x.len();

            //---
            // Extract terms and associated bounds
            //---

            let nl = experiment.n_peaks(&self.species, &self.ppm)*4;
            let mut terms = experiment.terms(&self.species, &self.ppm);

            if terms.len() == 0 {
                let error = SetupError::NoPeaks(self.ppm[0], self.ppm[1]);
                bail!(SetupError::error(error.into(), &self.name));
            }
            
            let general_bounds = self.general_bounds.clone();
            let offset_bounds = self.offset_bounds.clone();

            let bounds = experiment.bounds(
                    &self.species, &self.ppm, &offset_bounds, &general_bounds
                )
                .into_iter()
                .map(|x| (x[0], x[1]))
                .unzip();

            type Split = (Vec<f64>, Vec<f64>);
            let (mut lower_bounds, mut upper_bounds): Split = bounds;

            //---
            // Calculate ranges and split peaks
            //---

            let ranges = match ranges {
                Ranges::Manual(ref manual) => {
                    Fit1D::calculate_manual_ranges(
                        &x, &lower_bounds, &upper_bounds, manual.cutoff
                    )
                }
            };

            // It is also useful to have range values
            let range_values: Vec<_> = ranges.iter()
                .map(|range| {
                    [x[range.end - 1], x[range.start]]
                })
                .collect();

            // Assign each peak to a range.
            // (Note that peaks may not be sorted)

            // TODO: Check for NANs
            let positions: Vec<_> = terms.iter()
                .step_by(4)
                .collect();

            // Peaks are assigned to ranges one at a time
            let assignments: Vec<_> = positions.iter()
                .map(|p| {
                    range_values.iter()
                        .enumerate()
                        .find(|(_, r)| (**p >= r[0]) && (**p <= r[1]))
                        .expect("Error generating ranges")
                        .0
                })
                .collect();

            //---
            // Trim x/y data to just ranges
            //---

            let mut x_trimmed: Vec<f64> = Vec::new();
            let mut y_trimmed: Vec<Complex<f64>> = Vec::new();

            let mut offset = 0;
            let ranges = ranges.into_iter()
                .map(|mut range| {
                    x_trimmed.extend(&x[range.clone()]);
                    y_trimmed.extend(&y[range.clone()]);

                    let n = range.len();
                    range.start = 0 + offset;
                    range.end = n + offset;
                    offset += n;

                    range
                })
                .collect();

            let mut x = x_trimmed;
            let mut y = y_trimmed;

            //---
            // Extend terms and bounds with baseline data
            //---

            // First baseline
            let basis = experiment.baseline().basis(&x, &range_values);

            let baseline_terms = experiment.baseline().terms(&range_values);
            let nb = baseline_terms.len()/2;

            let baseline_lower_bounds = vec![
                experiment.baseline().bounds().lower; nb*2
            ];
            let baseline_upper_bounds = vec![
                experiment.baseline().bounds().upper; nb*2
            ];
            
            terms.extend(baseline_terms);
            lower_bounds.extend(baseline_lower_bounds);
            upper_bounds.extend(baseline_upper_bounds);

            //---
            // Scale x 0-1 and unfold complex data
            //---
            
            let mut xmin = f64::INFINITY;
            let mut xmax = f64::NEG_INFINITY;
            for x in x.iter() {
                xmin = xmin.min(*x);
                xmax = xmax.max(*x);
            }

            let xcoeff = xmax - xmin;

            let f1 = |x: f64| (x-xmin)/xcoeff;
            let g1 = |x: f64| x*xcoeff + xmin;

            let f2 = |x: f64| x/xcoeff;
            let g2 = |x: f64| x*xcoeff;

            // First, scale the x data
            for x in x.iter_mut() {
                *x = f1(*x);
            }

            // Then terms... starting with position (and bounds)
            let iterator = lower_bounds.iter_mut()
                .zip(upper_bounds.iter_mut())
                .zip(terms.iter_mut().take(nl));

            for ((l, u), x) in iterator.step_by(4) {
                *l = f1(*l);
                *u = f1(*u);
                *x = f1(*x);
            }

            // Followed by width
            let iterator = lower_bounds.iter_mut()
                .zip(upper_bounds.iter_mut())
                .zip(terms.iter_mut().take(nl));

            for ((l, u), x) in iterator.skip(1).step_by(4) {
                *l = f2(*l);
                *u = f2(*u);
                *x = f2(*x);
            }

            // Unfold the y data
            let y: Vec<_> = y.into_iter()
                .map(|y| vec![y.re, y.im])
                .flatten()
                .collect();

            //---
            // Add phase, scaling as necessary
            //---

            let np = *experiment.phase().np();
            let phase_terms = experiment.phase().terms().clone();
            let mut phase_terms = phase_terms[0..np].to_vec();

            if np > 1 {
                phase_terms[0] = phase_terms[0] + phase_terms[1]*xmin;
                phase_terms[1] = g2(phase_terms[1]);
            }

            terms.extend(&phase_terms[..np]);

            // Lower and upper bounds are set at extremes of +/- 
            // (with first order term not allowed to contribute more than
            // +/- pi over full length of x data)
            let [xmin_full, xmax_full] = experiment.domain();
            let xcoeff_full = xmax_full - xmin_full;

            let phase_bounds = *experiment.phase().bounds();

            let order_0 = phase_bounds/xcoeff_full*xcoeff;
            let order_1 = PI/xcoeff_full*xcoeff;

            let lower = phase_terms[0] - order_0;
            let phase_lower_bounds = vec![lower, -order_1];
            lower_bounds.extend(&phase_lower_bounds[..np]);

            let upper = phase_terms[0] + order_0;
            let phase_upper_bounds = vec![upper, order_1];
            upper_bounds.extend(&phase_upper_bounds[..np]);

            //---
            // Generating fit object
            //---

            // Initializing nlopt object
            let n = x.len();
            let x = Array::from_shape_vec(( n, ), x).unwrap();
            let y = Array::from_shape_vec((n, 2), y).unwrap();

            let fit = fitting::Fit1D::new(
                x, y, ranges, assignments, basis, nl, nb, np, tol
            );

            type F = fn(
                    &[f64], Option<&mut [f64]>, &mut fitting::Fit1D
                ) -> f64;

            let mut opt: Nlopt<F, fitting::Fit1D>;

            opt = Nlopt::new(
                Algorithm::Slsqp, 
                //Algorithm::Cobyla, 
                nl + nb*2 + np, 
                fitting::Fit1D::obj, 
                Target::Minimize, 
                fit
            );

            opt.set_lower_bounds(&lower_bounds[..])
                .expect("Unexpected error setting lower bounds.");
            opt.set_upper_bounds(&upper_bounds[..])
                .expect("Unexpected error setting upper bounds.");

            // Updating constraints
            let mut constraints = experiment.constraints(
                &self.species, &self.ppm
            );

            for c in constraints.iter_mut() {
                c.update(&self.constraint_leeways);
            }

            // Equality constraints
            let parameter_equalities = PeakEquality::new(
                &constraints, nl + nb*2 + np, xcoeff
            );
            let m = parameter_equalities.m();
            let tolerance: Vec<f64> = iter::repeat(1e-4).take(m).collect();

            opt.add_equality_mconstraint(
                m, 
                PeakEquality::mobj, 
                parameter_equalities, 
                &tolerance
            ).expect("Error formulating peak constraints.");

            // Inequality constraints
            let parameter_inequalities = PeakInequality::new(
                &constraints, nl + nb*2 + np, xcoeff
            );
            let m = parameter_inequalities.m();
            let tolerance: Vec<f64> = iter::repeat(1e-4).take(m).collect();

            opt.add_inequality_mconstraint(
                m, 
                PeakInequality::mobj, 
                parameter_inequalities, 
                &tolerance
            ).expect("Error formulating peak constraints.");

            // Phase constraints are needed for 1st order correction
            if np > 1 {

                let phase_inequalities = PhaseInequality::new(
                    phase_bounds.clone(), nl+nb*2+np
                );
                let tolerance: Vec<f64> = iter::repeat(1e-4).take(4).collect();
                
                opt.add_inequality_mconstraint(
                    4, 
                    PhaseInequality::mobj, 
                    phase_inequalities, 
                    &tolerance
                ).expect("Error formulating phase constraints.");
            }

            // Basic constraints
            opt.set_maxtime(10.0).unwrap();
            opt.set_xtol_rel(tol).unwrap();

            //---
            // Run the fit
            //---

            // Run the optimization
            let _out = opt.optimize(&mut terms[..]);

            println!("{:?}", _out);

            //---
            // Update values
            //---
            
            // Phase requires a nonlinear transformation.
            // Phase was calculated at xmin set to 0 and xmax set to 1.
            // This needs to be remapped to a + bx format.
            // However... xmin is on the right and xmax on the left.
            // So when the slope is applied, it is actually at maximum on the left.
            
            let mut phase_right = terms[nl+nb*2];
            let mut phase_left = phase_right;

            if np > 1 {
                phase_left = phase_right + terms[nl+nb*2+1];
            }

            let slope = (phase_left - phase_right)/xcoeff;
            let intercept = phase_right - slope*xmin;

            experiment.phase_mut().set_terms(
                &[intercept, slope]
            );

            // Update the baseline
            experiment.baseline_mut().set_terms(
                &terms[nl..(nl + nb*2)], &range_values
            );

            // Snip off the baseline/phase
            let mut terms: Vec<_> = terms.into_iter().take(nl).collect();

            // Rescaling from 0 -> 1
            for x in terms.iter_mut().step_by(4) {
                *x = g1(*x);
            }

            for x in terms.iter_mut().skip(1).step_by(4) {
                *x = g2(*x);
            }

            experiment.replace(&self.species, &terms, &self.ppm);
        }

        Ok( () )
    }


    //-------------------------------------------------------------------------
    // Split chemical shift into ranges
    fn calculate_manual_ranges (
            x: &[f64], lower_bounds: &[f64], upper_bounds: &[f64], cutoff: f64
        ) -> Vec<Range<usize>> {

        // Quick order check
        let x_min = *x.last().expect("Missing data.");
        let x_max = *x.first().expect("Missing data.");
        assert!(x_min < x_max, "Unexpected sort order.");

        // Common work arrays
        let mut p_temp = vec![0.0f64; 4];
        let mut y_temp = vec![0.0f64; x.len()*2];

        // Calculate lower and upper range of each peak based on bounds
        let mut ranges = lower_bounds.iter()
            .chunks(4)
            .into_iter()
            .zip(upper_bounds.iter().chunks(4).into_iter())
            .map(|(mut lb, mut ub)| {

                // We need both low and high positions
                let p_low = lb.next().unwrap();
                let p_high = ub.next().unwrap();
                p_temp[0] = *p_high;

                // Height and width come from upper bounds
                let _ = lb.next().unwrap();
                let h = ub.next().unwrap();
                p_temp[1] = *h;

                let _ = lb.next().unwrap();
                let w = ub.next().unwrap();
                p_temp[2] = *w;

                // Fraction comes from lower bounds
                let f = lb.next().unwrap();
                let _ = ub.next().unwrap();
                p_temp[3] = *f; 

                fitting::f_peak(&p_temp, x, &mut y_temp, 1e-4);

                // TODO: make the cutoff an option
                let y_cutoff = h*cutoff;

                let mut range = x.iter()
                    .zip(y_temp.iter().step_by(2))
                    .skip_while(|(_, y)| **y < y_cutoff)
                    .take_while(|(_, y)| **y >= y_cutoff)
                    .fold((x_min, x_max), | a, (x, _)| {
                        (a.0.max(*x), a.1.min(*x), )
                    });
                
                // The generated range is for the upper
                // position, so it must be extended for the lower
                // position.
                range.1 -= p_high - p_low;

                range
            })
            .sorted_by(|a, b| a.0.total_cmp(&b.0))
            .rev();

        // Popping off first range for aggregation
        let first = ranges.next().expect("Error identifying ranges");

        // Combining ranges (relying on sort order)
        // If we fall into last range, extend it; otherwise, push a new entry.
        ranges.into_iter()
            .fold(vec![first], |mut a, r| {

                let n = a.len();

                if r.0 > a[n - 1].1 {
                    a[n - 1].1 = a[n - 1].1.min(r.1);
                } else {
                    a.push(r);
                }

                a
            })
            // Then onvert start, stop indexes into formal ranges
            .into_iter()
            .map(|r| {
                let i = x.iter()
                    .enumerate()
                    .find(|(i, x)| **x <= r.0)
                    .map(|(i, _)| i)
                    .expect("Error calculating ranges");

                let j = x.iter()
                    .enumerate()
                    .find(|(i, x)| **x <= r.1)
                    .map(|(i, _)| i)
                    .expect("Error calculating ranges");

                i..(j+1)
            })
            .collect()

    }
}



//=============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct Fit1DTemplate {

    fit: String,

    #[serde(default)]
    species: Globs,
    #[serde(default)] 
    experiments: Globs,

    #[serde(default = "default_ppm")]
    ppm: [f64; 2],

    #[serde(default)] 
    constraint_leeways: ConstraintLeeways,

    #[serde(default = "default_param")]
    width: f64,
    #[serde(default = "default_param")]
    height: f64,
    #[serde(default = "default_param")]
    fraction: f64,

    #[serde(default)] 
    general_bounds: Bounds,
    #[serde(default)] 
    offset_bounds: Bounds,

    #[serde(default)] 
    outputs: Vec<Output>,
}

//-----------------------------------------------------------------------------

fn default_ppm() -> [f64; 2] {
    [f64::NEG_INFINITY, f64::INFINITY]
}

fn default_param() -> f64 {
    f64::NAN
}
