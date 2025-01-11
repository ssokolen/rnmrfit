use anyhow::{Context, Result, anyhow, bail};
use cascade::cascade;
use getset::{Getters, Setters};
use itertools::Itertools;
use rustfft::num_complex::Complex;
use serde::Deserialize;
use std::{
    iter,
    collections::HashMap,
};

use crate::processing;
use crate::processing::{
    Bounds, Constraint, ConstraintLeeways, Parameters, 
    Job1D, PeaksTemplate, Peaks, Resonance,
};

use crate::fitting::{f_peak, f_peak_area};

// Unbounded ppm range
const UNBOUNDED: [f64; 2] = [f64::NEG_INFINITY, f64::INFINITY];


//=============================================================================
// Species1D

#[derive(Clone, Debug, Default, Getters, Setters)]
#[getset(get = "pub")]
pub struct Species1D {

    #[getset(set = "pub")]
    name: String,

    #[getset(skip)]
    resonances: Vec<Resonance>,
    proportions: Vec<f64>,

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

impl Species1D {


    //-------------------------------------------------------------------------
    pub fn from_template(
        template: Species1DTemplate,
        job: &Job1D
    ) -> Result<Self> {

        // Extract out peaks, names, proportions
        let (peaks, names, proportions) = match template.resonances {
            ResonanceTemplate::Ok(p, n, a) => (p, n, a),
            ResonanceTemplate::Err(message) => bail!(message)
        };

        // Peaks are converted first
        let peaks: Vec<Result<Resonance>> = peaks
            .iter()
            .map(|x| {
                match x {
                    PeaksTemplate::Label(ref label) => {
                        job.resonances().get(label)
                            .ok_or(
                                anyhow!(
                                    "Resonance \"{}\" is not defined", 
                                    label 
                                )
                            )
                            .and_then(|x| Ok(x.clone()))
                    }
                    _ => {
                        Peaks::from_template(x.clone())
                            .and_then(|x| Ok( Resonance::from_peaks(x) ))
                    }
                }
                .context("Invalid peak definition")
            })
            .collect();

        let mut resonances: Vec<Resonance> = peaks
            .into_iter()
            .collect::<Result<Vec<Resonance>>>()?;

        // Proportion are handled later...

        // Rename resonances if applicable
        if names.len() > 0 {
            for (r, n) in resonances.iter_mut().zip(names.into_iter()) {
                r.set_name(n);
            }
        }

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

        let species = cascade!{
            Species1D::from_resonances(resonances);
                ..set_name(template.species);
                ..set_width(width);
                ..set_height(height);
                ..set_fraction(fraction);
                ..set_proportions(proportions);
                ..set_constraint_leeways(template.constraint_leeways);
                ..set_general_bounds(template.general_bounds);
                ..set_offset_bounds(template.offset_bounds);
        };

        Ok( species )
    }
}


//=============================================================================
// Constructors

impl Species1D {

    //-------------------------------------------------------------------------
    pub fn from_resonances(resonances: Vec<Resonance>) -> Species1D {

        Species1D {
            resonances: resonances,
            ..Default::default()
        }
    }
}


//=============================================================================
// Simple getters/setters

impl Species1D {

    //-------------------------------------------------------------------------
    pub fn resonances<'a>(
        &'a self,
        frequency: f64,
        ppm: &'a[f64; 2]
    ) -> impl Iterator<Item = &'a Resonance> {
        
        self.resonances.iter()
            .filter(move |x| x.is_in(frequency, ppm))

    }


    //-------------------------------------------------------------------------
    pub fn resonances_mut<'a>(
        &'a mut self,
        frequency: f64,
        ppm: &'a[f64; 2]
    ) -> impl Iterator<Item = &'a mut Resonance> {
        
        self.resonances.iter_mut()
            .filter(move |x| x.is_in(frequency, ppm))

    }


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

    pub fn set_proportions(&mut self, proportions: Vec<f64>) {
        let n_proportions = proportions.len();
        let n_resonances = self.resonances.len();
        
        if n_proportions == 0 {
            return
        } else if n_proportions == 1 {
            self.proportions = iter::repeat(proportions[0])
                .take(n_resonances)
                .collect();
        } else if n_proportions != n_resonances {
            panic!("Number of proportions must match number of resonances.");
        } else {
            self.proportions = proportions;
        }
    }
}


//==============================================================================
// Methods

impl Species1D {

    //-------------------------------------------------------------------------
    pub fn n_peaks(&self, frequency: f64, ppm: &[f64; 2]) -> usize {

        self.resonances(frequency, ppm)
            .map(|x| x.n_peaks())
            .sum()
    }


    //-------------------------------------------------------------------------
    pub fn terms(&self, frequency: f64, ppm: &[f64; 2]) -> Vec<f64> {

        self.resonances(frequency, ppm)
            .map(|resonance| {
                let mut resonance = resonance.clone();
                resonance.update(frequency, self.parameters());
                resonance.terms(frequency)
            })
            .flatten()
            .collect()
    }

    //-------------------------------------------------------------------------
    pub fn replace(
        &mut self, 
        terms: &[f64], 
        frequency: f64, 
        ppm: &[f64; 2]
    ) {

        if terms.len()/4 != self.n_peaks(frequency, ppm) {
            panic!("Invalid number of terms.")
        }

        let mut i = 0;

        for r in self.resonances_mut(frequency, ppm) {
            let j = i + r.n_peaks()*4;
            r.replace(&terms[i..j], frequency);
            i = j;
        }
    }


    //-------------------------------------------------------------------------
    pub fn bounds(
        &self, 
        frequency: f64,
        ppm: &[f64; 2]
    ) -> Vec<[f64; 2]> {

        let mut offset_bounds = self.offset_bounds().bounds();
        offset_bounds[1][0] /= frequency;
        offset_bounds[1][1] /= frequency;
        let offset_bounds = IntoIterator::into_iter(offset_bounds);

        let mut general_bounds = self.general_bounds().bounds();
        general_bounds[1][0] /= frequency;
        general_bounds[1][1] /= frequency; 
        let general_bounds = IntoIterator::into_iter(general_bounds);

        // Generate offset bounds then update with general ones
        let bounds: Vec<[f64; 2]> = self.terms(frequency, ppm)
            .into_iter()
            .zip(offset_bounds.cycle())
            .map(|(x, [l, u])| [x+l, x+u])
            .zip(general_bounds.cycle())
            .map(|([l1, u1], [l2, u2])| {
                [f64::max(l1, l2), f64::min(u1, u2)]
            })
            .collect();

        // Update bounds with resonance information
        bounds
            .into_iter()
            .zip(
                self.resonances(frequency, ppm)
                    .map(|x| x.bounds(frequency))
                    .flatten()
            )
            .map(|([l1, u1], [l2, u2])| {
                [f64::max(l1, l2), f64::min(u1, u2)]
            })
            .collect()
    }


   //--------------------------------------------------------------------------
   pub fn constraints(
        &self, 
        frequency: f64,
        ppm: &[f64; 2],
        mut offset: usize,
    ) -> Vec<Constraint> {

        // New area constraints have to be built up from lhs/rhs indexes
        let mut indexes: Vec<Vec<usize>> = Vec::new();
        let mut proportions: Vec<f64> = Vec::new();

        // Updated resonance constraints can be added directly
        let mut constraints: Vec<Constraint> = Vec::new();

        // Pull out resonance constraints and indexes
        for (i, resonance) in self.resonances(frequency, ppm).enumerate() {
            
            let n_peaks = resonance.n_peaks();

            constraints.extend(resonance.constraints(frequency, offset));

            // If there is a proportional area, generate indexes
            if let Some(proportion) = self.proportions.get(i) {
                indexes.push( (offset .. (offset + n_peaks)).collect() );
                proportions.push(*proportion);
            }

            offset += n_peaks;
        }

        // Update resonance constraints
        for constraint in constraints.iter_mut() {
            constraint.update(&self.constraint_leeways);
        }

        // Then add area constraints
        if proportions.len() > 0 {
            let leeways = &self.constraint_leeways;

            let a0 = proportions[0];
            let j0 = &indexes[0];

            for i in 1..proportions.len() {
                let a = proportions[i]/a0;
                let j1 = &indexes[i];

                constraints.push(
                    Constraint::area((j0.clone(), j1.clone()), a, leeways.area)
                );
            }
        }

        constraints
    }

    //-------------------------------------------------------------------------
    pub fn peak_names(&self) -> Vec<String> {

        self.resonances
            .iter()
            .map(|x| x.peak_names().clone())
            .flatten()
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn resonance_names(&self) -> Vec<String> {

        self.resonances
            .iter()
            .map(|x| x.name().clone() )
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn peak_lineshapes(
        &self, 
        x: &[f64], 
        frequency: f64
    ) -> Vec<Vec<Complex<f64>>> {

        let (xmin, xmax) = x.iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |acc, y| {
                (acc.0.min(*y), acc.1.max(*y))
            });

        self.terms(frequency, &[xmin, xmax])
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

        self.terms(frequency, &UNBOUNDED)
            .chunks(4)
            .map(|y| f_peak_area(&y[..], 1e-10) )
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn resonance_lineshapes(
        &self, 
        x: &[f64], 
        frequency: f64
    ) -> Vec<Vec<Complex<f64>>> {

        self.resonances
            .iter()
            .map(|y| y.resonance_lineshape(x, frequency) )
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn resonance_areas(&self,  frequency: f64) -> Vec<f64> {

        self.resonances
            .iter()
            .map(|y| y.resonance_area(frequency) )
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn species_lineshape(
        &self, 
        x: &[f64], 
        frequency: f64
    ) -> Vec<Complex<f64>> {

        let lineshapes = self.peak_lineshapes(x, frequency);
        processing::sum_across(lineshapes)
    }


    //-------------------------------------------------------------------------
    pub fn species_area(&self,  frequency: f64) -> f64 {

        self.peak_areas(frequency).iter().sum()
    }
}


//=============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct Species1DTemplate {

    species: String,
    resonances: ResonanceTemplate,

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


// Collections of peaks with names and proportions
#[derive(Clone, Debug, Deserialize)]
#[serde(from = "ResonanceTemp")] 
enum ResonanceTemplate {
    Ok(Vec<PeaksTemplate>, Vec<String>, Vec<f64>),
    Err(String)
}


impl From<ResonanceTemp> for ResonanceTemplate {

    fn from(temp: ResonanceTemp) -> Self {

        // TODO: finish patterns
        let resonance = match temp {
            /*
            NamesProportion(HashMap<String, (String, f64)>),
            ProportionNames(HashMap<String, (f64, String)>),
            Proportions(HashMap<String, f64>),
            */
            ResonanceTemp::Names(map) => {
                let (templates, names): (Vec<_>, Vec<_>) = map.into_iter()
                    .unzip();

                let templates: Result<Vec<_>> = templates.into_iter()
                    .map(|x| PeaksTemplate::from_string(x))
                    .collect();

                match templates {
                    Ok(templates) => {
                        ResonanceTemplate::Ok(
                            templates, names, vec![]
                        )
                    },
                    Err(message) => {
                        ResonanceTemplate::Err(message.to_string())
                    }
                }
            },

            ResonanceTemp::Multiple(templates) => {
                let names: Vec<String> = templates.iter()
                    .map(|x| x.to_string())
                    .collect();
                ResonanceTemplate::Ok(
                    templates, names, vec![]
                )
            },
            
            ResonanceTemp::Single(template) => {
                let name = template.to_string();
                ResonanceTemplate::Ok(
                    vec![template], vec![name], vec![]
                )
            }
            _ => ResonanceTemplate::Err("Cannot parse resonances".to_string())
        };

        resonance
    }
}


#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)] 
enum ResonanceTemp {
    NamesProportion(HashMap<String, (String, f64)>),
    ProportionNames(HashMap<String, (f64, String)>),
    Proportions(HashMap<String, f64>),
    Names(HashMap<String, String>),
    Multiple(Vec<PeaksTemplate>),
    Single(PeaksTemplate),
}




