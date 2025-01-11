use anyhow::{Result, bail};
use cascade::cascade;
use getset::{Getters, Setters};
//use ndarray::prelude::*;
//use rustfft::num_complex::Complex;
use serde::Deserialize;
use std::{
    f64::consts::PI,
    iter,
    convert::TryInto
};

use crate::processing::Bound;


//=============================================================================
// Template

#[derive(Clone, Debug, Default, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct Phase1DTemplate {

    #[serde(default)]
    order: Option<usize>,
    #[serde(default)]
    bounds: Option<f64>,
}


impl Phase1DTemplate {

    pub fn update(&mut self, template: &Phase1DTemplate) {

        if self.order.is_none() {
            self.order = *template.order();
        }

        if self.bounds.is_none() {
            self.bounds = *template.bounds();
        }
    }

}


// Sensible global defaults
pub const DEFAULT_PHASE: Phase1DTemplate = Phase1DTemplate {
    
    order: Some(1),
    bounds: Some(360.0),
};


//=============================================================================
// Phase1D
//
// Note, terms are stored rather strangely as theta at ppm[0] and additive
// theta at ppm[1], i.e., (theta[0], theta[1] - theta[0]). This works because 
// each experiment does a regression on calculated angles anyway.

#[derive(Clone, Debug, Getters, Setters)]
pub struct Phase1D {

    #[getset(get = "pub")]
    terms: Vec<f64>,

    #[getset(get = "pub", set = "pub")]
    bounds: f64,

    #[getset(get = "pub")]
    np: usize,
}


//=============================================================================
// Constructors

impl Phase1D {

    //-------------------------------------------------------------------------
    pub fn new(order: usize) -> Result<Phase1D> { 

        if order > 2 {
            bail!("Phase order may only be 0 or 1.")
        }

        let terms: Vec<f64> = iter::repeat(0.0)
            .take(order + 1)
            .collect();

        let phase = Phase1D {
            terms: terms,
            bounds: PI/4.0,
            np: order + 1,
        };

        Ok(phase)
    }
}


//=============================================================================
// Methods

impl Phase1D {

    //-------------------------------------------------------------------------
    pub fn set_terms(&mut self, terms: &[f64]) {

        if terms.len() != 2 {
            panic!("Invalid number of terms.")
        }

        self.terms = terms.to_vec();
    }
}
