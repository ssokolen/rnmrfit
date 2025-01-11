use anyhow::{Context, Result, anyhow, bail};
use glob::Pattern;
use log::error;
use serde::{Deserialize};
use std::{
    fs,
    ops::Add,
    path::PathBuf,
};
use thiserror::Error;
use yaml_rust::{Yaml, YamlLoader, YamlEmitter};

mod baseline;
mod area;
mod constraints;
mod data;
mod experiment;
mod fit;
mod job;
mod parameters;
mod peaks;
mod phase;
mod plot;
mod resonance;
mod species;

pub use baseline::{Baseline1D, Baseline1DTemplate, DEFAULT_BASELINE};
pub use area::{AreaCalculation1D, AreaCalculation1DTemplate};
pub use constraints::{
    Constraint, ConstraintLeeways, PeakEquality, PeakInequality, PhaseInequality
};
pub use data::{Data1D, Data1DTemplate, Processing1DTemplate};
pub use experiment::{Experiment1D};
pub use fit::{Fit1D, Fit1DTemplate, Output};
pub use job::{Job1D, Ranges};
pub use parameters::{ParameterList1D, ParameterList1DTemplate};
pub use peaks::{PeaksTemplate, Peaks};
pub use phase::{Phase1D, Phase1DTemplate, DEFAULT_PHASE};
pub use plot::{Plot1D, Plot1DTemplate};
pub use resonance::{Resonance, ResonanceTemplate};
pub use species::{Species1D, Species1DTemplate};


// Unbounded ppm range
pub const UNBOUNDED: [f64; 2] = [f64::NEG_INFINITY, f64::INFINITY];


//=============================================================================
// Common errors

//-----------------------------------------------------------------------------
#[derive(Error, Debug)]
pub enum ParsingError {
    #[error("Entry does not match any recognized pattern: \n\n{0}\n\n")]
    NoEntry(String),

    #[error("Error parsing {0} entry. {1}. \n\n{2}\n\n")]
    MalformedEntry(String, anyhow::Error, String),

    #[error("Cannot define default settings after any other entry: \n\n{0}\n.")]
    SettingsOrder(String),

    #[error("Cannot define resonance, species, or data after fit, plot, calculate, or list operations: \n\n{0}\n.")]
    SpectraOrder(String),
}


// The different functions help identify where the error came from
impl ParsingError {

    pub fn no_entry(yaml: String) -> Self { 

        let error = ParsingError::NoEntry(yaml);
        error!("{}", &error);

        error
    }

    pub fn malformed_entry(
            entry: &str, error: anyhow::Error, yaml: String
        ) -> Self { 

        let error = ParsingError::MalformedEntry(
            entry.to_string(), error, yaml
        );
        error!("{}", &error);

        error
    }

    pub fn settings_order(yaml: String) -> Self { 

        let error = ParsingError::SettingsOrder(yaml);
        error!("{}", &error);

        error
    }

    pub fn spectra_order(yaml: String) -> Self { 

        let error = ParsingError::SpectraOrder(yaml);
        error!("{}", &error);

        error
    }
}


//-----------------------------------------------------------------------------
#[derive(Error, Debug)]
pub enum SetupError {
    #[error("Error in YAML setup.")]
    IncorrectEntry(#[from] anyhow::Error),

    #[error("No peaks within ppm bounds {0}-{1}.")]
    NoPeaks(f64, f64),
}

// Unlike ParsingError, all errors come from Fit, so only function needed
impl SetupError {

    /// Logging is essential for users not using a terminal
    pub fn error(error: anyhow::Error, name: &String) -> Self { 

        let error = anyhow!("Fit \"{}\" -- {}", name, error);

        error!("{:?}", &error);
        SetupError::IncorrectEntry(error)
    }
}



//=============================================================================
// Helper functions

//-----------------------------------------------------------------------------
pub fn path_to_yaml_docs(path: &PathBuf) -> Result<Vec<Yaml>> {

    let path_string = path.to_string_lossy();

    let contents = fs::read_to_string(path)
        .with_context(|| format!("\"{}\" could not be read.", &path_string) )?;

    YamlLoader::load_from_str(&contents[..])
        .with_context(|| format!("\"{}\" could not be parsed.", &path_string))
}


//-----------------------------------------------------------------------------
pub fn yaml_to_string(yaml: &Yaml) -> String {

    let mut doc_string = String::new();
    let mut emitter = YamlEmitter::new(&mut doc_string);
    emitter.dump(yaml).expect("Unexpected error writing YAML");
    doc_string
}


//-----------------------------------------------------------------------------
pub fn sum_across<T: Add<Output = T> + Copy>(mut x: Vec<Vec<T>>) -> Vec<T> {
    
    if x.len() == 0 { 
        return vec![] 
    }

    let init = x.pop().unwrap();
    if x.len() == 0 { 
        return init 
    }

    (0..x[0].len())
        .map(|i| {
            x.iter().fold(init[i], |acc, line| acc + line[i])
        })
        .collect()
}


//=============================================================================
// Helper structs

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum Globs {
  One(String),
  Multiple(Vec<String>),
}


impl Default for Globs {

    fn default() -> Globs {
        Globs::One("*".to_string())
    }
}


impl Globs {

    pub fn patterns(self) -> Result<Vec<Pattern>> {

        let strings = match self {
            Globs::One(string) => {
                if string == "" { vec![] }
                else { vec![string] }
            },
            Globs::Multiple(strings) => strings
        };

        strings.iter()
            .map(|x| {
                Pattern::new(x).map_err(|_| {
                    anyhow!(
                        "\"{}\" is not a valid \"glob\" pattern.", x 
                    )
                })
            })
            .collect()
    }
}


//=============================================================================
// Common components


//-----------------------------------------------------------------------------
#[derive(Clone, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Scaffold {
    Peak,
    Resonance,
    Species,
}


impl Default for Scaffold {

    fn default() -> Self {
        Scaffold::Species
    }
}


//-----------------------------------------------------------------------------
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameters {

    #[serde(default = "default_nan")]
    pub position: f64,
    #[serde(default = "default_nan")]
    pub width: f64,
    #[serde(default = "default_nan")]
    pub height: f64,
    #[serde(default = "default_nan")]
    pub fraction: f64,

}


fn default_nan() -> f64 {
    f64::NAN
}


impl Default for Parameters {

    fn default() -> Self {
        Parameters {
            position: f64::NAN,
            width: f64::NAN,
            height: f64::NAN,
            fraction: f64::NAN,
        }
    }
}


// Sensible global defaults
pub const DEFAULT_PARAMETERS: Parameters = Parameters {
    
    // absolute ppm position cannot be guessed a priori
    position: f64::NAN,

    // 1 Hz is a reasonable starting width
    width: 1.0,

    // default heights are appied last
    height: f64::NAN,

    // 0.2 is a reasonable value as most peaks will be closer to Lorentz
    fraction: 0.2
};


impl Parameters {

    pub fn from_position(position: f64) -> Parameters {

        Parameters {
            position: position,
            ..Default::default()
        }
    }


    pub fn update(&mut self, parameters: &Parameters) {

        if self.position.is_nan() {
            self.position = parameters.position
        }

        if self.width.is_nan() {
            self.width = parameters.width
        }

        if self.height.is_nan() {
            self.height = parameters.height
        }

        if self.fraction.is_nan() {
            self.fraction = parameters.fraction
        }
    }


    pub fn terms(&self, frequency: f64) -> Vec<f64> {
        vec![self.position, self.width/frequency, self.height, self.fraction]
    }


    pub fn update_terms(&mut self, terms: &[f64]) -> Result<()> {

        if terms.len() != 4 {
            bail!("Invalid number of terms.")
        }

        self.position = terms[0];
        self.width = terms[1];
        self.height = terms[2];
        self.fraction = terms[3];

        Ok( () )
    }
}


//-----------------------------------------------------------------------------
#[derive(Clone, Copy, Debug, Default, Deserialize)]
#[serde(default)] 
pub struct Bounds {

    pub position: Bound,
    pub width: Bound,
    pub height: Bound,
    pub fraction: Bound,

}

// Sensible global defaults
pub const DEFAULT_GENERAL_BOUNDS: Bounds = Bounds {
    
    // absolute ppm position cannot be reasonably bounded
    position: Bound { lower: f64::NEG_INFINITY, upper: f64::INFINITY },

    // reasonable peaks should be wider than 0.05 Hz and less than 10 Hz
    // (although even this upper bound should really be more restrictive)
    width: Bound { lower: 0.05, upper: 10.0 }, 

    // since intensity is scaled to a maximum of one for the tallest peak,
    // no peak should technically exceed one, but some leeway is provided
    // for baseline impact
    height: Bound { lower: 0.0, upper: 2.0 },

    // prevent fraction from hitting 1.0
    fraction: Bound { lower: 0.0, upper: 1.0 - 1e-5 },
};


pub const DEFAULT_OFFSET_BOUNDS: Bounds = Bounds {
    
    // Only position can be offset bounded
    position: Bound { lower: -0.2, upper: 0.2 },

    // All the rest are left as-is
    width: Bound { lower: f64::NEG_INFINITY, upper: f64::INFINITY },
    height: Bound { lower: f64::NEG_INFINITY, upper: f64::INFINITY },
    fraction: Bound { lower: f64::NEG_INFINITY, upper: f64::INFINITY },
};


impl Bounds {


    //-------------------------------------------------------------------------
    pub fn offset(&self, parameters: &Parameters) -> Bounds {

        Bounds {
            position: self.position.offset(parameters.position),
            width: self.width.offset(parameters.width),
            height: self.height.offset(parameters.height),
            fraction: self.fraction.offset(parameters.fraction),
        }
    }


    //-------------------------------------------------------------------------
    pub fn update(&mut self, bounds: &Bounds) {

        self.position.update(&bounds.position);
        self.width.update(&bounds.width);
        self.height.update(&bounds.height);
        self.fraction.update(&bounds.fraction);

    }


    //-------------------------------------------------------------------------
    pub fn bounds(&self) -> [[f64; 2]; 4] {
        
        [
            self.position.bounds(),
            self.width.bounds(),
            self.height.bounds(),
            self.fraction.bounds(),
        ]
    }
}


//-----------------------------------------------------------------------------
#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(default)] 
pub struct Bound {

    lower: f64,
    upper: f64,

}


impl Default for Bound {

    fn default() -> Self {
        Bound {
            lower: f64::NEG_INFINITY,
            upper: f64::INFINITY,
        }
    }
}


impl Bound {

    //-------------------------------------------------------------------------
    pub fn new(lower: f64, upper: f64) -> Bound {
        
        Bound {
            lower: lower,
            upper: upper
        }
    }   


    //-------------------------------------------------------------------------
    pub fn offset(&self, parameter: f64) -> Bound {

        Bound {
            lower: parameter + self.lower,
            upper: parameter + self.upper,
        }
    }


    //-------------------------------------------------------------------------
    pub fn update(&mut self, bound: &Bound) {

        self.lower = self.lower.max(bound.lower);
        self.upper = self.upper.min(bound.upper);

    }


    //-------------------------------------------------------------------------
    pub fn bounds(&self) -> [f64; 2] {

        let mut lower = self.lower;

        if lower.is_nan() {
            lower = f64::NEG_INFINITY;
        }

        let mut upper = self.upper;

        if upper.is_nan() {
            upper = f64::INFINITY;
        }

        [lower, upper]
    }
}

