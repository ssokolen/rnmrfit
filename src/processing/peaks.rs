use anyhow::{Context, Result, anyhow, bail};
use lazy_static::lazy_static;
use regex::Regex;
use serde::{Deserialize};
use statrs::function::factorial;
use std::{
    iter,
    convert::{TryFrom, TryInto},
};

use crate::processing::{Parameters};

//=============================================================================
// Peaks

#[derive(Clone, Debug)]
pub enum Peaks {
    Pattern(String),
    Parameters(Vec<Parameters>, Vec<f64>),
}


//-----------------------------------------------------------------------------
impl Default for Peaks {

    fn default() -> Peaks {
        Peaks::Parameters(vec![Parameters::from_position(0.0)], vec![1.0])
    }
}


//=============================================================================
// YAML parsing

//-----------------------------------------------------------------------------
impl Peaks {

    //-------------------------------------------------------------------------
    pub fn from_template(peaks: PeaksTemplate) -> Result<Self> {

        match peaks {
            PeaksTemplate::Label(_) => {
                bail!("Peaks cannot be generated via reference label.")
            },
            PeaksTemplate::Pattern(text) => {
                Peaks::from_pattern(text.clone())
            },
            PeaksTemplate::Positions(positions) => {
                Peaks::from_positions(positions.clone(), vec![])
            },
            PeaksTemplate::Areas((positions, areas)) => {
                Peaks::from_positions(positions.clone(), areas.clone())
            }
        }
    }
}


//=============================================================================
// Constructors

//-----------------------------------------------------------------------------
impl Peaks {

    //-------------------------------------------------------------------------
    pub fn from_pattern(pattern: String) -> Result<Self> {

        // Check to make sure that the pattern is valid by converting to areas
        let peaks = Peaks::Pattern(pattern.clone());
        peaks.expand(600.0)
            .with_context(|| {
                format!("Could not generate Peak from {}", &pattern)
            })?;

        Ok( Peaks::Pattern(pattern) )
    }


    //-------------------------------------------------------------------------
    pub fn from_positions(
        positions: Vec<f64>, 
        mut areas: Vec<f64>
    ) -> Result<Self> {

        if areas.len() == 0 {
            areas = iter::repeat(1.0).take(positions.len()).collect();
        } else if positions.len() != areas.len() {
            bail!("Peak position and areas must be same length.");
        }

        let parameters: Vec<Parameters> = positions
            .into_iter()
            .map(|x| Parameters::from_position(x))
            .collect();

        Ok( Peaks::Parameters(parameters, areas) )
    }
}


//=============================================================================
// Methods

//-----------------------------------------------------------------------------
impl Peaks {

    //--------------------------------------------------------------------------
    pub fn n_peaks(&self) -> usize {

        self.extract(0.0)
            .expect("Coupling could not be converted into peak parameters.")
            .0
            .len()
    }


    //-------------------------------------------------------------------------
    pub fn terms(&self, frequency: f64) -> Vec<f64> {

        
        let (parameters, areas) = self.extract(frequency)
            .expect("Coupling could not be converted into peak parameters.");

        parameters
            .into_iter()
            .zip(areas.into_iter())
            .map(|(p, a)| {

                let mut terms = p.terms(frequency);
                terms[2] *= a;

                terms
            })
            .flatten()
            .collect()
    }

    //-------------------------------------------------------------------------
    // Note, updating missing parameters requires patterns to be expanded
    pub fn update(&mut self, frequency: f64, defaults: &Parameters) {
        
        let (mut parameters, areas) = self.extract(frequency)
            .expect("Coupling could not be converted into peak parameters.");

        for p in parameters.iter_mut() {
            p.update(defaults)
        }

        *self = Peaks::Parameters(parameters, areas)
    }

    //-------------------------------------------------------------------------
    pub fn replace(&mut self, terms: &[f64], frequency: f64) {

        if terms.len()/4 != self.n_peaks() {
            panic!("Invalid number of terms.")
        }

        // Once updated, patterns are permanently transformed into parameters
        let (_, areas) = self.extract(frequency)
            .expect("Coupling could not be converted into peak parameters.");

        // Since terms are multiplied by areas, replacements are divided
        let parameters: Vec<Parameters> = terms
            .chunks(4)
            .zip(areas.iter())
            .map(|(x, a)| {
                Parameters {
                    position: x[0],
                    width: x[1]*frequency,
                    height: x[2]/a,
                    fraction: x[3]
                }
            })
            .collect();

        *self = Peaks::Parameters(parameters, areas);
    }


    //-------------------------------------------------------------------------
    pub fn extract(
        &self, 
        frequency: f64
    ) -> Result<(Vec<Parameters>, Vec<f64>)> {

        if let Peaks::Parameters(parameters, areas) = self.expand(frequency)? {
           Ok( (parameters, areas) )
        } else {
            panic!("Unreachable.")
        }
    }


    //-------------------------------------------------------------------------
    /// Expand the multiplet pattern to a set of peaks and relative areas
    /// based on NMR magnet frequency.
    pub fn expand(&self, frequency: f64) -> Result<Self> {

        // Return self if already in expanded form
        let pattern = match self {
            Peaks::Pattern(pattern) => {
                pattern.clone()
            },
            Peaks::Parameters(_, _) => {
                return Ok( self.clone() )
            }
        };

        // All of the regex
        lazy_static! {
            static ref RE_BRACKETS: Regex = 
                Regex::new(r"[(){}|\[\]]").unwrap();
            static ref RE_SPACES: Regex = 
                Regex::new(r"[ ,:#$%_/-]+").unwrap();
            static ref RE_TRIPLET: Regex = 
                Regex::new(r"(^[0-9.]+)\s([a-z]\D*)(.*$)").unwrap();
            static ref RE_NUMBERS: Regex = 
                Regex::new(r"^\d+$").unwrap();
        }

        //----------------------------------------
        // Clean up

        let pattern = pattern.to_lowercase();

        // No brackets and only spaces as separators
        let pattern = RE_BRACKETS.replace_all(&pattern, "");
        let pattern = RE_SPACES.replace_all(&pattern, " ");
        let pattern = pattern.trim();

        // Initializing generic error
        let context_msg = format!(
            "\"{}\" could not be parsed as an NMR resonance.", 
            &pattern
        );

        let f_context = || {
            context_msg.clone()
        };

        // Split into number -- letters -- numbers triplet
        let triplet = RE_TRIPLET.captures(pattern)
            .with_context(f_context)?;

        let position_str = triplet.get(1)
            .with_context(f_context)?
            .as_str().trim();
        let splits_str = triplet.get(2)
            .with_context(f_context)?
            .as_str().trim();
        let frequency_str = triplet.get(3)
            .with_context(f_context)?
            .as_str().trim();

        //----------------------------------------
        // Position 

        // Position just has to be a number
        let position: f64 = position_str.parse()
            .with_context(f_context)
            .context("Invalid ppm value.")?;

        //----------------------------------------
        // Pattern

        let translations = vec![
            ( "quintet", "5" ),
            ( "pentet", "5" ),
            ( "sextet", "6" ),
            ( "septet", "7" ),
            ( "heptet", "7" ),
            ( "quint", "5" ),
            ( "octet", "8" ),
            ( "nonet", "9" ),
            ( "pent", "5" ),
            ( "sext", "6" ),
            ( "sept", "7" ),
            ( "hept", "7" ),
            ( "pnt", "5" ),
            ( "qui", "5" ),
            ( "qnt", "5" ),
            ( "sxt", "6" ),
            ( "spt", "7" ),
            ( "hpt", "7" ),
            ( "oct", "8" ),
            ( "non", "9" ),              
            ( "s", "1" ),
            ( "d", "2" ),
            ( "t", "3" ),
            ( "q", "4" ),
        ];

        let mut splits = splits_str.to_string();

        for ( word, number ) in translations {
            splits = splits.replace(word, number);
        }

        splits = splits.replace(" ", "");

        // The end result should be a string containing only numbers
        if ! RE_NUMBERS.is_match(&splits) {
            return Err(anyhow!(context_msg).context("Invalid split pattern."))
        }

        // Singlets should be removed from splits
        let splits: Vec<usize> = splits.chars()
            .map(|x| x.to_digit(10).unwrap() as usize )
            .filter(|&x| x > 1 )
            .collect(); 

        //----------------------------------------
        // Frequency is just split by spaces
        
        let mut frequencies: Vec<f64> = Vec::new();

        if frequency_str != "" {

            frequencies = frequency_str.split(" ")
                .map(|x| x.parse().context("Could not convert to number") )
                .collect::<Result<Vec<_>>>()
                .with_context(f_context)
                .context("Invalid frequency.")?;
        }

        //----------------------------------------
        // Putting everything together
        
        // Each split should have a corresponding frequency
        if frequencies.len() < splits.len() {
            let err = anyhow!(context_msg)
                .context("Not enough split frequencies.");
            return Err(err)
        }

        if frequencies.len() > splits.len() {
            let err = anyhow!(context_msg)
                .context("Too many split frequencies.");
            return Err(err)
        }

        // Start with singlet at specified position
        let mut positions: Vec<f64> = vec![position];
        let mut areas: Vec<f64> = vec![1.0];

        let iterator = splits.iter().zip(frequencies.iter());
        for (number, coupling_frequency) in iterator {
            
            let coupling_frequency = *coupling_frequency/frequency;

            let offsets = itertools_num::linspace::<f64>(-0.5, 0.5, *number)
                .map(|x| x * (*number as f64 - 1.0) * coupling_frequency );

            let number: u64 = (*number).try_into()
                .expect("Unexpected integer conversion error.");

            let mut ratios: Vec<f64> = (0..(number))
                .map(|x| factorial::binomial(number - 1, x) )
                .collect();

            let ratio_sum: f64 = ratios.iter().sum();
            ratios.iter_mut().map(|x| *x = *x / ratio_sum).count();

            let ratios = ratios.iter();

            positions = positions.iter()
                .map(|x| {
                    offsets.clone().map(|y| x + y).collect::<Vec<f64>>()
                })
                .flatten()
                .collect();

            areas = areas.iter()
                .map(|x| {
                    ratios.clone().map(|y| x * y).collect::<Vec<f64>>()
                })
                .flatten()
                .collect();
        }


        Peaks::from_positions(positions, areas)
    }
}


//=============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize)]
#[serde(from = "PeaksTemp")] 
pub enum PeaksTemplate {
    Label(String),
    Pattern(String),
    Positions(Vec<f64>),
    Areas((Vec<f64>, Vec<f64>)),
}


impl PeaksTemplate {

    pub fn to_string(&self) -> String {

        match self {
            PeaksTemplate::Label(label) => label.clone(),
            PeaksTemplate::Pattern(pattern) => pattern.clone(),
            PeaksTemplate::Positions(positions) => {
                positions.iter()
                    .map(|x| format!("{:.4}", x))
                    .collect::<Vec<String>>()
                    .join(" ")
            },
            PeaksTemplate::Areas((positions, _)) => {
                positions.iter()
                    .map(|x| format!("{:.4}", x))
                    .collect::<Vec<String>>()
                    .join(" ")
            }
        }
    }


    pub fn from_string(text: String) -> Result<Self> {

        let first = text.clone().chars().next().unwrap();

        if first.is_digit(10) {
            return Ok( PeaksTemplate::Pattern(text) )
        } else {
            return Ok( PeaksTemplate::Label(text) )
        }

        bail!("Cannot iterpret \"{}\" as resonance.", text);
    }
}


impl From<PeaksTemp> for PeaksTemplate {

    fn from(temp: PeaksTemp) -> Self {

        match temp {
            PeaksTemp::Float(content) => {
                PeaksTemplate::Positions(vec![content])
            },
            PeaksTemp::Text(text) => {

                let first = text.clone().chars().next().unwrap();

                if first.is_digit(10) {
                    PeaksTemplate::Pattern(text)
                } else {
                    PeaksTemplate::Label(text)
                }
            },
            PeaksTemp::Positions(content) => {
                PeaksTemplate::Positions(content)
            },
            PeaksTemp::Areas(content) => {
                PeaksTemplate::Areas(content)
            }
        }
    }
}


//=============================================================================
// Intermediate parsing structure

#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)] 
enum PeaksTemp {
    Float(f64),
    Text(String),
    Positions(Vec<f64>),
    Areas((Vec<f64>, Vec<f64>)),
}









