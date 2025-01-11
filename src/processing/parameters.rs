use anyhow::{Result};
use getset::{Getters, MutGetters, Setters};
use itertools::{Itertools, iproduct, multiunzip};
use polars::prelude::*;
use serde::Deserialize;
use std::{collections::HashSet, iter};


use crate::processing::{
    Globs, Job1D,
};


// Unbounded ppm range
const UNBOUNDED: [f64; 2] = [f64::NEG_INFINITY, f64::INFINITY];


//=============================================================================
// ParameterList1D


#[derive(Debug, Default, Getters, MutGetters, Setters)]
#[getset(get = "pub")]
pub struct ParameterList1D {

    #[getset(set = "pub")]
    name: String,

    experiments: Vec<String>,
    species: Vec<String>,
}


//=============================================================================
// YAML parsing

impl ParameterList1D {

    //-------------------------------------------------------------------------
    pub fn from_template(
        template: ParameterList1DTemplate, 
        job: &Job1D
    ) -> Result<Self> {

        // Ensuring that species and experiments exist
        let species = template.species.patterns()?;
        job.check_species_globs(&species)?;

        let species = job.get_species_from_globs(&species);

        let experiments = template.experiments.patterns()?;
        job.check_experiment_globs(&experiments)?;

        let experiments = job.get_experiments_from_globs(&experiments);

        let parameter_list = ParameterList1D {
            name: template.list_parameters,
            experiments,
            species,
        };

        Ok( parameter_list )
    }
}


//=============================================================================
// Listing 

type ParameterRow = (
    String, String, String, String, 
    f64, f64, f64, f64
);

type ParameterColumns = (
    Vec<String>, Vec<String>, Vec<String>, Vec<String>, 
    Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>
);


impl ParameterList1D {

    //-------------------------------------------------------------------------
    pub fn list(&self, job: &Job1D) -> DataFrame {

        // Initialize parameter rows
        let mut rows: Vec<ParameterRow> = vec![];

        // Loop over all combinations of experiments and species
        for (e, s) in iproduct!(self.experiments.iter(), self.species.iter()) {

            let experiment = job.experiments().get(e)
                .expect("Missing experiment.");

            let frequency = *experiment.frequency();

            let fitted = experiment.species().contains_key(s);

            let species = experiment.species().get(s)
                .or(job.species().get(s))
                .expect("Missing species");

            // Peak calculations are performed at the resonance level
            for resonance in species.resonances(frequency, &UNBOUNDED) {

                let peak_names = resonance.peak_names();
                
                // Only calculated peak areas if the peak has been fit,
                // otherwise generate NAN values
                let peak_terms: Vec<f64> = match fitted {
                    true => resonance.terms(frequency),
                    false => {
                        iter::repeat(f64::NAN)
                            .take(peak_names.len()*4)
                            .collect()
                    }
                };

                // Chunk parameters into position/width/height/fraction
                let peak_terms = peak_terms.into_iter()
                    .chunks(4);

                // Add entries as new rows
                let iterator = peak_names.into_iter()
                    .zip(peak_terms.into_iter());

                for (peak, mut terms) in iterator {
                    rows.push(
                        (
                            e.clone(), 
                            s.clone(),
                            resonance.name().clone(),
                            peak.clone(),
                            terms.next().unwrap(),
                            // Width needs to be converted into Hz
                            terms.next().unwrap()*frequency,
                            terms.next().unwrap(),
                            terms.next().unwrap()
                        )
                    )
                }
            }
        }

        let columns: ParameterColumns = multiunzip(rows);

        macro_rules! column {
            ( $index:tt ) => {
                columns.$index.into_iter().collect::<Series>()
            };
        }

        // Generate dataframe
        let df = df!(
            "Experiment" => column!(0),
            "Species" => column!(1),
            "Resonance" => column!(2),
            "Peak" => column!(3),
            "Position" => column!(4),
            "Width" => column!(5),
            "Height" => column!(6),
            "Fraction" => column!(7),
        ).unwrap();

        let columns = ["Experiment", "Species", "Resonance", "Peak"];

        df.lazy()
            .sort(columns, Default::default())
            .collect()
            .unwrap()
    }
}


//=============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct ParameterList1DTemplate {

    list_parameters: String,

    #[serde(default)] 
    experiments: Globs,
    #[serde(default)]
    species: Globs,
}
