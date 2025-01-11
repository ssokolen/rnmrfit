use anyhow::{Result};
use getset::{Getters, MutGetters, Setters};
use itertools::{iproduct, multiunzip};
use polars::prelude::*;
use serde::Deserialize;
use std::{collections::HashSet, iter};


use crate::processing::{
    Globs, Job1D, Scaffold
};


// Unbounded ppm range
const UNBOUNDED: [f64; 2] = [f64::NEG_INFINITY, f64::INFINITY];


//=============================================================================
// AreaCalculation1D


#[derive(Debug, Default, Getters, MutGetters, Setters)]
pub struct AreaCalculation1D {

    #[getset(get = "pub")]
    name: String,

    experiments: Vec<String>,
    species: Vec<String>,
    
    sum_by: Scaffold,
    normalize_by: String,

    summarize: bool,
}


//=============================================================================
// YAML parsing

impl AreaCalculation1D {

    //-------------------------------------------------------------------------
    pub fn from_template(
        template: AreaCalculation1DTemplate, 
        job: &Job1D
    ) -> Result<Self> {

        // Ensuring that species and experiments exist
        let species = template.species.patterns()?;
        job.check_species_globs(&species)?;

        let species = job.get_species_from_globs(&species);

        let experiments = template.experiments.patterns()?;
        job.check_experiment_globs(&experiments)?;

        let experiments = job.get_experiments_from_globs(&experiments);

        let area_calculation = AreaCalculation1D {
            name: template.calculate_areas,
            experiments,
            species,
            sum_by: template.sum_by,
            normalize_by: template.normalize_by,
            summarize: template.summarize
        };

        Ok( area_calculation )
    }
}


//=============================================================================
// Calculation

type AreaRow = (
    String, String, String, String, f64
);

type AreaColumns = (
    Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<f64>
);

impl AreaCalculation1D {

    //-------------------------------------------------------------------------
    pub fn calculate(&self, job: &Job1D) -> DataFrame {

        // Initialize area rows
        let mut rows: Vec<AreaRow> = vec![];

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
                
                // Only calculate peak areas if the peak has been fit,
                // otherwise generate NAN values
                let peak_areas: Vec<f64> = match fitted {
                    true => resonance.peak_areas(frequency),
                    false => {
                        iter::repeat(f64::NAN)
                            .take(peak_names.len())
                            .collect()
                    }
                };

                // Add entries as new rows
                let iterator = peak_names.into_iter()
                    .zip(peak_areas.into_iter());

                for (peak, area) in iterator {
                    rows.push(
                        (
                            e.clone(), 
                            s.clone(),
                            resonance.name().clone(),
                            peak.clone(),
                            area 
                        )
                    )
                }
            }
        }

        let columns: AreaColumns = multiunzip(rows);

        macro_rules! column {
            ( $index:tt ) => {
                columns.$index.into_iter().collect::<Series>()
            };
        }

        // Generate dataframe
        let mut df = df!(
            "Experiment" => column!(0),
            "Species" => column!(1),
            "Resonance" => column!(2),
            "Peak" => column!(3),
            "Area" => column!(4),
        ).unwrap();

        // Sum
        let groups = match self.sum_by {
            Scaffold::Species => vec!["Species"],
            Scaffold::Resonance => vec!["Species", "Resonance"],
            Scaffold::Peak => vec!["Species", "Resonance", "Peak"],
        };

        let mut columns = groups.clone();
        columns.insert(0, "Experiment");

        df = df.lazy()
            .group_by(&columns)
            .agg([ col("Area").sum() ])
            .collect()
            .unwrap();

        // Output name may be Area or Ratio depending on normalization
        let mut output = "Area";

        // Normalize
        if self.normalize_by != "" {

            output = "Ratio";

            let normalize_by = self.normalize_by.clone();

            // Convert empty strings into nulls
            let mut columns: Vec<_> = groups.iter()
                .map(|x| col(x.clone()) )
                .collect();

            let norm = df.clone().lazy()
                .with_column(
                    concat_str(&columns, "::", true).alias("Reference")
                )
                .filter( col("Reference").eq(lit(normalize_by.clone())) )
                .select([
                    col("Experiment"), 
                    col("Area").alias("Reference")]
                )
                .collect()
                .unwrap();

            // Extend selection
            let mut columns: Vec<_> = groups.iter()
                .map(|x| col(x.clone()))
                .collect();
            columns.insert(0, col("Experiment"));
            columns.push(( col("Area") / col("Reference") ).alias("Ratio"));

            // Join and normalize
            df = df.lazy()
                .join(
                    norm.lazy(),
                    [col("Experiment")],
                    [col("Experiment")],
                    JoinArgs::new(JoinType::Left),
                )
                .select(&columns)
                .collect()
                .unwrap();

        }

        // Add basic summary across all datasets
        if self.summarize {

            let mut columns: Vec<_> = groups.iter()
                .map(|x| col(x.clone()))
                .collect();

            columns.insert(0, col("Dataset"));

            df = df.lazy()
                .with_column(
                    col("Experiment")
                        .str()
                        .split(lit("::"))
                        .list()
                        .first()
                        .alias("Dataset")
                )
                .group_by(columns)
                .agg([
                    col(&output)
                        .mean()
                        .alias("mean"),
                    col(&output)
                        .std(1)
                        .alias("std"),
                    col(&output)
                        .min()
                        .alias("min"),
                    col(&output)
                        .median()
                        .alias("med"),
                    col(&output)
                        .max()
                        .alias("max"),
                ])
                .select([
                    col("Dataset"),
                    cols(&groups),
                    col("mean")
                        .alias(format!("{} mean", output).as_str()),
                    col("std")
                        .alias(format!("{} std. dev.", output).as_str()),
                    ( col("std") / col("mean") * lit(100.0) )
                        .alias(format!("{} coeff. var. (%)", output).as_str()),
                    col("min")
                        .alias(format!("{} minimum", output).as_str()),
                    col("med")
                        .alias(format!("{} median", output).as_str()),
                    col("max")
                        .alias(format!("{} maximum", output).as_str()),
                ])
                .collect()
                .unwrap();
        }

        // Sort and output
        let mut columns = groups.clone();
        if self.summarize {
            columns.insert(0, "Dataset");       
        } else {
            columns.insert(0, "Experiment");
        }

        df.lazy()
            .sort(&columns, Default::default())
            .collect()
            .unwrap()
    }
}


//=============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct AreaCalculation1DTemplate {

    calculate_areas: String,

    #[serde(default)] 
    experiments: Globs,
    #[serde(default)]
    species: Globs,

    #[serde(default = "default_sum_by")]
    sum_by: Scaffold,
    #[serde(default)] 
    normalize_by: String,

    #[serde(default)] 
    summarize: bool,
}


//-----------------------------------------------------------------------------
fn default_sum_by() -> Scaffold {
    Scaffold::Resonance
}
