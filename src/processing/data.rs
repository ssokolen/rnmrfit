use anyhow::{Result, bail};
use getset::{Getters, MutGetters};
use serde::Deserialize;
use std::{
    f64::consts::PI,
    iter,
    path::PathBuf,
};

use crate::processing::*;
use crate::importing::dir_names;


//=============================================================================
// Data1D

#[derive(Clone, Debug, Default, MutGetters)]
pub struct Data1D {

    #[getset(get_mut = "pub")]
    experiments: Vec<Experiment1D>,
}


impl Data1D {


//=============================================================================
// YAML parsing

    //------------------------------------------------------------------------
    pub fn from_template(
        template: Data1DTemplate, job: &Job1D,
    ) -> Result<Data1D> {

        let dataset = template.dataset;

        // First, check that provided path exists
        let path = PathBuf::from(template.path);

        if ! path.exists() {
            bail!("Cannot access \"{}\"", path.display());
        }

        // Generate glob patterns
        let patterns = template.experiments.patterns()?;

        // Generate file names and run them against patterns
        let filenames: Vec<String> = dir_names(&path)?
            .into_iter()
            .filter(|x| {
                patterns.iter()
                    .any(|y| y.matches(x)) 
            })
            .collect();

        let mut processing = template.processing;
        processing.update(job.processing());

        // Extract which procs dir to use (if any)
        let procs = match &processing {
            Processing1DTemplate::None => Some("".to_string()),
            Processing1DTemplate::Automatic(dir) => Some(dir.clone()),
            Processing1DTemplate::Manual(_) => None,
        };

        type Finisher = Box<dyn Fn(&mut Experiment1D)>;

        // List of finishing functions 
        let mut finishers: Vec<Finisher> = Vec::new();

        // General processing (if not using TopSping procs)
        if let Processing1DTemplate::Manual(mut processing) = processing {

            processing.update(&DEFAULT_PROCESSING);

            let phase = processing.phase.unwrap();
            if phase {
                finishers.push(Box::new(move |experiment: &mut Experiment1D| {
                    experiment.autophase();
                }));
            }

            let zf = processing.zf.unwrap();

            finishers.push(Box::new(move |experiment: &mut Experiment1D| {
                experiment.zero_fill(zf);
            }));

            let lb = processing.lb.unwrap();

            finishers.push(Box::new(move |experiment: &mut Experiment1D| {
                experiment.line_broaden(lb);
            }));

            let mut reference = processing.reference;
            reference.update(&DEFAULT_REFERENCE);
            let ppm = reference.ppm.unwrap();
            let threshold = reference.threshold.unwrap();

            finishers.push(Box::new(move |experiment: &mut Experiment1D| {
                experiment.reference_align(ppm, threshold);
            }));
        }

        // Baseline 
        let mut baseline = template.baseline;

        // First, update with general settings, then defaults
        baseline.update(job.baseline());
        baseline.update(&DEFAULT_BASELINE);

        let span = baseline.span().unwrap();
        let degree = baseline.degree().unwrap();
        let bounds = baseline.bounds().clone();

        finishers.push(Box::new(move |experiment: &mut Experiment1D| {
            experiment.set_baseline(span, degree, bounds);
        }));

        // Phase
        let mut phase = template.phase;

        // First, update with general settings, then defaults
        phase.update(job.phase());
        phase.update(&DEFAULT_PHASE);

        let order = phase.order().unwrap();
        // Converting to radians from degrees
        let bounds = phase.bounds().unwrap() * PI/180.0; 

        finishers.push(Box::new(move |experiment: &mut Experiment1D| {
            experiment.set_phase(order, bounds);
        }));

        let experiments: Vec<Experiment1D> = filenames.into_iter()
            .map(|filename| {
                let mut path = path.clone();
                path.push(&filename);

                // Ensure that path exists
                if ! path.exists() {
                    bail!("Cannot access \"{}\"", path.display())
                }

                let mut experiment = Experiment1D::from_path(&path, &procs)
                    .and_then(|mut x| {
                        x.set_name( format!("{}::{}", &dataset, &filename) );
                        Ok( x )
                    });

                // If everything loaded fine, apply finishers
                if let Ok(experiment) = experiment.as_mut() {

                    for f in finishers.iter() {
                        f(experiment);
                    }
                }

                experiment
            })
            .collect::<Result<Vec<Experiment1D>>>()?
            .into_iter()
            .collect();

        let data = Data1D::new(experiments);

        Ok( data )
    }
}


//=============================================================================
// Constructors

impl Data1D {

    //-------------------------------------------------------------------------
    pub fn new(experiments: Vec<Experiment1D>) -> Data1D {

        Data1D {
            experiments: experiments,
            ..Default::default()
        }
    }
}


//=============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct Data1DTemplate {

    dataset: String,
    path: String,

    #[serde(default)] 
    experiments: Globs,
    
    #[serde(default)] 
    processing: Processing1DTemplate,

    #[serde(default)]
    baseline: Baseline1DTemplate,
    #[serde(default)]
    phase: Phase1DTemplate,
}




//-----------------------------------------------------------------------------
#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum Processing1DTemplate {
    None,
    Automatic(String),
    Manual(ManualProcessing1DTemplate)
}


impl Default for Processing1DTemplate {

    fn default() -> Self {
        Processing1DTemplate::None
    }
}


impl Processing1DTemplate {

    pub fn update(&mut self, template: &Processing1DTemplate) {

        if matches!(self, Processing1DTemplate::None) {
            *self = template.clone()
        }
    }
}


#[derive(Clone, Debug, Default, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct ManualProcessing1DTemplate {

    #[serde(default)]
    lb: Option<f64>,
    #[serde(default)]
    zf: Option<f64>,
    #[serde(default)]
    phase: Option<bool>,
    #[serde(default)]
    reference: Reference1DTemplate,
}


impl ManualProcessing1DTemplate {

    pub fn update(&mut self, template: &ManualProcessing1DTemplate) {

        if self.lb.is_none() {
            self.lb = *template.lb();
        }

        if self.zf.is_none() {
            self.zf = *template.zf();
        }

        if self.phase.is_none() {
            self.phase = *template.phase();
        }

        self.reference.update(&template.reference);
    }

}


// Sensible global defaults
pub const DEFAULT_PROCESSING: ManualProcessing1DTemplate = 
    ManualProcessing1DTemplate {
        
        lb: Some(0.0),
        zf: Some(0.0),
        phase: Some(true),
        reference: DEFAULT_REFERENCE
};



//-----------------------------------------------------------------------------
#[derive(Clone, Debug, Default, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct Reference1DTemplate {

   ppm: Option<f64>,
   threshold: Option<f64>,

}


impl Reference1DTemplate {

    pub fn update(&mut self, template: &Reference1DTemplate) {

        if self.ppm.is_none() {
            self.ppm = *template.ppm();
        }

        if self.threshold.is_none() {
            self.threshold = *template.threshold();
        }
    }
}


// Sensible global defaults
pub const DEFAULT_REFERENCE: Reference1DTemplate = 
    Reference1DTemplate {
    
        ppm: Some(f64::NAN),
        threshold: Some(0.2),
};
