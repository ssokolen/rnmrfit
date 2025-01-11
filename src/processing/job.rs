use anyhow::{Context, Result, bail};
use glob::Pattern;
use getset::{Getters, MutGetters, Setters}; 
use itertools::Itertools;
use log::{info, warn};
use polars::prelude::*;
use serde::Deserialize;
use std::{
    fs, mem,
    collections::HashMap,
    fs::File,
    ops::Not,
    path::PathBuf,
};

use crate::processing;
use crate::processing::{*, Bounds};


//=============================================================================

#[derive(Debug, Default, Getters, MutGetters, Setters)]
#[getset(get = "pub")]
pub struct Job1D {

    name: String,

    resonances: HashMap<String, Resonance>,
    species: HashMap<String, Species1D>,
    #[getset(get_mut = "pub")]
    experiments: HashMap<String, Experiment1D>,

    processing: Processing1DTemplate,

    baseline: Baseline1DTemplate,
    phase: Phase1DTemplate,

    n_fits: usize,
    n_plots: usize,
    n_area_calculations: usize,
    n_lists: usize,

    parameters: Parameters,
    general_bounds: Bounds,
    offset_bounds: Bounds,
    ranges: Ranges,
    
    tasks: Vec<Task>,
    trace: Trace,
    errors: bool,
}


//=============================================================================
// Constructor

impl Job1D {

    //-------------------------------------------------------------------------
    pub fn from_yaml_file(path: &PathBuf) -> Result<Self> {

        info!("Parsing \"{}\"", path.to_string_lossy());

        let name = path.file_stem()
            .expect("Unexpected error parsing job name.");

        let mut job = Job1D {
            name: name.to_string_lossy().to_string(),
            parameters: DEFAULT_PARAMETERS,
            general_bounds: DEFAULT_GENERAL_BOUNDS,
            offset_bounds: DEFAULT_OFFSET_BOUNDS,
            ..Default::default()
        };

        let entry_strings = processing::path_to_yaml_docs(path)?
            .into_iter()
            .map(|doc| processing::yaml_to_string(&doc));

        for yaml in entry_strings {
            job.parse(&yaml);
        }

        if job.errors {
            bail!("One or more errors in YAML definitions.");
        }

        Ok( job )
    }
}


//=============================================================================
// General parsing

impl Job1D {

    //-------------------------------------------------------------------------
    pub fn parse(&mut self, yaml: &String) {

        let entry: Result<Entry, _> = serde_yaml::from_str(yaml);

        if entry.is_err() {
            ParsingError::no_entry(yaml.to_string());
            self.errors = true;
            return
        }

        // What gets parsed depends on the last entry
        match self.trace {
            Trace::None => self.parse_settings(yaml),
            Trace::Settings => self.parse_settings(yaml),
            Trace::Spectra => self.parse_spectra(yaml),
            Trace::Task => self.parse_task(yaml),
        }
    }

    //-------------------------------------------------------------------------
    // Parse settings includes settings, spectra, and tasks
    pub fn parse_settings(&mut self, yaml: &String) {

        let entry: Result<Entry, _> = serde_yaml::from_str(yaml);

        let mut settings = match entry.unwrap() {
            Entry::Settings(settings) => settings,
            Entry::Spectra(_) => {
                self.parse_spectra(yaml);
                return
            },
            Entry::Task(_) => {
                self.parse_task(yaml);
                return
            }
        };

        self.processing = settings.processing;

        settings.parameters.update(&self.parameters);
        self.parameters = settings.parameters;

        self.general_bounds.update(&settings.general_bounds);
        self.offset_bounds.update(&settings.offset_bounds);

        self.ranges = settings.ranges;

        self.baseline = settings.baseline;
        self.phase = settings.phase;

        self.trace = Trace::Settings;
    }


    //-------------------------------------------------------------------------
    // Parse spectra includes spectra, and tasks
    pub fn parse_spectra(&mut self, yaml: &String) {

        let entry: Entry = serde_yaml::from_str(yaml).unwrap();

        let spectra = match entry {
            Entry::Settings(_) => {
                let _ = ParsingError::settings_order(yaml.to_string());
                self.errors = true;
                return
            },
            Entry::Spectra(spectra) => spectra,
            Entry::Task(_) => {
                self.parse_task(yaml);
                return
            }
        };

        self.trace = Trace::Spectra;

        // Spectra includes resonances, species, and data
        match spectra {

            //----------------------------------------
            // Resonance
            SpectraEntry::Resonance(resonance) => {

                let resonance = match Resonance::from_template(resonance) {
                    Ok(resonance) => resonance,
                    Err(e) => {
                        ParsingError::malformed_entry(
                            "Resonance", e, yaml.to_string()
                        );
                        self.errors = true;
                        return
                    }
                };

                let name = resonance.name().clone();
                if self.resonances.contains_key(&name) {
                    warn!(
                        "Resonance \"{}\" is already defined, ignoring.", &name
                    );
                } else {
                    self.resonances.insert(name, resonance);
                }

            },

            //----------------------------------------
            // Species
            SpectraEntry::Species1D(species) => {

                let mut species = match 
                    Species1D::from_template(species, self) {

                    Ok(species) => species,
                    Err(e) => {
                        ParsingError::malformed_entry(
                            "Species", e, yaml.to_string()
                        );
                        self.errors = true;
                        return
                    }
                };

                // Update default parameters and bounds
                let mut parameters = species.parameters().clone();
                parameters.update(&self.parameters);
                species.set_parameters(parameters);

                let mut general_bounds = species.general_bounds().clone();
                general_bounds.update(&self.general_bounds);
                species.set_general_bounds(general_bounds);

                let mut offset_bounds = species.offset_bounds().clone();
                offset_bounds.update(&self.offset_bounds);
                species.set_offset_bounds(offset_bounds);

                let name = species.name().clone();
                if self.species.contains_key(&name) {
                    warn!("Overwriting species \"{}\".", &name);
                }

                self.species.insert(name, species);

                
            },

            //----------------------------------------
            // Data 
            SpectraEntry::Data1D(data) => {
                
                let mut data = match Data1D::from_template(data, self) {
                    Ok(data) => data,
                    Err(e) => {
                        ParsingError::malformed_entry(
                            "Dataset", e, yaml.to_string()
                        );
                        self.errors = true;
                        return
                    }
                };
                
                for experiment in data.experiments_mut().drain(..) {
                    
                    let name = experiment.name().clone();

                    if self.experiments.contains_key(&name) {
                        warn!("Overwriting experiment \"{}\".", &name);
                    }

                    self.experiments.insert(name, experiment);
                }

            }
        }

    }

    //-------------------------------------------------------------------------
    // Parse tasks includes only tasks
    pub fn parse_task(&mut self, yaml: &String) {

        let entry: Entry = serde_yaml::from_str(yaml).unwrap();

        let task = match entry {
            Entry::Settings(_) => {
                let _ = ParsingError::settings_order(yaml.to_string());
                self.errors = true;
                return
            },
            Entry::Spectra(_) => {
                let _ = ParsingError::spectra_order(yaml.to_string());
                self.errors = true;
                return
            },
            Entry::Task(task) => task 
        };

        self.trace = Trace::Task;

        // Task includes fit, plot, calculation, list
        match task {

            //----------------------------------------
            // Fit 
            TaskEntry::Fit1D(fit) => {
                let fit = match Fit1D::from_template(fit, self) {
                    Ok(fit) => fit,
                    Err(e) => {
                        ParsingError::malformed_entry(
                            "Fit", e, yaml.to_string()
                        );
                        self.errors = true;
                        return
                    }
                };

                self.n_fits += 1;
                self.tasks.push(Task::Fit(fit));
            },

            //----------------------------------------
            // Plot
            TaskEntry::Plot1D(plot) => {
                let plot = match Plot1D::from_template(plot, self) {
                    Ok(plot) => plot,
                    Err(e) => {
                        ParsingError::malformed_entry(
                            "Plot", e, yaml.to_string()
                        );
                        self.errors = true;
                        return
                    }
                };

                self.n_plots += 1;
                self.tasks.push(Task::Plot(plot));
            },

            //----------------------------------------
            // Calculation
            TaskEntry::CalculateArea1D(calc) => {

                let calc = match AreaCalculation1D::from_template(calc, self) {
                    Ok(calc) => calc,
                    Err(e) => {
                        ParsingError::malformed_entry(
                            "Area calculation", e, yaml.to_string()
                        );
                        self.errors = true;
                        return
                    }
                };

                self.n_area_calculations += 1;
                self.tasks.push(Task::CalculateArea(calc));
            }

            //----------------------------------------
            // List 
            TaskEntry::ListParameters1D(list) => {

                let list = match ParameterList1D::from_template(list, self) {
                    Ok(list) => list,
                    Err(e) => {
                        ParsingError::malformed_entry(
                            "Parameter list", e, yaml.to_string()
                        );
                        self.errors = true;
                        return
                    }
                };

                self.n_lists += 1;
                self.tasks.push(Task::ListParameters(list));
            }
        }
    }
}


//=============================================================================
// Glob-related functions

impl Job1D {

    //-------------------------------------------------------------------------
    pub fn get_species_from_globs(
        &self, 
        patterns: &Vec<Pattern>
    ) -> Vec<String> {

        patterns.iter()
            .map(|x| {
                self.species.keys()
                    .filter(move |y| x.matches(y))
            })
            .flatten()
            .unique()
            .map(|x| x.clone() )
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn get_experiments_from_globs(
        &self, 
        patterns: &Vec<Pattern>
    ) -> Vec<String> {

        patterns.iter()
            .map(|x| {
                self.experiments.keys()
                    .filter(move |y| x.matches(y))
            })
            .flatten()
            .unique()
            .map(|x| x.clone() )
            .collect()
    }


    //-------------------------------------------------------------------------
    // Identify missing species and bail if there are any
    pub fn check_species_globs(
        &self, 
        patterns: &Vec<Pattern>
    ) -> Result<()> {

        let missing: Vec<String> = patterns.iter()
            .filter(|x| {
                self.species()
                    .keys()
                    .any(|y| x.matches(y))
                    .not()
            })
            .map(|x| x.to_string())
            .collect();

        if missing.len() > 0 {
            let missing = missing.join(", ");
            bail!("Unmatched species: \"{}\"", missing)
        }

        Ok( () )
    }


    //-------------------------------------------------------------------------
    // Identify missing experiments and bail if there are any
    pub fn check_experiment_globs(
        &self, 
        patterns: &Vec<Pattern>
    ) -> Result<()> {

        let missing: Vec<String> = patterns.iter()
            .filter(|x| {
                self.experiments.keys()
                    .any(|y| x.matches(y))
                    .not()
            })
            .map(|x| x.to_string())
            .collect();

        if missing.len() > 0 {
            let missing = missing.join(", ");
            bail!("Unmatched experiments: \"{}\"", missing)
        }

        Ok( () )
    }
}


//=============================================================================
// Main processing code

impl Job1D {

    //-------------------------------------------------------------------------
    pub fn process(&mut self, path: &PathBuf, tol: f64) -> Result<()> {

        let mut path = path.clone();
        path.push(&self.name);

        fs::create_dir_all(&path)
            .with_context(|| {
                format!(
                    "Could not create output directory for \"{}\"",
                    &self.name
                )
            })?;

        let tasks = mem::replace(&mut self.tasks, Vec::new());
        
        for mut task in tasks.into_iter() {

            match task {
                Task::Fit(ref mut fit) => {
                    fit.fit(self, tol)?;
                },
                Task::CalculateArea(ref calculation) => {

                    let mut calculation_path = path.clone();
                    calculation_path.push(
                        format!("{}.csv", calculation.name())
                    );

                    let mut file = File::create(calculation_path)
                        .expect("Error creating csv file");

                    // Generate DataFrame
                    let mut calculation = calculation.calculate(&self);

                    // and write it to file
                    CsvWriter::new(&mut file)
                        .include_header(true)
                        .with_separator(b',')
                        .finish(&mut calculation)
                        .expect("Error writing csv file");

                },
                Task::ListParameters(ref list) => {

                    let mut list_path = path.clone();
                    list_path.push(
                        format!("{}.csv", list.name())
                    );

                    let mut file = File::create(list_path)
                        .expect("Error creating csv file");

                    // Generate DataFrame
                    let mut list = list.list(&self);

                    // and write it to file
                    CsvWriter::new(&mut file)
                        .include_header(true)
                        .with_separator(b',')
                        .finish(&mut list)
                        .expect("Error writing csv file");

                },
                Task::Plot(ref plot) => {

                    // The following must account for multiple possible plots
                    let mut paths: Vec<PathBuf> = Vec::new();
                    let plot_name = plot.name();
                        
                    // File paths can be generated before plots 
                    if false { //*plot.split_experiments() {

                        /*
                        for name in plot.get_experiment_names() {

                            let mut plot_path = path.clone();
                            plot_path.push(
                                format!("{}_{}.html", &plot_name, &name)
                            );

                            paths.push(plot_path);
                        }
                        */

                    } else {
                        
                        let mut plot_path = path.clone();
                        plot_path.push(
                            format!("{}.html", plot_name)
                        );
                        
                        paths.push(plot_path);
                    }

                    let plots = plot.plot(&self);

                    for (plot, path) in plots.into_iter().zip(paths.iter()) {
                        plot.to_html(path);
                    }
                }
            }
        }

        Ok( () )
    }
}


//=============================================================================
// Tasks stored by job

#[derive(Debug)]
enum Task {
    Plot(Plot1D),
    Fit(Fit1D),
    CalculateArea(AreaCalculation1D),
    ListParameters(ParameterList1D),
}


//=============================================================================
// Yaml template for generic entry

#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)] 
enum Entry {
    Settings(SettingsEntry),
    Spectra(SpectraEntry),
    Task(TaskEntry)
}


// We need to keep track of which entry was parsed last to ensure order
#[derive(Clone, Debug)]
enum Trace {
    None,
    Settings,
    Spectra,
    Task
}

impl Default for Trace {

    fn default() -> Self {
        Trace::None
    }
}


//-----------------------------------------------------------------------------
// Spectra and Task entries are just collections of existing templates

#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum SpectraEntry {
    Resonance(ResonanceTemplate),
    Species1D(Species1DTemplate),
    Data1D(Data1DTemplate)
}


#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum TaskEntry {
    Fit1D(Fit1DTemplate),
    Plot1D(Plot1DTemplate),
    CalculateArea1D(AreaCalculation1DTemplate),
    ListParameters1D(ParameterList1DTemplate)
}


//-----------------------------------------------------------------------------
// Settings entry has to be manually defined

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SettingsEntry {

    #[serde(default)] 
    processing: Processing1DTemplate,

    #[serde(default)] 
    parameters: Parameters,
    #[serde(default)] 
    general_bounds: Bounds,
    #[serde(default)] 
    offset_bounds: Bounds,
    #[serde(default)]
    ranges: Ranges,

    #[serde(default)]
    baseline: Baseline1DTemplate,
    #[serde(default)]
    phase: Phase1DTemplate,
}


#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum Ranges {
    Manual(ManualRanges)
}


impl Default for Ranges {

    fn default() -> Self {
        Ranges::Manual(ManualRanges::default())
    }
}


#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ManualRanges {

    #[serde(default = "default_cutoff")] 
    pub cutoff: f64,
}


impl Default for ManualRanges {

    fn default() -> Self {
        ManualRanges { cutoff: default_cutoff() }
    }
}


fn default_cutoff() -> f64 {
    5e-2
}



