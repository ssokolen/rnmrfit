use anyhow::{Result, anyhow, bail};
use cascade::cascade;
use colorbrewer::*;
use css_color_parser::Color as CssColor;
use getset::{Getters, MutGetters, Setters};
use rgb::RGB;
use rustfft::{FftPlanner, num_complex::Complex};
use serde::Deserialize;
use std::convert::TryFrom;
use std::f64::consts::PI;

use plotly::common::{Line, Mode, Title};
use plotly::common::color::Rgb;
use plotly::layout::{Axis, Layout};
use plotly::{Plot, Scatter};

use crate::processing::{
    Globs, Parameters, Job1D, //Scaffold
};





//=============================================================================
// Misc structs

#[derive(Clone, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Component {
    R,
    I,
}


impl Default for Component {

    fn default() -> Self {
        Component::R
    }
}


#[derive(Clone, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Apodization {
    Identity,
    Sine,
}


impl Default for Apodization {

    fn default() -> Self {
        Apodization::Identity
    }
}


impl Apodization {

    fn apodize(&self, y: &mut Vec<Complex<f64>>) {

        match self {
            Apodization::Identity => {},
            Apodization::Sine => {
                Apodization::fft(y);

                let n: f64 = (y.len() - 1) as f64;

                for (i, value) in y.iter_mut().enumerate() {
                    let scalar: f64 = (i as f64)/n;
                    *value = *value * (scalar * PI).sin();
                }

                Apodization::ifft(y);
            }
        }
    }

    fn fft(y: &mut Vec<Complex<f64>>) {

        // Take FFT
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(y.len());
        fft.process(y);
    }

    fn ifft(y: &mut Vec<Complex<f64>>) {

        // Take FFT
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_inverse(y.len());
        fft.process(y);
    }
}


//=============================================================================
// Helper functions

fn parse_color(color: String, n: usize) -> Result<Vec<RGB<u8>>> {

    // First, try parsing as a single color
    let single = color.parse::<CssColor>();

    if let Ok(color) = single {
        let color = RGB::new(color.r, color.g, color.b);
        return Ok(std::iter::repeat(color).take(n).collect())
    }

    // Then parse as palette
    let palette = color.parse::<Palette>()
        .map_err(|_| {
            anyhow!("\"{}\" can not be parsed as color or palette", color)
        })?;

    let n_u32 = u32::try_from(n).unwrap(); 

    get_color_ramp(palette, n_u32)
        .ok_or(anyhow!("Invalid number of colors ({}) for \"{}\"", n, color))

}


//=============================================================================
// Plot1D


#[derive(Debug, Default, Getters, MutGetters, Setters)]
#[getset(get = "pub")]
pub struct Plot1D {

    #[getset(set = "pub")]
    name: String,

    #[getset(set = "pub")]
    experiments: Vec<String>,
    #[getset(set = "pub")]
    species: Vec<String>,

    #[getset(set = "pub")]
    frequency: f64,
    #[getset(set = "pub")]
    ppm: [f64; 2],
    #[getset(set = "pub")]
    delta: f64,

    #[getset(set = "pub")]
    apodization: Apodization,

    #[getset(set = "pub")]
    split_experiments: bool,

    #[getset(set = "pub")]
    plot_data: bool,
    #[getset(set = "pub")]
    plot_residual: bool,
    #[getset(set = "pub")]
    plot_baseline: bool,
    #[getset(set = "pub")]
    plot_peaks: bool,

    //#[getset(get)]
    //group_peaks: Scaffold,
    #[getset(set = "pub")]
    component: Component,

    #[getset(set = "pub")]
    data_width: f64,
    #[getset(set = "pub")]
    residual_width: f64,
    #[getset(set = "pub")]
    baseline_width: f64,
    #[getset(set = "pub")]
    peaks_width: f64,

    data_colors: Vec<RGB<u8>>,
    residual_colors: Vec<RGB<u8>>,
    baseline_colors: Vec<RGB<u8>>,
    peaks_colors: Vec<RGB<u8>>,

    default_parameters: Parameters,
}


//=============================================================================
// YAML parsing

impl Plot1D {

    //-------------------------------------------------------------------------
    pub fn from_template(
        template: Plot1DTemplate, 
        job: &Job1D
    ) -> Result<Self> {

        // Starting with name
        let name = template.plot;

        // Ensuring that species and experiments exist
        let species = template.species.patterns()?;
        job.check_species_globs(&species)?;

        let species = job.get_species_from_globs(&species);

        let experiments = template.experiments.patterns()?;
        job.check_experiment_globs(&experiments)?;

        let experiments = job.get_experiments_from_globs(&experiments);

        // Parsing colours
        let color_names = vec![
            template.data_color.clone(), 
            template.residual_color.clone(), 
            template.baseline_color.clone(), 
            template.peaks_color.clone()
        ];

        let mut n = experiments.len();
        if n == 0 {
            n = species.len();
        }

        let mut colors = color_names.into_iter()
            .map(|x| {
                parse_color(x, n)
            })
            .collect::<Result< Vec<Vec<RGB<u8>>> >>()?;

        // Storing length for later before moving
        let n_experiments = experiments.len();

        // Override toggles if there are no species
        let mut plot_peaks = template.plot_peaks;
        let mut plot_residual = template.plot_residual;
        let mut plot_baseline = template.plot_baseline;
        let plot_data = template.plot_data;

        if species.len() == 0 {
            plot_peaks = false;
            plot_residual = false;
            plot_baseline = false;
        }

        // Initialize general plot
        let mut plot = cascade!{
            Plot1D::new(species, experiments);
                ..set_name(name);
                ..set_ppm(template.ppm);
                ..set_delta(template.delta);

                ..set_apodization(template.apodization);

                ..set_split_experiments(template.split_experiments);

                ..set_component(template.component);

                ..set_plot_peaks(plot_peaks);
                ..set_plot_residual(plot_residual);
                ..set_plot_baseline(plot_baseline);
                ..set_plot_data(plot_data);
                
                ..set_peaks_width(template.peaks_width);
                ..set_residual_width(template.residual_width);
                ..set_baseline_width(template.baseline_width);
                ..set_data_width(template.data_width);
                
                ..set_peaks_colors(colors.pop().unwrap())?;
                ..set_residual_colors(colors.pop().unwrap())?;
                ..set_baseline_colors(colors.pop().unwrap())?;                  
                ..set_data_colors(colors.pop().unwrap())
                    .expect("Error setting data colors.");
        };

        // Ensure frequency if there are no experiments
        if n_experiments == 0 {

            if let Some(frequency) = template.frequency {
                plot.set_frequency(frequency);
            } else {
                bail!(
                    "Default frequency must be set if no experiments defined."
                )
            }
        }

        Ok( plot )
    }
}


//=============================================================================
// Constructors

impl Plot1D {

    //--------------------------------------------------------------------------
    pub fn new(
        species: Vec<String>,
        experiments: Vec<String>, 
    ) -> Plot1D {

        // Setting sensible global defaults
        let default_parameters = Parameters {
            position: f64::NAN,
            width: 1.0,
            height: 0.5,
            fraction: 0.0,
        };

        Plot1D {
            species: species,
            experiments: experiments,
            default_parameters: default_parameters,
            ..Default::default()
        }
    }
}


//==============================================================================
// Compound getters/setters

impl Plot1D {

    //--------------------------------------------------------------------------
    pub fn set_data_colors(
        &mut self, 
        colors: Vec<RGB<u8>>
    ) -> Result<()> {
        
        if self.experiments.len() > 0 {
            if colors.len() != self.experiments.len() {
                bail!("Length of color vector must match experiments.")
            }
        } else {
            if colors.len() != self.species.len() {
                bail!("Length of color vector must match species.")
            }
        }

        self.data_colors = colors;
        Ok( () )
    }

    //--------------------------------------------------------------------------
    pub fn set_residual_colors(
        &mut self, 
        colors: Vec<RGB<u8>>
    ) -> Result<()> {
        
        if self.experiments.len() > 0 {
            if colors.len() != self.experiments.len() {
                bail!("Length of color vector must match experiments.")
            }
        }

        self.residual_colors = colors;
        Ok( () )
    }

    //--------------------------------------------------------------------------
    pub fn set_baseline_colors(
        &mut self, 
        colors: Vec<RGB<u8>>
    ) -> Result<()> {
        
        if self.experiments.len() > 0 {
            if colors.len() != self.experiments.len() {
                bail!("Length of color vector must match experiments.")
            }
        }

        self.baseline_colors = colors;
        Ok( () )
    }

    //--------------------------------------------------------------------------
    pub fn set_peaks_colors(
        &mut self, 
        colors: Vec<RGB<u8>>
    ) -> Result<()> {
        
        if self.experiments.len() > 0 {
            if colors.len() != self.experiments.len() {
                bail!("Length of color vector must match experiments.")
            }
        }

        self.peaks_colors = colors;
        Ok( () )
    }
}


//==============================================================================
// Plotting

impl Plot1D {

    //--------------------------------------------------------------------------
    pub fn plot(&self, job: &Job1D) -> Vec<Plot> {

        if self.experiments.len() == 0 {
            vec![self.plot_species(job)]
        }
        else if self.split_experiments {
            vec![]
            //self.plot_experiments_split(job)
        } else {
            vec![self.plot_experiments(job)]
        }
    }

    //--------------------------------------------------------------------------
    pub fn plot_species(&self, job: &Job1D) -> Plot {

        let f = match self.component {
            Component::R =>|x: &Complex<f64>| x.re,
            Component::I =>|x: &Complex<f64>| x.im,
        };

        // Initializing plot object
        let mut plot = new_plot(&self.name, self.ppm[1], self.ppm[0]);

        let n = ((self.ppm[1] - self.ppm[0])/self.delta).floor();
        let x: Vec<f64> = ( 0..(n as usize) )
            .map(|x| self.ppm[0] + (x as f64)/n*(self.ppm[1] - self.ppm[0]))
            .collect();

        for (i, species) in self.species.iter().enumerate() {

            let species = job.species()
                .get(species)
                .expect("Missing species.")
                .clone();

            let y = species 
                .species_lineshape(&x[..], self.frequency)
                .iter()
                .map(f)
                .collect::<Vec<f64>>();

            if self.plot_peaks {
                add_line(
                    &mut plot, &x, &y, 
                    &format!("{} ({})", species.name(), &self.name),
                    self.peaks_width, 
                    &self.peaks_colors[i],
                    &species.name()
                );
            }
        }

        plot
    }


    //-------------------------------------------------------------------------
    pub fn plot_experiments(&self, job: &Job1D) -> Plot {

        let f = match self.component {
            Component::R =>|x: &Complex<f64>| x.re,
            Component::I =>|x: &Complex<f64>| x.im,
        };

        // Max and min x values have to be found to set range
        let (xmin, xmax) = self.experiments.iter()
            .map(|x| {
                job.experiments()
                    .get(x)
                    .expect("Missing experiment.")
                    .data(Some(self.ppm), true)
                    .0
            })
            .flatten()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |acc, x| {
                (acc.0.min(x), acc.1.max(x))
            });

        // Initializing plot object
        let mut plot = new_plot(&self.name, xmin, xmax);

        for (i, experiment) in self.experiments.iter().enumerate() {

            // Experiment has to be cloned to enable parameter updating
            let experiment = job.experiments()
                .get(experiment)
                .expect("Missing experiment.")
                .clone();

            // TODO: Figure out what to do with missing parameters

            // For now, assume phase is always interpolated
            let (x, mut y) = experiment.data(Some(self.ppm), true);

            // Apodization
            self.apodization.apodize(&mut y);

            let y: Vec<f64> = y.iter()
                .map(f)
                .collect();
            
            let y_baseline: Vec<f64> = experiment.baseline().baseline(&x)
                .iter()
                .map(f)
                .collect();
            
            let mut y_residual: Vec<f64> = y.iter()
                .zip(y_baseline.iter())
                .map(|(y1, y2)| *y1 - *y2)
                .collect();

            let name = experiment.name();

            // Adding experimental data
            if self.plot_data {
                add_line(
                    &mut plot, &x, &y, 
                    &format!("Raw data ({})", &name),
                    self.data_width, 
                    &self.data_colors[i],
                    &name
                );
            }

            // Adding peaks
            for species in self.species.iter() {

                let mut y_fit: Vec<f64> = experiment.species_lineshape(species, &x)
                    .iter()
                    .map(f)
                    .collect();

                for (y1, y2) in y_residual.iter_mut().zip(y_fit.iter()) {
                    *y1 = *y1 - *y2;
                }

                for (y1, y2) in y_fit.iter_mut().zip(y_baseline.iter()) {
                    *y1 = *y1 + *y2;
                }

                if self.plot_peaks {
                    add_line(
                        &mut plot, &x, &y_fit, 
                        &format!("{} ({})", species, &name),
                        self.peaks_width, 
                        &self.peaks_colors[i],
                        &name
                    );
                }
            }

            if self.plot_residual {
                add_line(
                    &mut plot, &x, &y_residual, 
                    &format!("Residual ({})", &name),
                    self.residual_width, 
                    &self.residual_colors[i],
                    &name
                );
            }

            // Adding baseline
            if self.plot_baseline {
                add_line(
                    &mut plot, &x, &y_baseline, 
                    &format!("Baseline ({})", &name),
                    self.baseline_width, 
                    &self.baseline_colors[i],
                    &name
                );
            }
        }

        plot
    }


    //--------------------------------------------------------------------------
    /*
    pub fn plot_results_split(&mut self) -> Vec<Plot> {

        let f = match self.component {
            Component::R =>|x: &Complex<f64>| x.re,
            Component::I =>|x: &Complex<f64>| x.im,
        };

        let mut plots: Vec<Plot> = Vec::new();

        for (i, experiment) in self.experiments.iter().enumerate() {

            // Initializing plot object
            let title = format!("{} - {}", &self.name, &experiment.name());
            let mut plot = new_plot(&title);

            let (x, y) = experiment.get_data();

            let y_baseline: Vec<f64> = x.iter()
                .into_iter()
                .zip(self.baselines[i].iter())
                .filter(|(&x, _)| (x > self.ppm.0) && (x < self.ppm.1))
                .map(|(_, y)| f(y))
                .collect();

            let (x, y): (Vec<f64>, Vec<Complex<f64>>) = x.into_iter()
                .zip(y.into_iter())
                .filter(|(x, _)| (x > &self.ppm.0) && (x < &self.ppm.1))
                .unzip();

            let y = y.iter().map(f).collect::<Vec<f64>>();
            let mut y_residual = y.clone();

            // Adding experimental data
            if self.plot_data {
                add_line(
                    &mut plot, &x, &y, 
                    "Raw data",
                    self.data_width, 
                    &self.data_colors[i],
                    "Raw data",
                );
            }

            // Adding results
            for result in self.results[i].iter_mut() {

                let y_fit = result
                    .species_lineshape(
                        &x[..],
                        *experiment.frequency(),
                        Some(&self.default_parameters),
                    )
                    .iter()
                    .map(f)
                    .collect::<Vec<f64>>();

                for (y1, y2) in y_residual.iter_mut().zip(y_fit.iter()) {
                    *y1 = *y1 - *y2;
                }

                if self.plot_peaks {
                    add_line(
                        &mut plot, &x, &y_fit, 
                        result.name(),
                        self.peaks_width, 
                        &self.peaks_colors[i],
                        result.name(),
                    );
                }
            }

            if self.plot_residual {
                add_line(
                    &mut plot, &x, &y_residual, 
                    "Residual",
                    self.residual_width, 
                    &self.residual_colors[i],
                    "Residual"
                );
            }

            // Adding baseline
            if self.plot_baseline {
                add_line(
                    &mut plot, &x, &y_baseline, 
                    "Baseline",
                    self.baseline_width, 
                    &self.baseline_colors[i],
                   "Baseline" 
                );
            }

            plots.push(plot);
        }

        plots
    }
*/
}


//==============================================================================
// Helper functions


fn new_plot(title: &str, xmin: f64, xmax: f64) -> Plot {
    
    let mut plot = Plot::new();

    let axis = Axis::new()
        .range(vec![xmax, xmin]);

    let layout = Layout::new()
        .title(Title::new(title))
        .x_axis(axis);

    plot.set_layout(layout);

    plot
}


fn add_line(
    plot: &mut Plot, 
    x: &Vec<f64>,
    y: &Vec<f64>,
    title: &str,
    width: f64,
    color: &RGB<u8>,
    group: &str
) {

    let color: (u8, u8, u8) = color.clone().into();
    let color = Rgb::new(color.0, color.1, color.2);
    
    let line = Scatter::new(x.clone(), y.clone())
        .mode(Mode::Lines)
        .name(title)
        .line(
            Line::new()
                .width(width)
                .color(color)
        )
        .legend_group(group);

    plot.add_trace(line);

}


//==============================================================================
// YAML representation

#[derive(Clone, Debug, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct Plot1DTemplate {

    pub plot: String,

    #[serde(default)] 
    pub experiments: Globs,
    #[serde(default)]
    pub species: Globs,

    #[serde(default)] 
    pub component: Component,
    #[serde(default)]
    pub frequency: Option<f64>,
    #[serde(default = "default_ppm")]
    pub ppm: [f64; 2],
    #[serde(default = "default_delta")]
    pub delta: f64,

    #[serde(default)] 
    pub apodization: Apodization,

    #[serde(default = "default_false")]
    pub split_experiments: bool,

    #[serde(default = "default_true")]
    pub plot_data: bool,
    #[serde(default = "default_true")]
    pub plot_residual: bool,
    #[serde(default = "default_true")]
    pub plot_baseline: bool,
    #[serde(default = "default_true")]
    pub plot_peaks: bool,

    //#[serde(default)]
    //group_peaks: Scaffold,

    #[serde(default = "default_thin")]
    pub data_width: f64,
    #[serde(default = "default_thin")]
    pub residual_width: f64,
    #[serde(default = "default_thin")]
    pub baseline_width: f64,
    #[serde(default = "default_thick")]
    pub peaks_width: f64,

    #[serde(default = "default_data_color")]
    pub data_color: String,
    #[serde(default = "default_residual_color")]
    pub residual_color: String,
    #[serde(default = "default_baseline_color")]
    pub baseline_color: String,
    #[serde(default = "default_peaks_color")]
    pub peaks_color: String,
}


impl Default for Plot1DTemplate {

    fn default() -> Plot1DTemplate {

        Plot1DTemplate {

            plot: "".to_string(),

            experiments: Globs::default(),
            species: Globs::default(),

            frequency: None,
            ppm: default_ppm(),
            delta: default_delta(),
            
            apodization: Apodization::default(),

            split_experiments: false,

            plot_data: true,
            plot_residual: true,
            plot_baseline: true,
            plot_peaks: true,

            //#[serde(default)]
            //group_peaks: Scaffold,
            component: Component::default(),

            data_width: default_thin(),
            residual_width: default_thin(),
            baseline_width: default_thin(),
            peaks_width: default_thick(),

            data_color: default_data_color(),
            residual_color: default_residual_color(),
            baseline_color: default_baseline_color(),
            peaks_color: default_peaks_color(),
        }
    }
}


//-----------------------------------------------------------------------------
fn default_ppm() -> [f64; 2] {
    [f64::NEG_INFINITY, f64::INFINITY]
}

fn default_delta() -> f64 {
    0.001
}

fn default_true() -> bool {
    true
}

fn default_false() -> bool {
    false
}

fn default_thin() -> f64 {
    1.0
}

fn default_thick() -> f64 {
    2.0
}

fn default_data_color() -> String {
    "darkgrey".to_string()
}

fn default_residual_color() -> String {
    "darkred".to_string()
}

fn default_baseline_color() -> String {
    "darkcyan".to_string()
}

fn default_peaks_color() -> String {
    "cornflowerblue".to_string()
}
