use anyhow::{Context, Result, anyhow};
use chrono;
use core::iter::IntoIterator;
use fern;
use glob;
use log::{info, error};
use std::{
    fs, 
    path::PathBuf,
    time::Instant,
};
use structopt::StructOpt;

use rnmrfit::processing::Job1D;

//==============================================================================
// Command-line arguments

#[derive(Debug, StructOpt)]
#[structopt(about = "Robust NMR fitting tool (work in progress).")]
struct Opts {
    #[structopt()]
    job_files: Vec<String>,
    
    #[structopt(short = "s", long = "species")]
    species_files: Vec<String>,

    #[structopt(short = "o", long = "output", default_value = "./fit")]
    output: PathBuf,

    #[structopt(short = "t", long = "tolerance", default_value = "1e-6")]
    tol: f64,
}

//==============================================================================
// Main CLI entry point

fn main() -> Result<()> {

    let opts = Opts::from_args();

    //----------------------------------
    // Initializing output directory with log
    let output = opts.output.clone();

    fs::create_dir_all(&output)
        .with_context(|| {
            format!(
                "Could not create output directory \"{}\"",
                output.to_string_lossy()
            )
        })?;

    // Deleting old log and starting new one
    let mut log_file = output.clone();
    log_file.push("log.txt");
    fs::remove_file(&log_file).ok();

    fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{} -- {:5} -- {}",
                chrono::Local::now().format("[%H:%M:%S]"),
                record.level(),
                message,
            ))
        })
        .level(log::LevelFilter::Debug)
        .chain(std::io::stdout())
        .chain(fern::log_file(&log_file).unwrap())
        .apply()
        .unwrap();
    
    //----------------------------------
    
    info!("Identifying input files.");

    // Processing file name globs
    let job_files: Vec<PathBuf> = expand_glob(opts.job_files)
        .context("Could not parse job file names.")
        .map_err(|e| {
            error!("{:?}\n\n", e); 
            e
        })?;

    let _species_files: Vec<PathBuf> = expand_glob(opts.species_files)
        .context("Could not parse species file names.")
        .map_err(|e| {
            error!("{:?}\n\n", e); 
            e
        })?;

    //----------------------------------

    let jobs: Vec<Job1D> = job_files
        .into_iter()
        .map(|x| {
            Job1D::from_yaml_file(&x)
                .with_context(|| {
                    format!(
                        "Could not generate job from \"{}\"",
                        x.to_string_lossy()
                    )
                })
                .map_err(|e| {
                    error!("{:?}\n\n", e);
                    e
                })
        })
        .collect::<Result<Vec<_>>>()?;

    //----------------------------------
    // Once jobs created, process

    let tol = opts.tol;

    for mut job in jobs {
        let time = Instant::now();
        job.process(&output, tol)?;
        info!("Elapsed time: {:.2?}", time.elapsed());
    }

    Ok( () )
}


//------------------------------------------------------------------------------
fn expand_glob(patterns: Vec<String>) -> Result<Vec<PathBuf>> {

    // Convert patterns in globs and group errors
    let paths: Vec<_> = patterns
        .into_iter()
        .map(|x| { 
            glob::glob(&x).map_err(|e| anyhow!("{} in \"{}\"", e, &x))
        })
        .collect::<Result<Vec<_>>>()
        .context("Could not parse file pattern.")?
        .into_iter()
        .map(|x| {
            x.into_iter().map(|x| x.map_err(|e| anyhow!(e)))
        })
        .flatten()
        .collect::<Result<Vec<_>>>()
        .context("Could not process file pattern.")?;

    Ok( paths )
}
