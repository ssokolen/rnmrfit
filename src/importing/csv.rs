use anyhow::{Result, anyhow};
use cascade::cascade;
use itertools_num;
use roxmltree::Document;
use rustfft::num_complex::Complex;
use std::{
    fs,
    collections::HashMap,
    convert::TryInto,
    fs::File,
    io::prelude::*,
    path::PathBuf,
    str::FromStr,
};

use crate::importing::{Experiment1D, ImportingError};

//------------------------------------------------------------------------------
pub fn import_rs2d_processed_1d(path: &PathBuf, processing: u32) 
-> Result<Experiment1D, ImportingError> {

    // Storing experiment for later
    let experiment_path = path.clone();
    let experiment = path.file_name().unwrap().to_str().unwrap().to_string();

    // First, check if folder includes metadata file
    let mut path_serie = path.clone();
    path_serie.push("Serie.xml");
    
    if ! path_serie.exists() {
        return Err( ImportingError::Mismatch )
    }

    // Then, descend into Proc
    let mut path = path.clone();
    path.push("Proc");

    let proc_files = fs::read_dir(&path)
        .map_err(|_| ImportingError::Mismatch)?;

    let mut proc_files = proc_files.map(|x| {
            x.ok()
                .and_then(|y| Some(y.file_name()) )
                .and_then(|y| y.into_string().ok() )
                .and_then(|y| u32::from_str(&y).ok() )

        })
        .filter(|x| x.is_some() )
        .map(|x| x.unwrap() );

    // See if selected processing number exists
    if processing > 0 {
        if proc_files.any(|x| x == processing) {
            path.push(&processing.to_string());
        } else {
            let error = anyhow!("Processing \"{}\" not found.", processing);
            return Err(ImportingError::MalformedData(error))
        }
    } else {
        let min = proc_files.min();

        match min {
            None => {
                let error = anyhow!("No processing files found.");
                return Err(ImportingError::MalformedData(error));
            },
            Some(value) => {
                path.push(&value.to_string())
            }
        }
    }

    // Open and convert to complex representation
    let mut path_data = path.clone();
    path_data.push("data.dat");

    let mut file = File::open(path_data)
        .map_err(|_| {
            let error = anyhow!("Error opening data.dat file.");
            ImportingError::MalformedData(error)
        })?;

    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)
        .map_err(|_| {
            let error = anyhow!("Error reading data.dat file into bytes.");
            ImportingError::MalformedData(error)
        })?;

    let mut intensity: Vec<Complex<f64>> = buffer.chunks(4)
        .map(|x| {
            f32::from_be_bytes(
                x.try_into().unwrap()
            ) 
        })
        .collect::<Vec<f32>>()
        .chunks(2)
        .map(|x| Complex::<f64>::new(x[0].into(), x[1].into()) )
        .collect();

    // Scale y data to a max real height of 1
    let max = intensity.iter()
        .map(|x| x.re )
        .fold(0.0/0.0, f64::max);
    
    intensity.iter_mut().for_each(|x| *x = *x/max);

    let n = intensity.len();

    // Parsing header xml
    let mut path_serie = path.clone();
    path_serie.push("header.xml");

    let header: String = fs::read_to_string(path_serie)
        .map_err(|_| {
            let error = anyhow!("Error reading header.xml");
            ImportingError::MalformedData(error)
        })?;

    let root = Document::parse(&header)
        .map_err(|_| {
            let error = anyhow!("Error parsing header.xml as an xml file.");
            ImportingError::MalformedData(error)
        })?;

    // Extracting all the params entries rather than parsing
    let entries = root.descendants()
        .filter(|n| n.has_tag_name("entry"));

    let mut settings: HashMap<String, Vec<String>> = HashMap::new();

    for node in entries {
        let key = node.descendants()
            .find(|n| {
                n.has_tag_name("key")
            })
            .unwrap()
            .text()
            .unwrap()
            .to_string();

        let values: Vec<String> = node.descendants()
            .filter(|n| {
                n.parent_element().unwrap().has_tag_name("value") &&
                n.has_tag_name("value")
            })
            .map(|n| n.text().unwrap().to_string() )
            .collect();

        settings.insert(key, values);
    }

    let frequency = settings.get("OBSERVED_FREQUENCY").unwrap()[0].clone();
    let width = settings.get("SPECTRAL_WIDTH").unwrap()[0].clone();
    let nucleus = settings.get("OBSERVED_NUCLEUS").unwrap()[0].clone();

    // Offset depends on observed nucleus
    let mut channel = 0;

    for i in 1..5 {
        channel += 1;

        let key = format!("NUCLEUS_{}", i);
        let value = &settings.get(&key).unwrap()[0];

        if &nucleus == value {
            break
        }
    }

    let key = format!("OFFSET_FREQ_{}", channel);
    let offset = settings.get(&key).unwrap()[0].clone();

    // Converting frequency into MHz
    let frequency = frequency.parse::<f64>().unwrap()/1e6;
    let offset = offset.parse::<f64>().unwrap();
    let width = width.parse::<f64>().unwrap();

    // Generating ppm range
    let lower = (-width/2.0 + offset)/frequency;
    let upper = (width/2.0 + offset)/frequency;

    let chemical_shift = itertools_num::linspace(lower, upper, n)
        .collect();

    let experiment = cascade!{
        Experiment1D::new(chemical_shift, intensity)
                .map_err(|e| ImportingError::MalformedData(e) )?;
            ..set_name(experiment);
            ..set_path(experiment_path);
            ..set_frequency(frequency);
    };

    Ok( experiment )
}
