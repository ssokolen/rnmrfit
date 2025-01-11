use anyhow::{Result, anyhow};
use cascade::cascade;
use itertools_num;
use roxmltree::Document;
use rustfft::num_complex::Complex;
use std::{
    fs,
    collections::{HashSet, HashMap},
    convert::TryInto,
    path::PathBuf,
};

use crate::importing::{check_mismatch, dir_names, read_bytes, ImportingError};
use crate::processing::Experiment1D;

//------------------------------------------------------------------------------
pub fn import_rs2d_processed_1d(path: &PathBuf, processing: &Option<String>) 
-> Result<Experiment1D, ImportingError> {

    // Storing experiment for later
    let experiment = path.clone();
    let name: String = format!("{}", experiment.display());

    // First, check if folder includes metadata file
    check_mismatch(path, "Serie.xml")?;

    // Then, descend into Proc
    let mut path = path.clone();
    path.push("Proc");

    let proc_dirs: HashSet<String> = dir_names(&path)
        .map_err(|_| ImportingError::Mismatch)?
        .into_iter()
        .collect();

    if proc_dirs.len() == 0 {
        let e = anyhow!("No processed files found.");
        return Err( ImportingError::malformed(e, &name) );
    }

    // Select appropriate processing directory
    let processing = match processing {
        Some(processing) => processing,
        None => {
            let e = anyhow!("RS2D data cannot be processed manually.");
            return Err( ImportingError::malformed(e, &name) );
        }
    };

    if processing == "" {
        let dir = proc_dirs.iter().min().unwrap();
        path.push(dir);
    } else if proc_dirs.contains(processing) {
        path.push(processing);
    } else {
        let e = anyhow!("Cannot access processed data \"{}\".", processing);
        return Err( ImportingError::malformed(e, &name) );
    }

    // Load complex data, reverse and scale it
    let mut data_path = path.clone();
    data_path.push("data.dat");

    // Bytes are assumed to be big-endian f32
    let intensity: Vec<Complex<f64>> = read_bytes(&data_path)
        .map_err(|e| ImportingError::malformed(e, &name) )?
        .chunks(4)
        .map(|x| {
            f32::from_be_bytes(
                x.try_into()
                    .expect("Error reading bytes.")
            ) 
        })
        .collect::<Vec<f32>>()
        .chunks(2)
        .map(|x| Complex::<f64>::new(x[0].into(), x[1].into()) )
        .collect();

    let mut intensity: Vec<Complex<f64>> = intensity.into_iter()
        .collect();

    let n = intensity.len();

    // Parsing header xml
    let mut path_header = path.clone();
    path_header.push("header.xml");

    let header: String = fs::read_to_string(path_header)
        .map_err(|_| {
            let e = anyhow!("Error reading header.xml");
            ImportingError::malformed(e, &name)
        })?;

    let root = Document::parse(&header)
        .map_err(|_| {
            let e = anyhow!("Error parsing header.xml as an xml file.");
            ImportingError::malformed(e, &name)
        })?;

    // Extracting all the params entries rather than parsing
    let entries = root.descendants()
        .filter(|n| n.has_tag_name("entry"));

    // Settings can hold many irrelevant values so a separate map
    // is used to extract just the immediate values.
    let mut settings: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
    let mut setting_values: HashMap<String, Vec<String>> = HashMap::new();

    for node in entries {

        // Setting nodes are expected to have a key and value element,
        // where value holds other nodes
        let key = node.children()
            .find(|n| {
                n.has_tag_name("key")
            });

        let values = node.children()
            .find(|n| {
                n.has_tag_name("value")
            });

        let mut values_map: HashMap<String, Vec<String>> = HashMap::new();

        if let (Some(key), Some(values)) = (key, values) {

            let key = key.text()
                .unwrap()
                .to_string();

            for inner_node in values.children() {

                let tag = inner_node.tag_name()
                    .name()
                    .to_string();

                let value = inner_node.text()
                    .or(Some(""))
                    .unwrap()
                    .to_string();

                let map = values_map.entry(tag)
                    .or_insert(Vec::new());

                map.push(value);
            }


            let value = match values_map.get("value") {
                Some(values) => values.clone(),
                None => Vec::new()
            };

            setting_values.insert(key.clone(), value);
            settings.insert(key, values_map);
        }
    }

    let frequency = setting_values.get("OBSERVED_FREQUENCY")
        .unwrap()[0]
        .clone();
    let width = setting_values.get("SPECTRAL_WIDTH")
        .unwrap()[0]
        .clone();
    let nucleus = setting_values.get("OBSERVED_NUCLEUS")
        .unwrap()[0]
        .clone();

    // Converting frequency into MHz
    let frequency = frequency.parse::<f64>().unwrap()/1e6;

    // The width value may also be in Hz or ppm
    let mut width = width.parse::<f64>().unwrap();

    let width_type = &settings.get("SPECTRAL_WIDTH")
        .unwrap()
        .get("numberEnum")
        .unwrap()[0];

    match width_type.as_str() {
        "FrequencyOffset" => width /= frequency,
        "FrequencyShift" => {},
        // For older versions of RS2D files
        "SW" => width /= frequency,
        _ => panic!("Unrecognized offset type.")
    }

    // Offset depends on observed nucleus
    let mut channel = 0;

    for i in 1..5 {
        channel += 1;

        let key = format!("NUCLEUS_{}", i);
        let value = &setting_values.get(&key).unwrap()[0];

        if &nucleus == value {
            break
        }
    }

    let key = format!("OFFSET_FREQ_{}", channel);
    let offset = setting_values.get(&key).unwrap()[0].clone();

    // The offset value may be in Hz or ppm
    let mut offset = offset.parse::<f64>().unwrap();
    
    let offset_type = &settings.get(&key)
        .unwrap()
        .get("numberEnum")
        .unwrap()[0];
    
    match offset_type.as_str() {
        "FrequencyOffset" => offset /= frequency,
        "FrequencyShift" => {},
        _ => panic!("Unrecognized offset type.")
    }

    // Generating ppm range
    let lower = -width/2.0 + offset;
    let upper = width/2.0 + offset;

    let mut chemical_shift = itertools_num::linspace(upper, lower, n)
        .collect();

    let sampling_frequency = width*frequency;

    let experiment = Experiment1D::new(
            path,
            frequency,
            sampling_frequency,
            chemical_shift,
            intensity,
        )
        .map_err(|e| ImportingError::MalformedData(e) )?;

    Ok( experiment )
}

