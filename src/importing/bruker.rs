use anyhow::{Result, anyhow};
use cascade::cascade;
use itertools::Itertools;
use itertools_num;
use lazy_static::lazy_static;
use regex::Regex;
use rustfft::{FftPlanner, num_complex::Complex};
use std::{
    f64::consts::PI,
    fs::File,
    io::{Read},
    collections::{HashSet, HashMap},
    convert::TryInto,
    path::PathBuf,
};

use crate::importing::{check_mismatch, dir_names, read_bytes, ImportingError};
use crate::processing::Experiment1D;

//-----------------------------------------------------------------------------
// Helper function for converting strings to UTF-8
fn read_as_utf8(filename: &PathBuf) -> Result<String> {

    let mut file = File::open(filename)?;

    let mut contents = Vec::new();
    file.read_to_end(&mut contents)?;

    // Convert contents to UTF-8
    let utf8_contents: String = contents.into_iter()
        .map(|c| c as char)
        .collect();

    Ok(utf8_contents)
}


//-----------------------------------------------------------------------------
/// Import either raw FID or processed data depending on processing string.
pub fn import_bruker_1d(path: &PathBuf, processing: &Option<String>) 
        -> Result<Experiment1D, ImportingError> {

    match processing.as_ref() {
        None => import_bruker_fid_1d(path),
        Some(processing) => import_bruker_processed_1d(path, processing)
    }
}


//-----------------------------------------------------------------------------
fn import_bruker_processed_1d(path: &PathBuf, processing: &String) 
        -> Result<Experiment1D, ImportingError> {

    // Storing experiment for later
    let experiment = path.clone();
    let name: String = format!("{}", experiment.display());

    // First, check if folder includes metadata file
    check_mismatch(path, "acqus")?;

    // Then, descend into pdata 
    let mut path = path.clone();
    path.push("pdata");

    let proc_dirs: HashSet<String> = dir_names(&path)
        .map_err(|_| ImportingError::Mismatch)?
        .into_iter()
        .collect();

    if proc_dirs.len() == 0 {
        let e = anyhow!("No processed files found.");
        return Err( ImportingError::malformed(e, &name) );
    }

    // Select appropriate processing directory
    if processing == "" {
        let dir = proc_dirs.iter().min().unwrap();
        path.push(dir);
    } else if proc_dirs.contains(processing) {
        path.push(processing);
    }
    else {
        let e = anyhow!("Cannot access processed data \"{}\".", processing);
        return Err( ImportingError::malformed(e, &name) );
    }

    // Parsing procs
    let mut path_procs = path.clone();
    path_procs.push("procs");

    let procs: String = read_as_utf8(&path_procs)
        .map_err(|_| {
            let e = anyhow!("Error reading procs");
            ImportingError::malformed(e, &name)
        })?;

    // Extracting all the entries first to get at BYTORDP
    let mut settings: HashMap<String, Entry> = HashMap::new();

    for line in procs.split("\n") {
        if let Some( (key, value) ) = Entry::parse(line) {
            settings.insert(key, value);
        };
    }

    let frequency = settings.get("SF")
        .unwrap()
        .float()
        .unwrap();

    let width = settings.get("SW_p")
        .unwrap()
        .float()
        .unwrap();

    let offset = settings.get("OFFSET")
        .unwrap()
        .float()
        .unwrap();

    let byteord = settings.get("BYTORDP")
        .unwrap()
        .int()
        .unwrap();

    // Choosing conversion function based on byteord
    let f = 
        if byteord == 0 {
            i32::from_le_bytes
        } else if byteord == 1 {
            i32::from_be_bytes
        } else {
            let e = anyhow!("Unexpected BYTORDP value.");
            return Err( ImportingError::malformed(e, &name) );
        };

    // Load complex data, reverse and scale it
    let mut real_path = path.clone();
    real_path.push("1r");

    let real: Vec<i32> = read_bytes(&real_path)
        .map_err(|e| ImportingError::malformed(e, &name) )?
        .chunks(4)
        .map(|x| f(x.try_into().expect("Error reading bytes.")))
        .collect();

    let mut imag_path = path.clone();
    imag_path.push("1i");

    let imag: Vec<i32> = read_bytes(&imag_path)
        .map_err(|e| ImportingError::malformed(e, &name) )?
        .chunks(4)
        .map(|x| f(x.try_into().expect("Error reading bytes.")))
        .collect();

    // Flipping the imaginary domain sign is a bit of hack for now
    let intensity: Vec<Complex<f64>> = real.into_iter()
        .zip(imag.into_iter())
        .map(|(r, i)| Complex::<f64>::new(r.into(), (-i).into()))
        //.rev()
        .collect();

    let n = intensity.len();

    // Generating ppm range
    let lower = offset;
    let upper = offset - width/frequency;

    let chemical_shift = itertools_num::linspace(lower, upper, n)
        .collect();

    let sampling_frequency = (lower-upper) * frequency;

    let experiment = Experiment1D::new(
            experiment, 
            frequency, 
            sampling_frequency, 
            chemical_shift, 
            intensity
        )
        .map_err(|e| ImportingError::MalformedData(e) )?;

    Ok( experiment )
}


//-----------------------------------------------------------------------------
fn import_bruker_fid_1d(path: &PathBuf) 
-> Result<Experiment1D, ImportingError> {

    // Storing experiment for later
    let experiment = path.clone();
    let name: String = format!("{}", experiment.display());

    // First, check if folder includes metadata file
    check_mismatch(path, "acqus")?;

    // Parsing acqus
    let mut path_acqus = path.clone();
    path_acqus.push("acqus");

    let acqus: String = read_as_utf8(&path_acqus)
        .map_err(|_| {
            let e = anyhow!("Error reading acqus");
            ImportingError::malformed(e, &name)
        })?;

    // Extracting all the entries first to get at BYTORDA
    let mut settings: HashMap<String, Entry> = HashMap::new();

    for line in acqus.split("\n") {
        if let Some( (key, value) ) = Entry::parse(line) {
            settings.insert(key, value);
        };
    }

    let frequency = settings.get("SFO1")
        .unwrap()
        .float()
        .unwrap();

    let offset = settings.get("O1")
        .unwrap()
        .float()
        .unwrap();

    let width = settings.get("SW_h")
        .unwrap()
        .float()
        .unwrap();

    let data_type = settings.get("DTYPA")
        .unwrap()
        .int()
        .unwrap();

    let byte_order = settings.get("BYTORDA")
        .unwrap()
        .int()
        .unwrap();

    let n_expected = settings.get("TD")
        .unwrap()
        .int()
        .unwrap() as usize;

    let group_delay = settings.get("GRPDLY")
        .map(|x| x.float().unwrap() )
        .unwrap_or(0.0);
            
    // TODO: Add group_delay lookup
    /*
    let decim = settings.get("DECIM")
        .map(|x| x.int() )
        .unwrap_or(0);

    let dspvfs = settings.get("DSPFVS")
        .map(|x| x.int() )
        .uwnrap_or(0);
    */

    // Choosing byte sise based on data type
    let n_bytes = match data_type {
        0 => 4,
        2 => 8,
        _ => {
            let e = anyhow!("Unexpected DTYPA value.");
            return Err( ImportingError::malformed(e, &name) );
        }
    };

    // Choosing conversion function based on byte order
    let f = match (byte_order, data_type) {
        (0, 0) => |x: &[u8]| i32::from_le_bytes(x.try_into().unwrap()) as f64,
        (1, 0) => |x: &[u8]| i32::from_be_bytes(x.try_into().unwrap()) as f64,
        (0, 2) => |x: &[u8]| f64::from_le_bytes(x.try_into().unwrap()),
        (1, 2) => |x: &[u8]| f64::from_be_bytes(x.try_into().unwrap()),
        _ => {
            let e = anyhow!("Unexpected BYTORDA value.");
            return Err( ImportingError::malformed(e, &name) );
        }
    };

    // Load complex data and transform it
    let mut fid_path = path.clone();
    fid_path.push("fid");

    let mut fid: Vec<f64> = read_bytes(&fid_path)
        .map_err(|e| ImportingError::malformed(e, &name) )?
        .chunks(n_bytes)
        .map(f)
        .collect();

    // Remove potential padding
    if fid.len() >= n_expected {
        let diff = (fid.len() - n_expected)/2;

        for _ in 0 .. diff {
            if (fid[fid.len()-1] == 0.0) && (fid[fid.len()-2] == 0.0) {
                fid.pop();
                fid.pop();
            }
        }
    }

    if fid.len() != n_expected {
        let e = anyhow!("Unexpected fid length.");
        return Err( ImportingError::malformed(e, &name) );
    }

    let (real, imag): (Vec<f64>, Vec<f64>) = fid.into_iter()
        .tuples()
        .unzip();

    // Flipping the imaginary domain sign is a bit of hack for now
    let mut intensity: Vec<Complex<f64>> = real.into_iter()
        .zip(imag.into_iter())
        .map(|(r, i)| Complex::<f64>::new(r, i))
        .collect();

    // For debugging
    let complex = intensity.clone();

    let n = intensity.len();

    // Take FFT
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(n);
    fft.process(&mut intensity);

    // Shift negative values
    let positives: Vec<_> = intensity.drain(..((n+1)/2 + 1)).collect();
    intensity.extend(positives);

    // Phase correct
    let phase = group_delay;

    if phase > 0.0 {
        for (i, value) in intensity.iter_mut().enumerate() {
            let fraction = (i as f64)/(n as f64);
            *value *= Complex::new(0.0, 2.0 * PI * phase * fraction).exp();
        }
    }

    // Reverse
    intensity.reverse();

    // Generating ppm range
    let lower = offset/frequency + width/(2.0*frequency);
    let upper = offset/frequency - width/(2.0*frequency);

    let chemical_shift = itertools_num::linspace(lower, upper, n)
        .collect();

    let sampling_frequency = (lower-upper) * frequency;

    let experiment = Experiment1D::new(
            experiment, 
            frequency, 
            sampling_frequency, 
            chemical_shift, 
            intensity
        )
        .map_err(|e| ImportingError::MalformedData(e) )?;

    Ok( experiment )
}



//=============================================================================
// Basic proc entry enum

#[derive(Debug)]
enum Entry {
    Int(i32),
    Float(f64),
    String(String)
}


impl Entry {
    
    pub fn parse_int(line: &str) -> Option<(String, Entry)> {

        lazy_static! {
            static ref RE_INT: Regex = 
                Regex::new(r"##\$(.*?)= (-?\d+)").unwrap();
        }

        let captures = RE_INT.captures(line)?;
        let key: Option<String> = captures.get(1)
            .map(|x| x.as_str().to_string() );

        let value: Option<i32> = captures.get(2)
            .and_then(|x| x.as_str().parse().ok() );

        if key.is_some() && value.is_some() {
            Some( ( key.unwrap(), Entry::Int(value.unwrap()) ) )
        } else {
            None
        }
    }


    pub fn parse_float(line: &str) -> Option<(String, Entry)> {

        lazy_static! {
            static ref RE_FLOAT: Regex = 
                Regex::new(r"##\$(.*?)= (-?\d+\.\d+)").unwrap();
        }

        let captures = RE_FLOAT.captures(line)?;
        let key: Option<String> = captures.get(1)
            .map(|x| x.as_str().to_string() );

        let value: Option<f64> = captures.get(2)
            .and_then(|x| x.as_str().parse().ok() );

        if key.is_some() && value.is_some() {
            Some( ( key.unwrap(), Entry::Float(value.unwrap()) ) )
        } else {
            None
        }
    }


    pub fn parse_string(line: &str) -> Option<(String, Entry)> {

        lazy_static! {
            static ref RE_STRING: Regex = 
                Regex::new(r"##\$(.*?)= <(.*)>").unwrap();
        }

        let captures = RE_STRING.captures(line)?;
        let key: Option<String> = captures.get(1)
            .map(|x| x.as_str().to_string() );

        let value: Option<String> = captures.get(2)
            .map(|x| x.as_str().to_string() );

        if key.is_some() && value.is_some() {
            Some( ( key.unwrap(), Entry::String(value.unwrap()) ) )
        } else {
            None
        }
    }


    pub fn parse(line: &str) -> Option<(String, Entry)> {
        
        Entry::parse_float(line)
            .or_else(|| Entry::parse_int(line))
            .or_else(|| Entry::parse_string(line))
    }


    pub fn float(&self) -> Option<f64> {
        
        match self {
            Entry::Int(int) => Some(*int as f64),
            Entry::Float(float) => Some(*float),
            _ => None
        }
    }


    pub fn int(&self) -> Option<i32> {
        
        match self {
            Entry::Int(int) => Some(*int),
            _ => None
        }
    }
}
