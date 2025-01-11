use anyhow::{Result, anyhow, bail};
use log::{error, warn};
use std::{
    fs,
    fs::File,
    io::Read,
    path::PathBuf,
};
use thiserror::Error;

mod bruker;
mod rs2d;

pub use bruker::import_bruker_1d; 
pub use rs2d::import_rs2d_processed_1d;

//==============================================================================
// Utility functions

pub fn check_mismatch(
        path: &PathBuf, 
        dir: &str
    ) -> Result<(), ImportingError> {

    let mut path = path.clone();
    path.push(dir);

    if ! path.exists() {
        return Err( ImportingError::Mismatch );
    }

    Ok( () )
}


pub fn dir_names(path: &PathBuf) -> Result<Vec<String>> {

    // The same warning message may appear multiple times
    let warning = || {
        warn!(
            "Error reading entry in \"{}\". \
             This may be the result of an unreachable path \
             or a malformed file name.", 
            path.display()
        )
    };

    let paths: Vec<String> = fs::read_dir(path)
        .map_err(|_| anyhow!("Cannot access \"{}\"", path.display()))?
        // Only keeping accessible paths
        .filter_map(|x| {
            let path = x.ok();
            
            if path.is_none() { 
                warning(); 
            }
            
            path
        })
        // Only keep dirs
        .filter(|x| {
             match x.metadata() {
                 Ok(meta) => meta.is_dir(),
                 Err(_) => {
                    warning();
                    false
                 }
             }
        })
        // Keep names, generating warnings if there are issues
        .filter_map(|x| {

            let name = x.file_name()
                .into_string()
                .ok();

            if name.is_none() {
                warning();
            }

            name
        })
        .collect();

    Ok( paths )
}

pub fn read_bytes(path: &PathBuf) -> Result<Vec<u8>> {

    let name: String = format!("{}", path.display());

    if ! path.exists() {
        bail!("Cannot access data \"{}\".", &name)
    }

    let mut file = File::open(path)
        .map_err(|_| {
            anyhow!("Cannot open data \"{}\"", &name)
        })?;

    let mut buffer: Vec<u8> = vec![];
    file.read_to_end(&mut buffer)
        .map_err(|_| {
            anyhow!("Cannot read data from \"{}\"", &name)
        })?;

    Ok( buffer )
}


//==============================================================================
// Common errors

#[derive(Error, Debug)]
pub enum ImportingError {
    #[error("Error importing data.")]
    MalformedData(#[from] anyhow::Error),

    #[error("Data does not match importing method.")]
    Mismatch
}

impl ImportingError {

    /// Malformed entries are always logged so that multiple errors
    /// can be identified at once.
    pub fn malformed(error: anyhow::Error, name: &String) -> Self { 

        error!("Experiment \"{}\" -- {:?}", name, &error);
        ImportingError::MalformedData(error)
    }
}
