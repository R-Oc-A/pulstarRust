
use polars::{error::ErrString, prelude::*};
use super::*;
use std::fs::File;

/// This module contains the methods and functions to process the intensity grids. 
/// It has three parts.
/// 
/// The first part consist of methods to produce a [DataFrame] that stores the file name of the relevant
/// intensity grid files
/// 
/// The second part consist of methods to obtain a [DataFrame]s from the contents of the intensity grid files 
/// 
/// The third part has the methods used to filter and compute intensity and continuum energy fluxes out the intensity [DataFrame]s
/// 
/// the Energy fluxes. 
/// 
/// The implementations for the profile config structure defined in this module are used to 
/// extract a DataFrame that works as a database for the intensity grid files. 
/// 
/// With this approach searching for the appropriate intensity grid files are carried out faster. 
impl ProfileConfig{
    /// This function checks within the collection of intensity grid files provided by the toml file
    /// to see if all of the grid files are loaded into the directory.
    /// If they are not, it should return a ErrorKind::NotFound
    pub fn intensity_grids_are_loaded(&self)->Result<(),std::io::Error>{
        for grid in self.intensity_grids.iter(){
            grid.is_there_a_file(&self.path_to_grids)?;
        }
        Ok(())
    }

}

impl IntensityGrid{
    /// This function asses whether the purported intensity grid file is stored in a given directory.
    /// ### Argument:
    /// * path - a string that indicates the relative path to the directory of the intensity grid files.
    /// ### Returns:
    /// * Ok(File) if it's able to locate the intensity file and is able to open it.
    /// * Error - otherwise.
    fn is_there_a_file(&self,path:&str)->Result<File,std::io::Error>{
        let file_name = match self{
            Self::Joris{filename,
                log_gravity:_,
                temperature:_} => { filename.clone()}
            Self::Nadya{filename,
                metalicity:_,
                log_gravity:_,
                temperature:_} => { filename.clone()}
        };
        let full_name = format!("{}{}",String::from(path), file_name);
        File::open(&full_name)
    }
}

//----------------------------------------
//------parsing intensity grids-----------
//---------------------------------------- 
pub mod parse_intensity_grids;

//--------------------------------------------------
//-----Extracting Intensity and continuum fluxes----
//--------------------------------------------------
pub mod extract_intensity_fluxes;