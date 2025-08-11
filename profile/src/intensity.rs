
use polars::{error::ErrString, prelude::*};
use crate::intensity::parse_intensity_grids::IntensityDataFrames;

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

    /// This function unwrapps the collection of the "grid file identification data", 
    /// producing
    /// * an ordered vector with all of the temperatures,
    /// * an ordered vector with all of the log g,
    /// * an ordered vector with all of the file names.
    fn unwrap_grids(&self)->UnwrappedGrids{
        let mut temperature_vector:Vec<f64>=Vec::new();
        let mut logg_vector:Vec<f64>=Vec::new();
        let mut filename_vector:Vec<String>=Vec::new();
        for grid in self.intensity_grids.iter(){
            temperature_vector.push(grid.temperature);
            logg_vector.push(grid.log_gravity);
            filename_vector.push(
                format!("{}{}",
                self.path_to_grids,
                grid.filename.clone())
            );
        }
        UnwrappedGrids { temperatures: temperature_vector, 
            log_g: logg_vector, 
            filenames: filename_vector}
    }
    /// This function creates a Data frame out of the intensity grid collection stored
    /// in the profile_config structure. 
    /// 
    /// The columns have the following headers
    /// 
    /// `temperature`| `log_gravity`| `file name`
    fn get_intensity_grids_dataframe(&self)->DataFrame{
        self.unwrap_grids().create_dataframe_sorted()
    }

    /// This function returns an instance of [IntensityDataframes] that will hold all of the relevant information parsed from the intensity grid files. 
    /// 
    /// ### Arguments:
    /// * `wavelengths` - a &[[f64]] slice that holds the requested wavelengths
    /// * `maxval_rel_dopplershift` - the maximum value of the relative Doppler shift
    /// * `minval_rel_dopplershift` - the minimum value of the relative Doppler shift
    /// ### Returns:
    /// This function returns an instance of [IntensityDataFrames] that has 
    /// * a [Vec] with the temperatures associated to each dataframe
    /// * a [Vec] with the log g values associated to each dataframe
    /// * a [Vec] with the DataFrames produced with of each relevant intensity grid file
    pub fn get_filtered_intensity_dataframes(
        &self,
        wavelengths: &[f64],
        maxval_rel_dopplershift:f64,
        minval_rel_dopplershift:f64,
    )->IntensityDataFrames{
        
        parse_intensity_grids::parse_relevant_intensity_grids(
            self.get_intensity_grids_dataframe(),
            wavelengths,
            maxval_rel_dopplershift,
            minval_rel_dopplershift
        )
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
        let file_name = self.filename.clone();
        let full_name = format!("{}{}",String::from(path), file_name);
        File::open(&full_name)
    }


}

/// This struct is useful to create a the database that handles all of the intensity grid files.
/// The data base produced using this structure will be easier to query from.
struct UnwrappedGrids{
    pub temperatures:Vec<f64>,
    pub log_g:Vec<f64>,
    pub filenames:Vec<String>,
}

impl UnwrappedGrids{
    /// This function creates a data frame with the columns "temperature","log_gravity","file name"
    /// that are Float64, Float64, &[str] in that order. 
    /// 
    /// The data frame is ordered from lower to greater both in temperature and gravity. 
    fn create_dataframe_sorted(self)->DataFrame{
        let filenames_str:Vec<&str> = self.filenames.iter().map(|s| s.as_str()).collect();

        let temperature_series = Series::new("temperature".into(),self.temperatures);
        let log_gravity_serie=Series::new("log_gravity".into(),self.log_g);
        let file_name_serie = Series::new("file name".into(),filenames_str);
        
        let df = DataFrame::new(vec![temperature_series.into(),
            log_gravity_serie.into(),
            file_name_serie.into()]).unwrap()
            .sort(["temperature","log_gravity"], 
            SortMultipleOptions::default()
            .with_order_descending_multi([false,false]));
 
       df.unwrap()
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