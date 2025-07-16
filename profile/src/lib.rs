use polars::prelude::*;
use serde::Deserialize;
use temp_name_lib::type_def::{CLIGHT,N_FLUX_POINTS};//Velocity of light in m/s

use crate::{intensity::{get_flux_continuum, get_temp_logg_filenames}, interpolate::interpolate};
use std::fs;
//use std::sync::Arc;

mod intensity;
mod interpolate;
pub mod utils;

// This are the necessary parameters to run the profile program.
/// The initialization of the program requires to specify the wavelength range over which the intensity fluxes will be computed.
/// 
/// This structure also captures the maximum pulsation velocity, which is used to look within the intensity grids, and 
/// a collection of the intensity grid files that will be used for the interpolation. And the path to the intensity grids which are meant to be collected on the same folder. 
#[derive(Deserialize,Debug,PartialEq)]
pub struct ProfileConfig{
    pub wavelength_range:WavelengthRange,
    pub max_velocity:f64,
    pub path_to_grids: String,
    pub intensity_grids:Vec<IntensityGrid>,
}
/// The wave length range is defined in nanometers.
/// The start should be bigger than the end and the step should be reasonable enough
#[derive(Deserialize,Debug,PartialEq)]
pub struct WavelengthRange{
    pub start: f64,
    pub end: f64,
    pub step: f64,
}
/// The intensity grids are characterized by
/// 
/// - the file name stored as a string,
/// - the temperature in Kelvin
/// - the logarithm of the surface gravity
#[derive(Deserialize,Debug,PartialEq)]
pub struct IntensityGrid{
pub temperature: f64,
pub log_gravity: f64,
pub filename: String,
}


impl ProfileConfig {
    /// This function is used to fill the parameters required for the profile program to run out of the toml configuration file.
    /// #### Arguments:
    ///     `path_to_file` - this is a string that indicates the path to the `profile_input.toml` file
    /// #### Returns:
    ///      new instance of the profile config structure.
    pub fn read_from_toml(path_to_file:&str)->Self{
        let contents = match fs::read_to_string(path_to_file){
            Ok(c)=>c,
             Err(_) => { panic!("Could not read file {}",path_to_file)}
            };
        let params: ProfileConfig = match toml::from_str(&contents){
            Ok(d) => d,
            Err(_) => { panic!("Unable to load data from {}",path_to_file)}
        }; 

        match params.intensity_grids_are_loaded(){
            Ok(()) => { params }
            Err(e) => { panic!("Intensity files are not properly loaded into the directory, '{}' not found.
 Either load the file into the directory or check the toml file to see if there was a miss spelling",e)}
        }
    }

}

impl WavelengthRange{
    /// This method returns the wavelength vector out of  the range specified on the toml file
    /// It also checks if the step size is reasonable enough
    pub fn get_wavelength_vector(&self)->Vec<f64>{
        if self.end < self.start {panic!("Wave length range is ill defined. Please correct the toml file")}
        let capacity = ((self.end-self.start)/self.step).floor() as usize + 1usize;
        if capacity >= N_FLUX_POINTS as usize {panic!("Error, too many flux points requested.")}
	    
        let mut wavelength:Vec<f64> = Vec::with_capacity(capacity);
	    wavelength.push(self.start);
	    for i in 0..=capacity {//<- inclussive loop so wavelength[capacity]==λ_f.
	            wavelength.push( self.start + self.step * (i as f64) );
        }
        wavelength
    }
}


/// This function that inserts the relative Doppler wavelength shift to the pulstar lazyframe.
/// 
/// For a surface cell, the pulsation + rotation speed creates a doppler shift for all of the wavelengths given by 
/// 
/// delta_lambda = lambda * (1 - total_velocity/speed of light)
/// 
/// this last factor is inserted into the data frame so that it reduces the frequency of times it needs to be calculated
/// 
/// ### Arguments: 
/// 
/// `lf` - a LazyFrame created out of the rasterized_star.parquet file.
/// 
/// ### Returns:
/// 
/// a `Lazyframe` with the same headers as the `rasterized_star.parquet` with an extra column |relative shift| that holds the value
/// 
/// (1 - v/c)
pub fn insert_col_relative_dlambda(lf:LazyFrame)// I take ownership of the data frame since I will produce a new one and want the old one to be dropped after appending
->LazyFrame{
    let doppler_shift = lit(1.0) - col("velocity")/lit(CLIGHT)*lit(1.0e3);
    lf.clone().lazy().
    with_column(doppler_shift.clone().alias("relative shift"))
}

/// This function takes a polars data frame and returns all of the values from a given column that holds f64 values. 
/// 
/// ### Arguments: 
/// 
/// `column_name` - a string slice that holds the name of a column. The column should hold f64 values.
/// 
/// `df`- a polars DataFrame
/// 
/// ### Returns:
/// 
/// `Vec<f64>` - a vector that contains all of the values on the column.
fn extract_column_as_vectorf64(column_name: &str,df:&DataFrame)->Vec<f64>{
    let column = df.column(column_name).unwrap();
    column.f64().unwrap().into_iter().flatten().collect()
}

/// This function takes a polars data frame and returns all of the values from a given column that holds string values. 
/// 
/// ### Arguments: 
/// 
/// `column_name` - a string slice that holds the name of a column. The column should hold String values.
/// 
/// `df`- a polars DataFrame
/// 
/// ### Returns:
/// 
/// `Vec<String>` - a vector that contains all of the values on the column.
fn extract_column_as_vector_string (column_name: &str, df:&DataFrame)->Vec<String>{
    let column = df.column(column_name).unwrap();
    let vec_str:Vec<&str> = column.str().unwrap().into_iter().flatten().collect();

    let vecc:Vec<String> = vec_str.iter().map(|s| s.to_string()).collect();

    vecc
}


/// This function that returns the intensity flux and the continuum flux interpolated from the intensity grids
/// for each surface cell of the rasterized sphere.
/// 
/// ### Arguments:
/// 
/// `coschi` -  Projection of the normal vector of the surface cell with the unit vector towards the observer.
/// 
/// `temperature` - temperature value over the surface cell of the rasterized star (in kelvin).
/// 
/// `log_gravity` - log g value over the surface cell of the rasterized star.
/// 
/// `relative shift` - relative doppler wavelength shift, this is of course related to the velocity.
/// 
/// 'area' - area of the surface cell of the rasterized star.
/// 
/// 'wavelengths` - a borrowed vector containing the observed wavelengths.
/// 
/// `grid_id_lf` - a polars LazyFrame that contains the Database of the intensity grid files. 
/// 
/// ### Returns:
/// 
/// `(Flux,Cont)` - a tupple of two vectors computed for the surface cell, where
/// 
/// 'Flux' - is the Intensity flux vector
/// 
/// 'Cont' - is the Continuum flux vector
pub fn return_flux_for_cell_thetaphi(
    coschi:f64,
    temperature:f64,
    log_gravity:f64,
    relative_shift:f64,
    area:f64,
    wavelengths:&[f64],
    grid_id_lf:LazyFrame
    )->(Vec<f64>,Vec<f64>){
    
	//unwrapping the four relevant grid files creating  a vector that holds the temperatures, logg values and the file names.
    let grids_info = get_temp_logg_filenames(grid_id_lf, temperature, log_gravity);
    let grid_temperatures = grids_info.0;
    let grid_loggs = grids_info.1;
    let grid_names = grids_info.2;

    // Initializing the vectors that will hold the flux and continuum quantities.
    let mut flux_collection:Vec<Vec<f64>>=Vec::new();
    let mut cont_collection:Vec<Vec<f64>> = Vec::new();
    let mut wavelength_collection:Vec<Vec<f64>> = Vec::new();

    // Doppler shift the wavelengths.
    let shifted_wavelengths= get_doppler_shifted_wavelengths(relative_shift, wavelengths);

    // For each of the relevant grid files, extract the intensity and continuum values.
    for name in grid_names.into_iter(){
        let fluxes_from_grid = get_flux_continuum(name, 
            &shifted_wavelengths, 
            coschi).unwrap();
        flux_collection.push(fluxes_from_grid.0);
        cont_collection.push(fluxes_from_grid.1);
        wavelength_collection.push(fluxes_from_grid.2);
    }
	
    //function that linearly interpolates the intensity flux Ic and the continuum flux I from the 4 grids
    let fluxcont_interpolated = interpolate(&grid_temperatures,
         &grid_loggs, 
         &shifted_wavelengths, 
         &flux_collection, 
         &cont_collection, 
         &wavelength_collection, 
         temperature, 
         log_gravity);
        
    // Finally return the collected quantities of the surface.
    (multiply_vecf64_by_scalar(fluxcont_interpolated.0,area),//<-flux
        multiply_vecf64_by_scalar(fluxcont_interpolated.1, area))//<-continuum
}



/// This function applies binary search on an ordered vector to find the index where a reference value
/// would be inserted. 
/// 
/// This is useful on this program because it gives the first and last values of the vector that are less than and greather than respectively.
/// 
/// ### Arguments: 
/// 
/// `vector` - a borrowed vector of `f64` 
/// 
/// `key` - a reference value.
/// 
/// ### Returns:
/// 
/// `index` - a `usize` value witht the property that  `vector[index-1]<key<vector[index]`
fn search_geq(vector:&[f64],key:f64)-> usize {
    //non empty vector
    let first_element =vector.get(0);
    if let Some(first_value) = first_element{
        if key <= *first_value {0usize}
        else{
            let size = vector.len();
            let mut top = size - 1;
            let mut bottom = 0;
            let mut middle = bottom + (top - bottom) / 2;
            if key > vector[top]{panic!("key value is too large")}

            while bottom<top{
                if vector[middle] <key{
                    bottom = middle + 1;
                } else {
                    top = middle;
                }
                middle = bottom + (top - bottom)/2;
            }
            top
        }
    }
    else{ panic!("empty vector!")}
}

/// This function multiplies a vector of `f64` by a scalar
/// 
/// ### Arguments:
/// 
/// `vecf64` - a vector containing `f64` values
/// 
/// `scalar` - an `f64` value
/// 
/// ### Returns:
/// 
/// `vec_result` - The multiplied vector.
fn multiply_vecf64_by_scalar(vecf64:Vec<f64>,scalar:f64)->Vec<f64>{
    //let mut vec_result:Vec<f64> = Vec::new();
    let vec_result = vecf64.into_iter().map(|s| s * scalar).collect();
    vec_result
}


/// This function multiplies all of the wavelengths by the relative doppler wavelength shift
/// 
/// ### Arguments:
/// 
/// `wavelengths` - a vector containing the observed wavelengths (in nm)
/// 
/// `relative shift` - an `f64` value that holds the relative doppler shift
/// 
/// ### Returns:
/// 
/// `Vec<f64>` - The doppler shifted wavelenghts.
/// This function applies a doppler shift to the observed wavelengths. It is a multiplication by a scalar
fn get_doppler_shifted_wavelengths(relative_shift:f64,
        wavelengths:&[f64])->Vec<f64>{
            let mut shifted:Vec<f64>=Vec::new();
            for lambda in wavelengths.iter(){
                shifted.push(*lambda * relative_shift);
            }
            shifted
}