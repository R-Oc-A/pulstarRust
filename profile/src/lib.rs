use polars::prelude::*;
use temp_name_lib::type_def::{CLIGHT};

use crate::{intensity::{get_doppler_shifted_wavelengths, get_flux_continuum, get_temp_logg_filenames}, interpolate::interpolate};
//use std::fs::File;
//use std::sync::Arc;

pub mod intensity;
pub mod interpolate;
pub mod utils;

pub struct Config{
    pub lambda_0:f64,
    pub lambda_f:f64,
    pub delta_lbd:f64,
    pub v_max: f64,
    //pub n_phases:u32,
}


//function that appends doppler shift to the pulstar df
pub fn append_doppler_shift(df:LazyFrame)// I take ownership of the data frame since I will produce a new one and want the old one to be dropped after appending
->LazyFrame{
    let doppler_shift = lit(1.0) - col("velocity")/lit(CLIGHT)*lit(1.0e3);
    df.clone().lazy().
    with_column(doppler_shift.clone().alias("doppler shift"))
}

fn extract_column_as_vectorf64(column_name: &str,df:&DataFrame)->Vec<f64>{
    let column = df.column(column_name).unwrap();
    column.f64().unwrap().into_iter().flatten().collect()
}

fn extract_column_as_vector_string (column_name: &str, df:&DataFrame)->Vec<String>{
    let column = df.column(column_name).unwrap();
    let vec_str:Vec<&str> = column.str().unwrap().into_iter().flatten().collect();

    let vecc:Vec<String> = vec_str.iter().map(|s| s.to_string()).collect();

    vecc
}
/// This function that returns the intensity flux and the continuum flux interpolated from the intensity grids
/// for each cell.
/// INPUTS: coschi-> projection of the normal vector of the surface cell with the unit vector towards the observer.
///         temperature -> 
pub fn return_flux_for_cell_thetaphi(
    coschi:f64,
    temperature:f64,
    log_gravity:f64,
    doppler_shift:f64,
    area:f64,
    wavelengths:&[f64],
    grid_id_lf:LazyFrame
    )->(Vec<f64>,Vec<f64>){
    
	//function that selects 4 grids temperature +_ times gravity+_ ..and 
    let grids_info = get_temp_logg_filenames(grid_id_lf, temperature, log_gravity);
    let grid_temperatures = grids_info.0;
    let grid_loggs = grids_info.1;
    let grid_names = grids_info.2;

    // function that gets the I_lambda and I_continuum for each of the 4 grids
    let mut flux_collection:Vec<Vec<f64>>=Vec::new();
    let mut cont_collection:Vec<Vec<f64>> = Vec::new();
    let mut wavelength_collection:Vec<Vec<f64>> = Vec::new();
    let shifted_wavelengths= get_doppler_shifted_wavelengths(doppler_shift, wavelengths);
    for name in grid_names.into_iter(){
        let fluxes_from_grid = get_flux_continuum(name, 
            &shifted_wavelengths, 
            coschi).unwrap();
        flux_collection.push(fluxes_from_grid.0);
        cont_collection.push(fluxes_from_grid.1);
        wavelength_collection.push(fluxes_from_grid.2);
    }
	
    //function that linearly interpolates Ic and I from the 4 grids
    let fluxcont_interpolated = interpolate(&grid_temperatures,
         &grid_loggs, 
         &shifted_wavelengths, 
         &flux_collection, 
         &cont_collection, 
         &wavelength_collection, 
         temperature, 
         log_gravity);
    
    (multiply_vecf64_by_scalar(fluxcont_interpolated.0,area),//<-flux
        multiply_vecf64_by_scalar(fluxcont_interpolated.1, area))//<-continuum
}

///This function returns the index where vector[index-1]<key<vector[index]
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


fn multiply_vecf64_by_scalar(vecf64:Vec<f64>,scalar:f64)->Vec<f64>{
    //let mut vec_result:Vec<f64> = Vec::new();
    let vec_result = vecf64.into_iter().map(|s| s * scalar).collect();
    vec_result
}