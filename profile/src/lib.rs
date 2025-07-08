use polars::prelude::*;
use temp_name_lib::type_def::{CLIGHT};
//use std::fs::File;
//use std::sync::Arc;

pub mod intensity;
pub mod interpolate;
pub struct Config{
    pub lambda_0:f64,
    pub lambda_f:f64,
    pub delta_lbd:f64,
    pub v_max: f64,
    pub n_phases:u32,
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
//function that returns the flux for each cell
pub fn return_flux_for_cell_thetaphi(
    coschi:f64,
    temperature:f64,
    log_gravity:f64,
    doppler_shift:f64,
    area:f64,
    wavelengths:&[f64],
    grid_id_df:&DataFrame
    )->Vec<Vec<f64>>{
    
    
    
	//function that selects 4 grids temperature +_ times gravity+_ ..and 
    
    //a function that creates lazyframes for the intensity grids
    //filter values that are outside of wavelength range
    //construct a lambda array dopplershifted
    //look for a function that select lambdas that are close to the values of the lambda array
    
    //function that appends Ic and I to each data frame
	//function that linearly interpolates Ic and I from the 8 grids
    //

    vec![vec![0.0],vec![0.0]]
}


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