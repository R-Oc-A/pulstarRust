use std::fs::File;
use std::env;

//Load and parse the yaml file

//Profile input file structure. So far..

//On the Original File it has the number of time points.
//This can be inferred from the parquet file so It's no longer necessary.
//I should ask Joris about the number of time points (phases) are expected to simulate,
//so that I can make an estimation of how should the
//memory management be carried out.
//
//Wavelength range (in nm) two f64 values
//wavelength step size (in nm) one f64 value
//

//Upper limit of the maximum velocity
// an f64 value

//List of the intensity files
//an array of the temperature of the intensity files
//an array of log g values of the intensity files. 


//This will be the overall structure of the yaml file: 
/*

pub struct Config{
    pub lambda_0:f64,
    pub lambda_f:f64,
    pub delta_lbd:f64,
    pub v_max: f64,
    //pub n_phases:u32,
}
---


*/