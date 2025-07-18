//! Pulstar program is a binary that rasterizes a star and produces a polars DataFrame that contains
//! the (linear) variations on surface temperature, log g, and also the pulsation velocity components.
//! for each of the surface cells. 
use serde::Deserialize;
use std::fs;


/// This structure is necessary for starting the program. 
/// It contains `mode_data` which is a collection of the pulsation modes to be implemented, the `star_data` that characterizes the star, and the `time points` to be simulated. 
/// 
/// The pulsation mode contains: 
/// 
/// `l` - degree of the mode.
/// 
/// `m` - azimuthal order
/// 
/// `rel_dr` - relative displacement, that is dr/r0
/// 
/// `k` - correction factor
/// 
/// `frequency` - this is the frequency of oscillation in cycles per day
/// 
/// `phase offset` - this quantity is in rad
/// 
/// `rel_dtemp` - relative temperature diference, that is dT/T0
/// 
/// `phase_rel_dtemp` - 
/// 
/// `rel_dg` - relative gravity diference, that is dg/g0
/// 
/// `phase_rel_dg` -
/// 
/// The `star_data` contains
/// 
/// `mass` - the mass of the star
/// 
/// `radius` - the radius of the star
/// 
/// `effective_temperature` - the effective temperature of the star.
/// 
/// `v_sin_i` - the equatorial rotational velocity.
/// 
/// `inclination angle` - the inclination angle respective to the observer. 
/// 
/// 
/// The `time_points` contains
/// 
/// - a vector of all of the oscillation phases to be created. 
/// It could be provided as an `Explicit` collection where the individual terms are posted explicitly 
/// or as a `Uniform` collection where the array is characterized by a beggining, an end, and the number of time points. 
#[derive(Deserialize,Debug,PartialEq)]
pub struct PPulstarConfig{
    pub mode_data:Vec<PulsationMode>,
    pub star_data:StarData,
    pub time_points:TimeType,
    pub mesh: MeshConfig,
}

#[derive(Deserialize,Debug,PartialEq)]
pub struct PulsationMode{
    pub l: u16, 
    pub m: i32,
    pub rel_dr: f64,
    pub k: f64,
    pub frequency: f64,
    pub phase_offset: f64,
    pub rel_dtemp: f64,
    pub phase_rel_dtemp: f64,
    pub rel_dg: f64,
    pub phase_rel_dg: f64,
}

#[derive(Deserialize,Debug,PartialEq)]
pub struct StarData{
    pub mass: f64,
    pub radius: f64,
    pub effective_temperature: f64,
    pub v_sini: f64,
    pub inclination_angle: f64,
}

#[derive(Deserialize,Debug,PartialEq,Clone)]
pub enum  TimeType{
    Explicit{collection:Vec<f64>},
    Uniform{ start:f64, end:f64, step:f64}
}

#[derive(Deserialize,Debug,PartialEq)]
pub enum MeshConfig{
    Sphere{theta_step:f64,
           phi_step:f64},
    //Here maybe some other geometries may rise
}

//----------------------------------------
//----Parsing the pulstar_input.toml------
//----------------------------------------

impl PPulstarConfig {
    /// This function is used to fill the parameters required for the pulstar program to run out of the toml configuration file.
    /// #### Arguments:
    ///     `path_to_file` - this is a string that indicates the path to the `profile_input.toml` file
    /// #### Returns:
    ///      new instance of the profile config structure.
    pub fn read_from_toml(path_to_file:&str)->Self{
        let contents = match fs::read_to_string(path_to_file){
            Ok(c)=>c,
             Err(_) => { panic!("Could not read file {}",path_to_file)}
            };
        let params: PPulstarConfig = match toml::from_str(&contents){
            Ok(d) => d,
            Err(e) => {println!("Unable to load data from {}",path_to_file);
                panic!("error {}",e)}
        }; 
        params
    }

    /// This function extracts the time points from the configuration file of the pulstar code as a vector with elements of `f64` type
    pub fn get_time_points(&self)->Vec<f64>{
        match self.time_points.clone() {
            TimeType::Explicit { collection } =>{collection}
            TimeType::Uniform { start, end, step } =>{
                let mut time_vec:Vec<f64> = Vec::new();
                let steps = ((end - start)/step) as usize;
                for i in 0..=steps {
                    time_vec.push(start + i as f64 * step);
                }
                time_vec
            }
        }
    }

    /// This function extracts the mesh structure from the configuration file of the pulstar code. 
    pub fn get_mesh_structure(&self)->(f64,f64){
        match self.mesh{
            MeshConfig::Sphere { theta_step,
                 phi_step } => {(theta_step,phi_step)}
        }
    }
}


pub mod utils;
