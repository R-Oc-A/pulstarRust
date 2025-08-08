//! Pulstar program is a binary that rasterizes a star and produces a polars DataFrame that contains
//! the (linear) variations on surface temperature, log g, and also the pulsation velocity components.
//! for each of the surface cells. 
use serde::Deserialize;
use temp_name_lib::math_module::spherical_harmonics;
use temp_name_lib::utils::{MathErrors,MACHINE_PRECISION};
use temp_name_lib::type_def::PI;
use std::fs;
use nalgebra as na;


/// This structure is necessary for starting the program. 
/// It contains `mode_data` which is a [Vec] collection of the pulsation modes to be implemented, the `star_data` that characterizes the star, and the `time points` to be simulated. 
/// 
/// Each pulsation mode contains: 
/// * `l` - degree of the mode.
/// * `m` - azimuthal order
/// * `rel_dr` - relative displacement, that is dr/r0
/// * `k` - correction factor
/// * `frequency` - this is the frequency of oscillation in cycles per day
/// * `phase offset` - this quantity is in rad
/// * `rel_dtemp` - relative temperature diference, that is dT/T0
/// * `phase_rel_dtemp` - 
/// * `rel_dg` - relative gravity diference, that is dg/g0
/// * `phase_rel_dg` -
/// 
/// The `star_data` contains
/// * `mass` - the mass of the star
/// * `radius` - the radius of the star
/// * `effective_temperature` - the effective temperature of the star.
/// * `v_sin_i` - the equatorial rotational velocity.
/// * `inclination angle` - the inclination angle respective to the observer. 
/// 
/// The `time_points` contains
/// * a vector of all of the oscillation phases to be created. 
/// It could be provided as an `Explicit` collection where the individual terms are posted explicitly 
/// or as a `Uniform` collection where the array is characterized by a beggining, an end, and the number of time points. 
#[derive(Deserialize,Debug,PartialEq)]
pub struct PulstarConfig{
    /// A vector collection of all of the modes that will be analyzed
    pub mode_data:Vec<PulsationMode>,

    /// The parameters that describe the star
    pub star_data:StarData,

    /// A vector collection of all the time points to be analized in the range [0,1]
    pub time_points:TimeType,

    /// This structure indicates that the star analysis will be performed on a spherical surface parameterized by the colatitude and the azimuthal angles on a regular grid with spacing Δθ Δφ
    pub mesh: MeshConfig,
}

/// This structure parameterizes a pulsation mode
#[derive(Deserialize,Debug,PartialEq)]
pub struct PulsationMode{
    /// The degree of the mode
    pub l: u16, 

    /// The azimuthal order of the mode
    pub m: i16,

    /// The relative radial displacement Δr/r_0
    pub rel_dr: f64,
    
    /// The correction factor k
    pub k: f64,

    /// The frequency of oscillation in cycles per day
    pub frequency: f64,

    /// The phase offset
    pub phase_offset: f64,

    /// The relative temperature difference ΔT/T_0
    pub rel_dtemp: f64,

    /// The phase offset of temperature
    pub phase_rel_dtemp: f64,

    /// The relative gravity difference Δg/g0
    pub rel_dg: f64,

    /// The phase offset of gravity
    pub phase_rel_dg: f64,

    /// Current phase of the displacement pulsation
    pub phase: f64,

    /// Current phase of the Temperature variation
    pub phase_temp:f64,

    ///Current phase of the log g variation
    pub phase_logg:f64,
}   

/// This structure parameterizes the star
#[derive(Deserialize,Debug,PartialEq)]
pub struct StarData{
    /// The mass of the star in solar units
    pub mass: f64,

    /// The radius of the star in solar units
    pub radius: f64,

    /// The effective temperature in K
    pub effective_temperature: f64,

    /// The equatorial rotational velocity
    pub v_omega: f64,

    /// The inclination angle in degrees
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
    //[Ricardo:]Here maybe some other geometries may rise
}

//----------------------------------------
//----Parsing the pulstar_input.toml------
//----------------------------------------

impl PulstarConfig {
    /// This function is used to fill the parameters required for the pulstar program to run out of the toml configuration file.
    /// #### Arguments:
    /// * `path_to_file` - this is a string that indicates the path to the `profile_input.toml` file
    /// #### Returns:
    /// * new instance of the profile config structure.
    pub fn read_from_toml(path_to_file:&str)->Self{
        let contents = match fs::read_to_string(path_to_file){
            Ok(c)=>c,
             Err(_) => { panic!("Could not read file {}",path_to_file)}
            };
        let params: PulstarConfig = match toml::from_str(&contents){
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

/// This module contains the functions and methods used for input/output
pub mod utils;

/// This module contains the functions, methods, and structures that are used for describing the 
/// geometry of the star, such as the aproximate deformation of a surface cell due to pulsations, the 
/// lagrangian displacement. It also contains methods to change between cartesian and spherical coordinates.
pub mod reference_frames;


/// This module contains the functions and methods used to compute the normal component of the pulsation velocity 
/// over a specific surface cell (e.g. coordinates (θ,φ) and size Δθ×Δφ).
pub mod local_pulsation_velocity;

/// This module contains the functions and methods used to compute the variation of temperature and gravity 
/// over  a specific surface cell (e.g. coordinates (θ,φ) and size Δθ×Δφ).
pub mod local_temperature_and_gravity;

pub trait ConvertToRad {
    fn convert_to_radians(&mut self);
}

impl ConvertToRad for PulsationMode{
    fn convert_to_radians(&mut self) {
        self.phase_rel_dg= self.phase_rel_dg.to_radians();
        self.phase_rel_dtemp = self.phase_rel_dtemp.to_radians();
    }
}

impl ConvertToRad for PulstarConfig{
    fn convert_to_radians(&mut self) {
        for mode in self.mode_data.iter_mut(){
            mode.convert_to_radians();
        }
    }
}


pub trait AdvanceInTime {
    fn advance_in_time(&mut self,time_point:f64);
}

impl AdvanceInTime for PulsationMode{
    fn advance_in_time(&mut self,time_point:f64) {
        self.phase = 2.0 * PI *(self.frequency * time_point 
            + self.phase_offset);

        self.phase_temp = self.phase + self.phase_rel_dtemp;
        
        self.phase_logg = self.phase + self.phase_rel_dg;

    }
}
impl AdvanceInTime for PulstarConfig{
    /// For each of the pulsation modes this method computes the current phase of pulsation
    /// 
    /// ### Arguments: 
    /// * `time_point` - a [f64] value that will be used to compute the phase of the pulsation
    /// ### Returns:
    /// * This function updates the phase parameter of the [PulsationMode]s.
    fn advance_in_time(&mut self,time_point:f64) {
        for mode in self.mode_data.iter_mut(){
            mode.advance_in_time(time_point);
        }
    }
}