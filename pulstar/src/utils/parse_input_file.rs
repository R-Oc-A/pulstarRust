use serde::Deserialize;
use crate::{MeshConfig, ParsingFromToml, PulsationMode, StarData, TimeType};
use crate::PI;
use std::fs;
#[derive(Deserialize,Debug,PartialEq)]
pub struct InputParameters{
    /// A vector collection of all of the modes that will be analyzed, No phases for velocity, temperatuer and gravity are added.
    pub mode_data:Vec<PulsationModeNoPhases>,

    /// The parameters that describe the star
    pub star_data:StarData,

    /// A vector collection of all the time points to be analized in the range [0,1]
    pub time_points:TimeType,

    /// This structure indicates that the star analysis will be performed on a spherical surface parameterized by the colatitude and the azimuthal angles on a regular grid with spacing Δθ Δφ
    pub mesh: MeshConfig,
}

/// This structure parameterizes a pulsation mode
#[derive(Deserialize,Debug,PartialEq)]
pub struct PulsationModeNoPhases{
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
}   




impl ParsingFromToml for InputParameters {
    /// This function is used to read the parameters required for the pulstar program to run out of the toml configuration file.
    /// #### Arguments:
    /// * `path_to_file` - this is a string that indicates the path to the `profile_input.toml` file
    fn read_from_toml(contents:&str)->Self {
        let params: Self = match toml::from_str(contents){
            Ok(d) => d,
            Err(e) => {println!("Unable to load data from {}",contents);
                panic!("error {}",e)}
        }; 
        params
    }
}

impl PulsationModeNoPhases{
    pub fn get_initial_phases(no_phases_vec:Vec<Self>) -> Vec<PulsationMode>{
        let mut mode_data:Vec<PulsationMode> = Vec::new();
        for mode in no_phases_vec.into_iter(){
            mode_data.push(
                PulsationMode { l: mode.l,
                    m:mode.m, 
                    rel_dr:mode.rel_dr,
                    k:mode.k,
                    frequency:mode.frequency,
                    phase_offset:mode.phase_offset,
                    rel_dtemp:mode.rel_dtemp,
                    phase_rel_dtemp:mode.phase_rel_dtemp,
                    rel_dg:mode.rel_dg,
                    phase_rel_dg:mode.phase_rel_dg,
                    phase:mode.phase_offset * 2.0 * PI,
                    phase_temp: mode.phase_offset * 2.0 *PI + mode.phase_rel_dtemp.to_radians(),
                    phase_logg: mode.phase_offset * 2.0 *PI + mode.phase_rel_dg.to_radians()}
            )
        }
        mode_data
    }
}