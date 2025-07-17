//! Pulstar program is a binary that rasterizes a star and produces a polars DataFrame that contains
//! the (linear) variations on surface temperature, log g, and also the pulsation velocity components.
//! for each of the surface cells. 
use serde::Deserialize;
pub mod utils;


/// This structure is necessary for starting the program. 
/// It contains 
#[derive(Deserialize,Debug,PartialEq)]
pub struct pulstar_config{
    pub mode_data:Vec<PulsationMode>,
    pub star_data:StarData,
    pub time_points:Vec<f64>,
}

#[derive(Deserialize,Debug,PartialEq)]
pub struct PulsationMode{
    pub l:f64, 
    pub m: f64,
    pub rel_dr: f64,
    pub k: f64,
    pub frequency: f64,
    pub phase_offset: f64,
    pub rel_dTemp: f64,
    pub phase_rel_dTemp: f64,
    pub rel_dg: f64,
    pub  phase_rel_dg: f64,
}

#[derive(Deserialize,Debug,PartialEq)]
pub struct StarData{
    pub mass: f64,
    pub radius: f64,
    pub effective_temperature: f64,
    pub v_sini: f64,
    pub inclination_angle: f64,
}
