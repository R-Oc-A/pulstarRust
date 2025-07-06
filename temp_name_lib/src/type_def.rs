use nalgebra as na;
use std::f64::consts::PI as PI;

pub const GRAVCONSTANT:f64 = 6.67259e-11;       // SI-units (m^3/s^2/kg) */
pub const MASSSUN:f64 = 1.9891e30;              // SI-units (kg)         */
pub const RADIUSSUN:f64 = 6.9599e8;             // SI-units (m)          */
pub const CLIGHT:f64 = 299792458.0;             // SI-units (m/s)        */

//pub const PI:f64 = 3.14159265358979;
pub const DEG2RAD:f64 = PI/180.0;             // Conversion from degrees to radians */
pub const RAD2DEG:f64 = 180.0/PI;             // Conversion from radians to degrees */
pub const SEC_IN_DAY:f64 = 86400.0;             // Number of seconds in a day         */
pub const CYCLI2RAD:f64 = 2.0*PI/SEC_IN_DAY;  // Conversion from cycli/day to rad/s */


pub const THETA_STEP:u16 = 4; //step in grid in the theta direction
pub const PHI_STEP:u16 = 8; // step in grid in the phi direction
pub const N_FLUX_POINTS:u16=10000; // number of points in one flux profile
pub const MAX_N_MODES:u16 = 10; // max number of pulsation modes allowed
pub const MAX_N_TIMES:u16 = 3000; //max number of time points 

//[Ricardo]: components of a 3d vector
//x,y,z cartesian corresponds to coords[1], coords[2], coords[3]
//r,theta,phi spherical corresponds to coords[1], coords[2], coords[3]
#[derive(PartialEq)]
pub enum VectorBase {
    Spherical,
    Cartesian,
}

pub struct VectorComponents{
    pub base:VectorBase,
    pub coords:na::Vector3<f64>,
}

impl VectorComponents{
    pub fn new_zero(base:VectorBase)->VectorComponents{
        VectorComponents { base: base, coords: na::Vector3::new(0.0,0.0,0.0) }
    }
}


//[Ricardo]: Configuration parameters of the simulation...I guess.
// In the c++ version this are called parameters.
#[derive(Debug,PartialEq)]
pub struct Config{
    pub n_modes:u16,
    pub l:Vec<u16>,
    pub m:Vec<i16>,
    pub rel_deltar:Vec<f64>,
    pub k:Vec<f64>,// K correction
    pub phase:Vec<f64>,
}

//[Ricardo]: Eigenfunctions of the linearized system...I guess. Since it's only two components
//(amplitud and phase difference(radians)) with respect to xir
//I presume It's from a fourier transformation. Later I might encapsulate the
//linear equation, eigenfunction, eigenvalues on the same data struct, since they all
//seem to be related
#[derive(Debug,PartialEq)]
pub struct Eigenfunctions{
    pub ampl: f64,
    pub phasedif: f64,
}
//[Ricardo]:This will be used in the function to Write the output of the pulstar code, I guess.
pub enum ToBeSaved {
    Velocity,
    Temperature,
    Gravity,
    CosineChi,
    Area,
}