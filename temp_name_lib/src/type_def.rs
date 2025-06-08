use nalgebra as na;

const THETA_STEP:u32 = 4; //step in grid in the theta direction
const PHI_STEP:u32 = 8; // step in grid in the phi direction
const N_FLUX_POINTS:u32=10000; // number of points in one flux profile
const MAX_N_MODES:u32 = 10; // max number of pulsation modes allowed
const MAX_N_TIMES:u32 = 3000; //max number of time points 

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
pub struct Config{
    pub n_modes:u16,
    pub l:Vec<u16>,
    pub m:Vec<i16>,
    pub rel_deltar:Vec<f64>,
    pub k:Vec<f64>,//[Ricardo]: It's the value along the k unit vector.I dont know anymore, it says kapa value.
    pub phase:Vec<f64>,
}

//[Ricardo]: Eigenfunctions of the linearized system...I guess. Since it's only two components
//(amplitud and phase difference(radians)) with respect to xir
//I presume It's from a fourier transformation. Later I might encapsulate the
//linear equation, eigenfunction, eigenvalues on the same data struct, since they all
//seem to be related
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