use std::{f64::consts::PI};

use ndarray::prelude::Array1;
use temp_name_lib::{math_module::spherical_harmonics::{norm_factor::ylmnorm, plmcos::plmcos}, utils::MathErrors};
use crate::{MeshConfig, PulsationMode, PulstarConfig, reference_frames::{Coordinates, ampl_r, ampl_t}};
use nalgebra as na;

pub const NPTS:usize = 1000usize; 


/// Contains the quantities obtained by solving the eigenvalue problem of the Laplace Tidal Differential equation,
/// which are useful for computing the colatitudinal factor of the lagrangian displacement calculated via the traditional approximation of rotation 
pub struct TARCollection{
/// An [Array1<f64>] that contains all of the equally spaced μ=cos(θ) values in the range [-1,1]
    mu_values: Array1<f64>,
/// the calculated eigenvalue of the Laplace Tidal Differential Operator.
    lambda:f64,
/// An [Array1<f64>] that contains the radial  hough function 𝚯_r
    pub h_r:Array1<f64>,
/// An [Array1<f64>] that contains the colatitudinal  hough function 𝚯_θ
    pub h_t:Array1<f64>,
/// An [Array1<f64>] that contains the   azimuthal function 𝚯_ɸ
    pub h_p:Array1<f64>,
/// An [Array1<f64>] that contains  the derivative with respect to μ of the radial  hough function 𝚯_r
    dh_r:Array1<f64>,
/// An [Array1<f64>] that contains  the derivative with respect to μ of the colatitudinal hough function 𝚯_θ
    dh_t:Array1<f64>,
/// An [Array1<f64>] that contains  the derivative with respect to μ of the azimuthal hough function 𝚯_ɸ
    dh_p:Array1<f64>,
/// a [usize] that contains the number of points in each of the [TARCollection] members
    npts:usize,
}


impl PulsationMode{

    /// This method returns an instance of the [TARCollection] for a given pulsation mode. 
    pub fn new_tar_collection(&self, pulsconfig: &PulstarConfig)->TARCollection{
        let q = self.get_spin_parameter(pulsconfig);
        println!("spin parameter is {}",q);
        let mut npts = (180.0/ match pulsconfig.mesh{
            MeshConfig::Sphere { theta_step, phi_step:_ }=>{theta_step}
        })as usize;
        if npts < NPTS{
            npts = NPTS;
        }
        let (lambda,
            mu_values,
            h_r,
            h_t,
            h_p,
            dh_r,
            dh_t,
            dh_p,
            )
            =temp_name_lib::math_module::hough::hough(q, self.l, self.m, npts, (self.l*(self.l+1)) as f64, true);
            //=temp_name_lib::math_module::hough::hough(q, self.l, self.m, npts, -(self.m.pow(2)) as f64, true);
        println!("lambda is {}, and l*(l+1) is {}",lambda, self.l*(self.l + 1));


        //----------------------------------------"
        //    Renormalizing Hough functions"
        //----------------------------------------"
        //Calculate the maximum of the associated Legendre polynomial evaluated in the mu=cos(θ) array;
        let max_plm = mu_values.iter().fold(
            plmcos(self.l, self.m.abs() as u16, (1.0-mu_values[0].powi(2)).sqrt(), mu_values[0]).abs(),
            |acc,x|{
                let plmcostheta = plmcos(self.l, self.m.abs() as u16, (1.0-x.powi(2)).sqrt(), *x).abs();
                if plmcostheta > acc{plmcostheta} else{acc}
            }
        );



        TARCollection { mu_values: Array1::from_vec(mu_values),
            lambda: lambda,
            h_r: Array1::from_vec(h_r)*1.0/max_plm,
            h_t: Array1::from_vec(h_t)*1.0/max_plm,
            h_p: Array1::from_vec(h_p)*1.0/max_plm,
            dh_r: Array1::from_vec(dh_r)*1.0/max_plm,
            dh_t: Array1::from_vec(dh_t)*1.0/max_plm,
            dh_p: Array1::from_vec(dh_p)*1.0/max_plm,
            npts:npts}
    }


}

/// Compute the Lagrangian displacement vector in spherical coordinates for the TAR
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - The Colatitude coordinate (θ in rads)
/// * `dtheta` - 
/// * `phi`   - azimuthal coordinate (φ  in rads)
/// * `radial_amplitude`     - amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::math_module::spherical_harmonics::norm_factors])
/// * `tangential_amplitude` - amplitude in the tangential direction times the normalization factor  'Y_l^m' (see [temp_name_lib::math_module::spherical_harmonics::norm_factors])
/// * `hough_functions` - a call by reference to an instance of [TARCollection] that contains the hough functions and derivatives for all theta angles. 
/// ### Returns:
/// This function can return an [Ok] or [Err] variants of [Result] that will have the following values binded to them:
/// * `Ok(Coordinates::Spherical)` - an Ok  variant that has binded the spherical components of the displacement vector in the`r,θ,φ` order.
/// * `Err(DivisionByZero)` - an Err variant that has binded the error produced if the colatitude  coordinate (theta) is too small.
pub fn tar_displacement(
    mode: &PulsationMode,
    theta:f64,
    dtheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64,
    houghs_functions:&TARCollection)->Result<Coordinates,MathErrors>{
        let sintheta = theta.sin();
        match sintheta.abs()<=f64::EPSILON.sqrt(){
            true => {Err(MathErrors::DivisionByZero)}
            false =>{
                let index = construct_index(theta, dtheta);
                
                let h_r = houghs_functions.h_r[index];
                let h_p= houghs_functions.h_p[index];
                let h_t = houghs_functions.h_t[index];
                
                // Im taking this expressions from Townsend 2020.
                let delta_r = radial_amplitude * h_r * (mode.phase + phi * mode.m as f64).cos();
                let delta_theta = tangential_amplitude * h_t/sintheta * (mode.phase + phi * mode.m as f64).cos();
                let delta_phi = -(mode.m as f64) * tangential_amplitude * h_p / sintheta.powi(2) * (mode.phase + phi * mode.m as f64).sin();

                Ok(Coordinates::Spherical( (na::Vector3::new(delta_r, delta_theta, delta_phi)) ))
            }
        }
    }

///Computes the derivatives of Δr/r0 with respect to θ in the point with spherical
///coordinates θ,ϕ using the TAR approach
/// 
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - The colatitude angle (theta in rads)
/// * `dtheta` - The resolution of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// * `hough_functions` - a call by reference to an instance of [TARCollection] that contains the hough functions and derivatives for all theta angles. 
/// ### Returns:
/// * a [f64] value of the derivative with respect of θ of the radial component of the lagrangian displacement. 
pub fn tar_d_dr_rdtheta(
    mode: &PulsationMode,
	theta: f64,
    dtheta:f64,
	phi: f64,
    houghs_functions:&TARCollection
    ) -> f64{
    let index = construct_index(theta, dtheta);

    let dh_r = -houghs_functions.dh_r[index]*theta.sin();// dH_r is the derivative of H_r with respect to μ=cos(θ), so here I applied the chain rule.

    ampl_r(mode)*dh_r
    * (mode.phase + (mode.m as f64) * phi).cos()
}

///Computes the derivatives of Δθ with respect to θ in the point with spherical
///coordinates θ,ϕ using the TAR approach
/// 
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - The colatitude angle (theta in rads)
/// * `dtheta` - The resolution of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// * `hough_functions` - a call by reference to an instance of [TARCollection] that contains the hough functions and derivatives for all theta angles. 
/// ### Returns:
/// This function can return an [Ok] or [Err] variants of [Result] that will have the following values binded to them:
/// * `Ok(f64)` - an [Ok]  variant that has binded the derivative of the theta displacement with respect to θ.
/// * `Err(DivisionByZero)` - an Err variant that has binded the error produced if the colatitude  coordinate (theta) is too small.
pub fn tar_d_dtheta_dtheta(
    mode: &PulsationMode,
	theta: f64,
    dtheta:f64,
	phi: f64,
    houghs_functions:&TARCollection) -> Result<f64,MathErrors>{
    
    let sintheta=theta.sin();
    match sintheta.abs() <= f64::EPSILON.sqrt(){
        true => { Err(MathErrors::DivisionByZero)}
        false => {
            let index = construct_index(theta, dtheta);
            let dh_t= - houghs_functions.dh_t[index]*sintheta;
            let h_t = houghs_functions.h_t[index];
            Ok(ampl_t(mode)
            * (dh_t/sintheta - h_t/sintheta.powi(2)*theta.cos())
            * (mode.phase + (mode.m as f64) * phi).cos()
        )
        }
    }
}

///Computes the derivatives of Δr/r0 with respect to φ in the point with spherical
///coordinates θ,φ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - sine of the colatitude angle (theta in rads)
/// * `dtheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// * `hough_functions` - a call by reference to an instance of [TARCollection] that contains the hough functions and derivatives for all theta angles. 
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to φ 
pub fn tar_d_dr_rdphi(
    mode: &PulsationMode,
    theta: f64,
    dtheta: f64,
	phi: f64,
    houghs_functions:&TARCollection) -> f64{
        let index = construct_index(theta, dtheta);
        let h_r = houghs_functions.h_r[index];
        - ampl_r(mode) * h_r * mode.m as f64
            *(mode.phase + (mode.m as f64)*phi).sin()

}
   
///Computes the derivatives of Δφ with respect to φ in the point with spherical
///coordinates θ,φ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - sine of the colatitude angle (theta in rads)
/// * `dtheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// * `hough_functions` - a call by reference to an instance of [TARCollection] that contains the hough functions and derivatives for all theta angles. 
/// ### Returns:
/// This function returns a [Result] with the following variants:
/// * `Ok(f64)` - Where the binded value is the derivative of the displacement in φ with respect to φ 
/// * `Err(DivisionByZero)` - Where the binded error is returned to the calling function and indicates that the theta value was too small.
pub fn tar_d_dphi_dphi(
    mode: &PulsationMode,
    theta: f64,
    dtheta: f64,
	phi: f64,
    houghs_functions:&TARCollection) -> Result<f64,MathErrors>{
    let sintheta = theta.sin();
    match sintheta < f64::EPSILON.sqrt(){  
        false => {
            let index = construct_index(theta, dtheta);
            let h_p = houghs_functions.h_p[index];
        
            Ok(ampl_t(mode)* (-(mode.m.pow(2) as f64))
            * h_p
            * (mode.phase + (mode.m as f64) * phi).cos()
            /(sintheta.abs().powi(2)) )
        }

        true =>{Err(MathErrors::DivisionByZero)}
    }
}

/// Houghs functions are computed simultaneously on an array of theta values, thus it's necessary to provide an index to get the expected value for a given theta. 
/// This is the function that does that. 
/// ### Arguments:
/// * `theta` - the colatitude coordinate
/// * `dtheta` - The difference between anytwo consecutive theta values of the theta array. It must have the same units as theta (radians, degrees)
/// ### Returns:
/// * `index` - a [usize] value that indicates the position of a given theta in the theta array
pub fn construct_index(theta:f64,dtheta:f64)->usize{
    let mut npts = (PI/dtheta).floor() as usize;
        if npts < NPTS{
            npts = NPTS;
        }
    let index =(theta/dtheta).floor() as usize * npts/((PI/dtheta).floor() as usize);
    index
}