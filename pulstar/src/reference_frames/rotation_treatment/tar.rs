use ndarray::prelude::Array1;
use temp_name_lib::utils::MathErrors;
use crate::{MeshConfig, PulsationMode, PulstarConfig, reference_frames::Coordinates};
use nalgebra as na;


/// Contains the quantities obtained by solving the eigenvalue problem of the Laplace Tidal Differential equation,
/// which are useful for computing the colatitudinal factor of the lagrangian displacement calculated via the traditional approximation of rotation 
pub struct TARCollection{
/// An [Array1<f64>] that contains all of the equally spaced μ=cos(θ) values in the range [-1,1]
    mu_values: Array1<f64>,
/// the calculated eigenvalue of the Laplace Tidal Differential Operator.
    lambda:f64,
/// An [Array1<f64>] that contains the radial  hough function 𝚯_r
    H_r:Array1<f64>,
/// An [Array1<f64>] that contains the colatitudinal  hough function 𝚯_θ
    H_t:Array1<f64>,
/// An [Array1<f64>] that contains the   azimuthal function 𝚯_ɸ
    H_p:Array1<f64>,
/// An [Array1<f64>] that contains  the derivative with respect to μ of the radial  hough function 𝚯_r
    dH_r:Array1<f64>,
/// An [Array1<f64>] that contains  the derivative with respect to μ of the colatitudinal hough function 𝚯_θ
    dH_t:Array1<f64>,
/// An [Array1<f64>] that contains  the derivative with respect to μ of the azimuthal hough function 𝚯_ɸ
    dH_p:Array1<f64>,
/// a [usize] that contains the number of points in each of the [TARCollection] members
    npts:usize,
}


impl PulsationMode{
    /// This method returns the spin parameter defined as 2Ω/ω where Ω is the rotation frequency and ω is the pulsation frequency. Both in cycles per day. 
    fn get_spin_parameter(&self,pulsconfig: &PulstarConfig)->f64{
        let rotation_frequency = pulsconfig.get_rotation_frequency();//in cycles per day
        2.0*rotation_frequency/self.frequency
    }

    /// This method returns an instance of the [TARCollection] for a given pulsation mode. 
    fn new_tar_collection(&self, pulsconfig: &PulstarConfig)->TARCollection{
        let q = self.get_spin_parameter(pulsconfig);
        
        let npts = (180.0/ match pulsconfig.mesh{
            MeshConfig::Sphere { theta_step, phi_step:_ }=>{theta_step}
        })as usize;
        
        let (lambda,
            mu_values,
            H_r,
            H_t,
            H_p,
            dH_r,
            dH_t,
            dH_p,
            )
            =temp_name_lib::math_module::hough::hough(q, self.l, self.m, npts, (self.l*(self.l+1)) as f64, true);
        
        //Because of Macro values, this doesn't work. I'm going to pass arround [vectors f64] from hough.
        TARCollection { mu_values: Array1::from_vec(mu_values),
            lambda: lambda,
            H_r: Array1::from_vec(H_r),
            H_t: Array1::from_vec(H_t),
            H_p: Array1::from_vec(H_p),
            dH_r: Array1::from_vec(dH_r),
            dH_t: Array1::from_vec(dH_t),
            dH_p: Array1::from_vec(dH_p),
            npts:npts}
    }


}

/// Compute the Lagrangian displacement vector in spherical coordinates for the TAR
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - The Colatitude coordinate (θ in rads)
/// * `dtheta` - 
/// * `phi`   - azimuthal coordinate (ɸ  in rads)
/// * `radial_amplitude`     - amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::math_module::spherical_harmonics::norm_factors])
/// * `tangential_amplitude` - amplitude in the tangential direction times the normalization factor  'Y_l^m' (see [temp_name_lib::math_module::spherical_harmonics::norm_factors])
/// 
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
                let index:usize = (theta/dtheta) as usize;

                let H_r = houghs_functions.H_r[index];
                let H_p= houghs_functions.H_p[index];
                let H_t = houghs_functions.H_t[index];

                // Im taking this expressions from Townsend 2020.
                let delta_r = radial_amplitude * H_r * (-mode.phase + phi * mode.m as f64).cos();
                let delta_theta = tangential_amplitude * H_t * (-mode.phase + phi * mode.m as f64).cos()/sintheta;
                let delta_phi = tangential_amplitude * H_p * (-mode.phase + phi * mode.m as f64).sin()/sintheta;

                Ok(Coordinates::Spherical(na::Vector3::new(delta_r, delta_theta, delta_phi)))
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
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to θ
pub fn tar_d_dr_rdtheta(
    mode: &PulsationMode,
	theta: f64,
    dtheta:f64,
	phi: f64,
    houghs_functions:&TARCollection
    ) -> f64{
    let index = (theta/dtheta) as usize;

    let dh_r = -houghs_functions.dH_r[index]*theta.sin();// dH_r is the derivative of H_r with respect to μ=cos(θ), so here I applied the chain rule.

    mode.rel_dr*dh_r
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
/// * an `f64` - This value is the derivative of the displacement in θ with respect to θ
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
            let index = (theta/dtheta) as usize;
            let dh_t= - houghs_functions.dH_t[index]*sintheta;
            let h_t = houghs_functions.H_t[index];
            Ok(mode.rel_dr*mode.k
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
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to φ 
pub fn tar_d_dr_rdphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let r_dr = mode.rel_dr;
    let phase= mode.phase;
    let l = mode.l;
    let m= mode.m;
    
    r_dr * ylmnorm(l,m) * (-m as f64)
    * plmcos(l, m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).sin()
}

///Computes the derivatives of Δϕ with respect to ϕ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// This function returns a [Result] with the following variants:
/// * `Ok(f64)` - Where the binded value is the derivative of the displacement in φ with respect to φ 
/// * `Err(DivisionByZero)` - Where the binded error is returned to the calling function and indicates that the theta value was too small.
pub fn tar_d_dphi_dphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> Result<f64,MathErrors>{

    match sintheta < MACHINE_PRECISION{  
        false => {
        let r_dr = mode.rel_dr;
        let phase= mode.phase;
        let k= mode.k;
        let l = mode.l;
        let m= mode.m;

        Ok(r_dr * k * ylmnorm(l, m) * (-(m as f64).powi(2))
        * plmcos(l, m.abs() as u16, sintheta, costheta)
        * (phase + (m as f64) * phi).cos()
        /(sintheta.abs().powi(2)) )
        }

        true =>{
        Err(MathErrors::DivisionByZero) //will pass the error in order for the calling function to do something
        }
    }
}