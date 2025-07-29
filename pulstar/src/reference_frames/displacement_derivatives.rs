use crate::{PPulstarConfig, PulsationMode};
use super::{MathErrors,MACHINE_PRECISION};
use super::spherical_harmonics::{
    d_plmcos_dtheta::{
        deriv1_plmcos_dtheta as d_plmcos_dtheta,
        deriv2_plmcos_dtheta as d2_plmcos_dtheta}, norm_factor::ylmnorm, plmcos::plmcos
    }; 

//? This module contains the functions to calculate the derivatives of the lagrangian displacement vector over 
//? the surface of a star using spherical coordinates. 




///Computes the derivatives of Δr/r0 with respect to θ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to θ
pub fn d_dr_rdtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let r_dr = mode.rel_dr;
    let phase= mode.phase_offset;
    let l = mode.l;
    let m= mode.m;
                            
    r_dr*ylmnorm(l,m)
    * d_plmcos_dtheta(l,m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).cos()
}

///Computes the derivatives of Δθ with respect to θ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the displacement in θ with respect to θ
pub fn d_dtheta_dtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let r_dr = mode.rel_dr;
    let phase= mode.phase_offset;
    let k = mode.k;
    let l = mode.l;
    let m= mode.m;

    r_dr*ylmnorm(l,m)*k
    * d2_plmcos_dtheta(l,m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).cos()
}

///Computes the derivatives of Δr/r0 with respect to φ in the point with spherical
///coordinates θ,φ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to φ 
pub fn d_dr_rdphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let r_dr = mode.rel_dr;
    let phase= mode.phase_offset;
    let l = mode.l;
    let m= mode.m;
    
    r_dr * ylmnorm(l,m) * (-m as f64)
    * plmcos(l, m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).sin()
}

///Computes the derivatives of Δϕ with respect to ϕ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// This function returns a [Result] with the following variants:
/// * `Ok(f64)` - Where the binded value is the derivative of the displacement in φ with respect to φ 
/// * `Err(DivisionByZero)` - Where the binded error is returned to the calling function and indicates that the theta value was too small.
pub fn d_dphi_dphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> Result<f64,MathErrors>{

    match sintheta < MACHINE_PRECISION{  
        false => {
        let r_dr = mode.rel_dr;
        let phase= mode.phase_offset;
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