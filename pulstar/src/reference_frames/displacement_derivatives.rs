use crate::reference_frames::rotation_treatment::tar::TARCollection;
use crate::{PulsationMode, RotationRegime};
use super::{MathErrors};
use super::rotation_treatment::{self};
//? This module contains the functions to calculate the derivatives of the lagrangian displacement vector over 
//? the surface of a star using spherical coordinates. 




///Computes the derivatives of Δr/r0 with respect to θ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - sine of the colatitude angle (theta in rads)
/// * `dtheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to θ
pub fn d_dr_rdtheta(
    mode: &PulsationMode,
	theta: f64,
	phi: f64,
    tar_functions:&Option<TARCollection>) -> f64{

    match &mode.rotation_effects{
        RotationRegime::NonRotating =>{
            let sintheta = theta.sin(); 
            let costheta = theta.cos(); 
            rotation_treatment::
            non_rotating::non_rotating_d_dr_rdtheta(mode, sintheta, costheta, phi)},

        RotationRegime::PerturbativeCoriolis =>{ 
            let sintheta =theta.sin();
            let costheta = theta.cos();
            rotation_treatment::
            perturbative_coriolis::perturbative_d_dr_rdtheta(mode, sintheta, costheta, phi)},

        RotationRegime::Tar =>{ if let Some(houghs_functions) = tar_functions{
            rotation_treatment::
            tar::tar_d_dr_rdtheta(mode, theta,  phi,houghs_functions)}
            else{panic!("hough functions where not properly loaded.")}},

        RotationRegime::CentrifugalDeformation{coefficient_expansion:_} =>{ rotation_treatment::
            centrifugal_deformation::deformed_d_dr_rdtheta(mode, theta.sin(), theta.cos(), phi)},
    }
}

///Computes the derivatives of Δθ with respect to θ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - sine of the colatitude angle (theta in rads)
/// * `dtheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the displacement in θ with respect to θ
pub fn d_dtheta_dtheta(
    mode: &PulsationMode,
	theta: f64,
	phi: f64,
    spin_parameter:f64,
    tar_functions:&Option<TARCollection>) -> f64{

    match &mode.rotation_effects{
        RotationRegime::NonRotating =>{
            let sintheta = theta.sin(); 
            let costheta = theta.cos(); 
            rotation_treatment::
            non_rotating::non_rotating_d_dtheta_dtheta(mode, sintheta, costheta, phi)},

        RotationRegime::PerturbativeCoriolis =>{ 
            let sintheta = theta.sin();
            let costheta = theta.cos();
            rotation_treatment::
            perturbative_coriolis::perturbative_d_dtheta_dtheta(mode, sintheta, costheta, phi,spin_parameter).unwrap()},

        RotationRegime::Tar =>{ if let Some(houghs_functions) = tar_functions{
            rotation_treatment::
            tar::tar_d_dtheta_dtheta(mode, theta,  phi,houghs_functions).unwrap()}
            else{panic!("hough functions where not properly loaded.")}},

        RotationRegime::CentrifugalDeformation{coefficient_expansion:_} =>{ rotation_treatment::
            centrifugal_deformation::deformed_d_dtheta_dtheta(mode, theta.sin(), theta.cos(), phi)}
    }
}

///Computes the derivatives of Δr/r0 with respect to φ in the point with spherical
///coordinates θ,φ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - sine of the colatitude angle (theta in rads)
/// * `dtheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to φ 
pub fn d_dr_rdphi(
    mode: &PulsationMode,
	theta: f64,
	phi: f64,
    tar_functions:&Option<TARCollection>) -> f64{

    match &mode.rotation_effects{
        RotationRegime::NonRotating =>{
            let sintheta = theta.sin(); 
            let costheta = theta.cos(); 
        rotation_treatment::
            non_rotating::non_rotating_d_dr_rdphi(mode, sintheta,costheta, phi)},

        RotationRegime::PerturbativeCoriolis =>{ 
            let sintheta = theta.sin();
            let costheta = theta.cos();
            rotation_treatment::
            perturbative_coriolis::perturbative_d_dr_rdphi(mode, sintheta, costheta, phi)},

        RotationRegime::Tar =>{ if let Some(houghs_functions) = tar_functions{
            rotation_treatment::
            tar::tar_d_dr_rdphi(mode, theta,  phi,houghs_functions)}
            else{panic!("hough functions where not properly loaded.")}},

        RotationRegime::CentrifugalDeformation{coefficient_expansion:_} =>{ rotation_treatment::
            centrifugal_deformation::deformed_d_dr_rdphi(mode, theta.sin(), theta.cos(), phi)}
    }
}

///Computes the derivatives of Δϕ with respect to ϕ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - sine of the colatitude angle (theta in rads)
/// * `dtheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// This function returns a [Result] with the following variants:
/// * `Ok(f64)` - Where the binded value is the derivative of the displacement in φ with respect to φ 
/// * `Err(DivisionByZero)` - Where the binded error is returned to the calling function and indicates that the theta value was too small.
pub fn d_dphi_dphi(
    mode: &PulsationMode,
	theta: f64,
	phi: f64,
    spin_parameter:f64,
    tar_functions:&Option<TARCollection>) -> Result<f64,MathErrors>{

    match &mode.rotation_effects{
        RotationRegime::NonRotating =>{
            let sintheta = theta.sin(); 
            let costheta = theta.cos(); 
            rotation_treatment::
            non_rotating::non_rotating_d_dphi_dphi(mode,sintheta, costheta, phi)},

        RotationRegime::PerturbativeCoriolis =>{
            let sintheta =theta.sin();
            let costheta = theta.cos();
            rotation_treatment::
            perturbative_coriolis::perturbative_d_dphi_dphi(mode, sintheta, costheta, phi,spin_parameter)},

        RotationRegime::Tar =>{ if let Some(houghs_functions) = tar_functions{
            rotation_treatment::
            tar::tar_d_dphi_dphi(mode, theta,  phi,houghs_functions)}
            else{panic!("hough functions where not properly loaded.")}},
            
        RotationRegime::CentrifugalDeformation{coefficient_expansion:_} =>{ rotation_treatment::
            centrifugal_deformation::deformed_d_dphi_dphi(mode, theta.sin(), theta.cos(), phi)}
    }
}