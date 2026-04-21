use crate::reference_frames::rotation_treatment::perturbative_coriolis::{amplitude_lm1, amplitude_lp1};

use super::reference_frames::Coordinates;
use super::na;
use temp_name_lib::utils::{MathErrors};
use temp_name_lib::math_module::spherical_harmonics::plmcos::plmcos;
use temp_name_lib::math_module::spherical_harmonics::d_plmcos_dtheta::{deriv1_plmcos_dtheta as d_plmcos_dtheta};
use temp_name_lib::math_module::spherical_harmonics::norm_factor::ylmnorm;

use super::*;

/// This function computes the components v_r,v_θ,v_φ of the pulsation velocity on a given surface cell of the star.
/// using the non-rotating approach of pulsation velocities
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude coordinate (theta in rads)
/// * 'costheta' - cosine of the colatitude coordinate (theta in rads)
/// * `phi_rad`   - azimuthal coordinate  in rads
/// * `velocity_amplitude`     - Amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::math_module::spherical_harmonics::norm_factors]) in km/s
/// 
/// ### Returns:
/// This function returns [Ok] or [Err] variants of [Result]
/// * Ok(`Coordinates::Spherical(pulsation_velocity)`) - Where pulsation_velocity  is a [na::Vector3]
/// * Err(DivisionByZero) - Where the error is pased to the calling function if the colatitude angle θ is too small. 
pub fn v_perturbative(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    velocity_amplitude:f64,
    spin_parameter:f64
)->Result<Coordinates,MathErrors>{

    // Spheroidal component of the velocity. 
    let spheroidal_part = local_pulsation_velocity::non_rotating::
    v_non_rotating(mode, sintheta, costheta, phi, velocity_amplitude)?;

    let phase = mode.phase;
    let l = mode.l;
    let lp1 = l+1;
    let lm1 = l-1;
    let m = mode.m.abs() as u16;
    let k = mode.k;


    let plp1m = plmcos(lp1,m, sintheta, costheta); 
    let dplp1m = d_plmcos_dtheta(lp1, m, sintheta, costheta);
    let plm1m = plmcos(lm1, m, sintheta, costheta); 
    let dplm1m = d_plmcos_dtheta(lm1, m, sintheta, costheta);

    // First toroidal component.
    let radial_velocity = velocity_amplitude * ylmnorm(mode.l, mode.m);
    let tangential_velocity = amplitude_lp1(radial_velocity,
         spin_parameter, mode.l as f64, mode.m as f64, k);
    
    let v_r = 0.0;
    let v_theta = -tangential_velocity/sintheta * plp1m * (mode.m as f64)
        * (phase + 0.5*PI + (mode.m as f64)*phi).cos();
    let v_phi   = tangential_velocity * dplp1m 
        * (phase + 0.5*PI + (mode.m as f64)*phi).sin();
    
    let first_toroidal_part = Coordinates::Spherical(
        na::Vector3::new(v_r,v_theta,v_phi));
        
    // Second toroidal component.
    let radial_velocity = velocity_amplitude * ylmnorm(mode.l, mode.m);
    let tangential_velocity = amplitude_lm1(radial_velocity,
         spin_parameter, mode.l as f64, mode.m as f64, k);
    
    let v_r = 0.0;
    let v_theta = -tangential_velocity/sintheta * plm1m * (mode.m as f64)
        * (phase + 0.5*PI + (mode.m as f64)*phi).cos();
    let v_phi   = tangential_velocity * dplm1m 
        * (phase + 0.5*PI + (mode.m as f64)*phi).sin();
    
    let second_toroidal_part = Coordinates::Spherical(
        na::Vector3::new(v_r,v_theta,v_phi));

    
    Ok( ((spheroidal_part + first_toroidal_part)?
        + second_toroidal_part)? )

}
