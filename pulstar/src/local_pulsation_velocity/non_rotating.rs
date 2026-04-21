use super::reference_frames::Coordinates;
use super::na;
use temp_name_lib::utils::{MACHINE_PRECISION,MathErrors};
use temp_name_lib::math_module::spherical_harmonics::plmcos::plmcos;
use temp_name_lib::math_module::spherical_harmonics::d_plmcos_dtheta::deriv1_plmcos_dtheta;
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
pub fn v_non_rotating(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi_rad:f64,
    velocity_amplitude:f64,
)->Result<Coordinates,MathErrors>{
    match sintheta.abs() <= MACHINE_PRECISION{
        true => { Err(MathErrors::DivisionByZero)}
        false => {
            let l=mode.l;
            let m = mode.m;
            let phase = mode.phase_offset;
            let k = mode.k;
            let v_r = velocity_amplitude * ylmnorm(l, m)
                * plmcos(l, m.abs() as u16, sintheta, costheta)
                * (phase + (m as f64) * phi_rad).sin();
            let v_theta = velocity_amplitude * k
                   * ylmnorm(l, m)
                   * deriv1_plmcos_dtheta(l, m.abs() as u16, sintheta, costheta)
                   * (phase + (m as f64) * phi_rad).sin();
            let v_phi = velocity_amplitude * k
                   * ylmnorm(l, m)
                   * (-(m as f64))
                   * plmcos(l, m.abs() as u16, sintheta, costheta)
                   * (phase * (m as f64) * phi_rad).cos()
                   / sintheta;
        Ok(Coordinates::Spherical(na::Vector3::new(v_r,v_theta,v_phi)))
        }
    }
}

