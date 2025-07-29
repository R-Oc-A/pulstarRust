use super::PPulstarConfig;
use super::reference_frames::Coordinates;
use super::na;
use temp_name_lib::type_def::{CYCLI2RAD, RADIUSSUN};
use temp_name_lib::utils::{MACHINE_PRECISION,MathErrors};
use temp_name_lib::joris_math::spherical_harmonics::plmcos::plmcos;
use temp_name_lib::joris_math::spherical_harmonics::d_plmcos_dtheta::deriv1_plmcos_dtheta;
use temp_name_lib::joris_math::spherical_harmonics::norm_factor::ylmnorm;

use super::*;

impl PPulstarConfig {
    /// This method calculates the amplitudes of the pulsation velocities per mode in km/s
    /// 
    /// ### Arguments: 
    /// * ` ` - 
    /// 
    /// ### Returns:
    /// * `velocity_amplitudes` - A [Vec] collection of the expected velocity amplitudes (with `f64` values) per mode. This collection is ordered in a way that there's a match with the pulsation mode.
    pub fn get_velocity_amplitudes(&self)->Vec<f64>{
        let mut velocity_amplitudes_allmodes:Vec<f64> = Vec::new();

        for mode in self.mode_data.iter(){
            let freqrad = mode.frequency*CYCLI2RAD;
            let velocity_amplitude = self.star_data.radius
                * freqrad * RADIUSSUN * 1.03e-3
                * mode.rel_dr;
            velocity_amplitudes_allmodes.push(velocity_amplitude);
        }

        velocity_amplitudes_allmodes
    }
}


/// This function computes the components v_r,v_θ,v_φ of the pulsation velocity on a given surface cell of the star.
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// * `sintheta` - sine of the colatitude coordinate (theta in rads)
/// * 'costheta' - cosine of the colatitude coordinate (theta in rads)
/// * `phi_rad`   - azimuthal coordinate  in rads
/// * `velocity_amplitude`     - Amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::joris_math::spherical_harmonics::norm_factors]) in km/s
/// 
/// ### Returns:
/// This function returns [Ok] or [Err] variants of [Result]
/// * Ok(`Coordinates::Spherical(pulsation_velocity)`) - Where pulsation_velocity  is a [na::Vector3]
/// * Err(DivisionByZero) - Where the error is passed to the calling function if the colatitude angle θ is too small. 
pub fn v_pulse_single_mode(
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

/// This function computes the projected(on the line of sight) pulsational velocity.
/// the result will have the same units as the velocity amplitudes.
/// 
/// ### Arguments:
/// * `parameters` - The data contained in [PPulstarConfig], here you find the parameters that describe the pulsation modes and the star.
/// * `theta_rad` - The colatitude angle in rads, must not be too small in order to avoid the poles.
/// * `phi_rad` - The azimuthal angle in rads
/// * `k` - Unit vector in directed towards the observer it is prefered to be in spherical coordinates.
/// * `velocity_amplitudes` - A [Vec] collection of the expected velocity amplitudes (with `f64` values) per mode. This collection is ordered in a way that there's a match with the pulsation mode in km/s.
/// 
/// ### Returns:
/// This function returns a [Result] with the following variants:
/// * `Ok(projected_amplitude)` - Where `projected_amplitude` is a `f64` value that contains the sum of all the pulsation velocities on a specific surface cell projected to the line of sight;
/// * 'Err(DivisionByZero)` - Where the error is passed to the calling function in case that the colatitude angle θ is too small;
pub fn observed_pulsation_velocity(
    parameters:&PPulstarConfig,
    theta_rad:f64,
    phi_rad:f64,
    k:&Coordinates,
    )->Result<f64,MathErrors>{
    
    let sintheta = theta_rad.sin();
    let costheta = theta_rad.cos();

    let mut collection_velocities:Vec<Coordinates>=Vec::new();

    let velocity_amplitudes = parameters.get_velocity_amplitudes();
    
    for (index,mode) in parameters.mode_data.iter().enumerate(){
        collection_velocities.push(v_pulse_single_mode(
            mode,
            sintheta,
            costheta,
            phi_rad,
            velocity_amplitudes[index])?);
    }
    let sum_velocities = collection_velocities.iter()
        .fold(collection_velocities[0],
         |accumulator,item| (accumulator + *item).unwrap());

    match k {
        Coordinates::Spherical(_)=>{
            Ok(sum_velocities.project_vector(&k)?)
        }
        Coordinates::Cartesian(_)=>{
            let k_spherical = k.transform(theta_rad, phi_rad);
            Ok(sum_velocities.project_vector(&k_spherical)?)
        }
    }
}