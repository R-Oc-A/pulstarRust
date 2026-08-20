use super::reference_frames::rotation_treatment::centrifugal_deformation::add_deformations_generic;

use super::*;

/// This function computes the components v_r,v_θ,v_φ of the pulsation velocity on a given surface cell of the star.
/// using the non-rotating approach of pulsation velocities
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude coordinate (theta in rads)
/// * 'costheta' - cosine of the colatitude coordinate (theta in rads)
/// * `phi`   - azimuthal coordinate  in rads
/// * `velocity_amplitude`     - Amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::math_module::spherical_harmonics::norm_factors]) in km/s
/// 
/// ### Returns:
/// This function returns [Ok] or [Err] variants of [Result]
/// * Ok(`Coordinates::Spherical(pulsation_velocity)`) - Where pulsation_velocity  is a [na::Vector3]
/// * Err(DivisionByZero) - Where the error is pased to the calling function if the colatitude angle θ is too small. 
pub fn v_deformed(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    velocity_amplitude:f64,
)->Result<Coordinates,MathErrors>{

    let displacement_func = 
    |mode:&PulsationMode,
    sintheta:f64,
    costheta: f64,
    phi:f64,
    velocity_amplitude: f64,
    _dummy: f64,|{v_non_rotating(mode, sintheta, costheta, phi, velocity_amplitude)};

    let deformed_velocity = add_deformations_generic(displacement_func,
    mode, sintheta, costheta, phi, velocity_amplitude, 0.0);

    deformed_velocity
}

