use super::reference_frames::Coordinates;
use crate::reference_frames::displacement;
use crate::reference_frames::rotation_treatment::tar::TARCollection;
use temp_name_lib::utils::MathErrors;

use super::*;

/// This function computes the components v_r,v_θ,v_φ of the pulsation velocity on a given surface cell of the star.
/// using the TAR approach of pulsation velocities
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` -  colatitude coordinate u θ in radians
/// * 'dtheta' - angular displacement in the colatitude coordinate from one neighbouring [SurfaceCell] to the next, it has the same units as θ.
/// * `phi`   - azimuthal coordinate  in rads
/// * `velocity_amplitude`     - Amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::math_module::spherical_harmonics::norm_factors]) in km/s
/// * `tar_functions` - An [Option] enum that has the following variants:
///     - [Some] variant that has binded a reference to a [TARCollection]
///     - [None] in case tar functions are not needed.
/// ### Returns:
/// This function returns [Ok] or [Err] variants of [Result]
/// * Ok(`Coordinates::Spherical(pulsation_velocity)`) - Where pulsation_velocity  is a [na::Vector3]
/// * Err(DivisionByZero) - Where the error is pased to the calling function if the colatitude angle θ is too small. 
pub fn v_tar(
    mode: &PulsationMode,
    theta: f64,
    dtheta: f64,
    phi:f64,
    velocity_amplitude:f64,
    tar_functions:&Option<TARCollection>,
)->Result<Coordinates,MathErrors>{
    let sintheta = theta.sin();
    match sintheta.abs() <= f64::EPSILON.sqrt(){
        true => { Err(MathErrors::DivisionByZero)}
        false => {
            let v_tangential = mode.k * velocity_amplitude;
            displacement(mode, theta, dtheta, phi, velocity_amplitude, v_tangential, tar_functions)
        }
    }
}
