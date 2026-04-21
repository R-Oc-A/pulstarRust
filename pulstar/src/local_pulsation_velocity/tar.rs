use super::reference_frames::Coordinates;
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
            if let Some(hough_functions) = tar_functions{
                let index = reference_frames::rotation_treatment::tar::construct_index(theta, dtheta);
                let h_r = hough_functions.h_r[index];
                let h_p = hough_functions.h_p[index];
                let h_t = hough_functions.h_t[index];

                let radial_velocity = velocity_amplitude;
                let tangential_velocity = mode.k*radial_velocity;

                let v_r = -radial_velocity * h_r * (mode.phase + phi * mode.m as f64).sin();
                let v_theta = -tangential_velocity * h_t/sintheta *(mode.phase + phi * mode.m as f64).sin();
                let v_phi = tangential_velocity * h_p/sintheta * (mode.phase + phi * mode.m as f64).cos();

                Ok(Coordinates::Spherical(na::Vector3::new(v_r,v_theta, v_phi)))
            }else{
                Err(MathErrors::FunctionNotFound)
            }
        }
    }
}
