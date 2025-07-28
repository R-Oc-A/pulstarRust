use nalgebra as na;
use super::*;

/// This enum binds the components of a 3D vector depending on the type of coordinates that are used
/// 
/// as of now there're the following:
/// * 'Spherical' - The component order of this variant is `(r,θ,φ)` 
/// * 'Cartesian' - The component order of this variant is `(x,y,z)`
pub enum Coordinates{
    Spherical(na::Vector3<f64>),
    Cartesian(na::Vector3<f64>),
}

/// Compute the Lagrangian displacement vector in spherical coordinates
/// 
/// ### Arguments:
/// * `pulse_params` - Parameters that characterize the pulsation modes and the star. See [crate::PPulstarConfig]
/// * 'index' - This unsigned integer indicates for which mode is the displacement calculated.
/// * `sintheta` - sine of the colatitude coordinate (theta in rad)
/// * 'costheta' - cosine of the colatitude coordinate (theta in rad)
/// * `phi`   - azimuthal coordinate  in rad
/// * `radial_amplitude`     - amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::joris_math::spherical_harmonics::norm_factors])
/// * `tangential_amplitude` - amplitude in the tangential direction times the normalization factor  'Y_l^m' (see [temp_name_lib::joris_math::spherical_harmonics::norm_factors])
/// 
/// ### Returns:
/// This function can return an [Ok] or [Err] variants of [Result] that will have the following values binded to them:
/// * `Ok(Coordinates::Spherical)` - an Ok  variant that has binded the spherical components of the displacement vector in the`r,θ,φ` order.
/// * `Err(DivisionByZero)` - an Err variant that has binded the error produced if the colatitude  coordinate (theta) is too small.
pub fn ddisplacement(pulse_params:&PPulstarConfig,
    index:usize,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64)->Result<Coordinates,MathErrors>{
        match sintheta.abs() <= f64::EPSILON.sqrt(){
            true => {Err(MathErrors::DivisionByZero)}
            false => {
                let phase = pulse_params.mode_data[index].phase_offset;
                let l = pulse_params.mode_data[index].l;
                let m =  pulse_params.mode_data[index].m;
                
                let plmcostheta = plmcos(l, m.abs() as u16, sintheta, costheta); 
                let dplmcostheta_dtheta = (- f64::from(l+1) * costheta * plmcostheta  // First derivative
                                        + f64::from((l as i16) - m + 1) 
                                        * plmcos(l+1, m.abs() as u16, sintheta, costheta))  
                                        / sintheta;

                let delta_r     = radial_amplitude * plmcostheta 
                                    * f64::cos(phase + f64::from(m)*phi);
                let delta_theta = tangential_amplitude * dplmcostheta_dtheta 
                                    * f64::cos(phase + f64::from(m)*phi);
                let delta_phi   = tangential_amplitude * f64::from(-m) * plmcostheta 
                                    * f64::sin(phase + f64::from(m)*phi) 
                                    / (sintheta.abs().powi(2));

                Ok(Coordinates::Spherical(na::Vector3::new(delta_r, delta_theta, delta_phi)))
            }
        }
    }


