use crate::{joris_math::spherical_harmonics::plmcos::plmcos,
     type_def::VectorComponents,
     type_def::VectorBase};
use nalgebra as na;

/// Compute the Lagrangian displacement vector in spherical coordinates
/// 
/// # Arguments:
/// # `phase` - omega*t + psi         [rad]
/// * `theta` - colatitude coordinate [rad]
/// * `phi`   - azimuthal coordinate  [rad]
/// * `l`     - degree of the spherical harmonic >= 0
/// * `m`     -  azimuthal number of the spherical harmonic, -l <= m <= l
/// * `ampl_radial`     - amplitude in the radial direction * Y_l^m
/// * `ampl_tangential` - amplitude in the tangential direction * Y_l^m
///
impl VectorComponents{
pub fn displacement(phase: &f64,
                    theta: &f64, 
                    phi: &f64,
                    l: &i16, 
                    m: &i16, 
                    ampl_radial: &f64, 
                    ampl_tangential: &f64,
                    ) -> VectorComponents {

    let sintheta = theta.sin();
    let costheta = theta.cos();
    let plmcostheta = plmcos(l, m.abs() as u16, sintheta, costheta); 
    let dplmcostheta_dtheta = (- f64::from(l+1) * costheta * plmcostheta            // First derivative
                               + f64::from((l as i16) - m + 1) 
                               * plmcos(l+1, m.abs() as u16, sintheta, costheta))  
                               / sintheta;

    let delta_r     = ampl_radial * plmcostheta 
                           * f64::cos(phase + f64::from(m)*phi);
    let delta_theta = ampl_tangential * dplmcostheta_dtheta 
                           * f64::cos(phase + f64::from(m)*phi);
    let delta_phi   = ampl_tangential * f64::from(-m) * plmcostheta 
                           * f64::sin(phase + f64::from(m)*phi) 
                           / (sintheta * sintheta);

    let displacement_coordinates = VectorComponents{
        base: VectorBase::Spherical,
        coords:na::Vector3::new(delta_r, delta_theta, delta_phi),
    };
    return displacement_coordinates;
}

}