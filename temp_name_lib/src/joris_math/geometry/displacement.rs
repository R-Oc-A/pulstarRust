use std::iter::Map;

use crate::{joris_math::spherical_harmonics::plmcos::plmcos,
     type_def::{Config, VectorBase, VectorComponents},
     utils::MathErrors};
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
pub fn displacement(parameters: &Config,
                    index: usize,
                    sintheta: f64, 
                    costheta: f64,
                    phi: f64,
                    ampl_radial: f64, 
                    ampl_tangential: f64,
                    ) -> Result<VectorComponents,MathErrors> {

    let machine_precision_value = sintheta.abs()< 1.0e-8;

    match machine_precision_value {
        
        false => { 
            let phase=parameters.phase[index];
            let l = parameters.l[index];
            let m=parameters.m[index];

            let plmcostheta = plmcos(l, m.abs() as u16, sintheta, costheta); 
            let dplmcostheta_dtheta = (- f64::from(l+1) * costheta * plmcostheta  // First derivative
                                    + f64::from((l as i16) - m + 1) 
                                    * plmcos(l+1, m.abs() as u16, sintheta, costheta))  
                                    / sintheta;

            let delta_r     = ampl_radial * plmcostheta 
                                * f64::cos(phase + f64::from(m)*phi);
            let delta_theta = ampl_tangential * dplmcostheta_dtheta 
                                * f64::cos(phase + f64::from(m)*phi);
            let delta_phi   = ampl_tangential * f64::from(-m) * plmcostheta 
                                * f64::sin(phase + f64::from(m)*phi) 
                                / (sintheta.abs().powi(2));

            Ok(VectorComponents{
                base: VectorBase::Spherical,
                coords:na::Vector3::new(delta_r, delta_theta, delta_phi),
            })
        }
        true => {
            Err(MathErrors::DivisionByZero)
        }
    }
}

pub mod derivatives;