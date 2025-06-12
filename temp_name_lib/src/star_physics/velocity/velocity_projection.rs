use crate::{joris_math::{geometry::projections::project_vector,
                        ref_frame_convrs::spherical_to_cartesian}, 
            type_def::{Config, VectorBase, VectorComponents}, utils::MathErrors, 
            };
use nalgebra as na;
use super::pulsation;

///This function computes the projected(on the line of sight) pulsational velocity.
/// The result will have the same units as the velocity_amplitude.
pub fn project_vpuls(parameters:&Config,
                     theta: f64,
                     phi: f64,
                     k: &VectorComponents,
                     velocity_amplitude: &[f64],
                    ) -> Result<f64,MathErrors>{
    
    let sintheta = theta.sin();
    let costheta = theta.cos();
    let mut velocity_pulsation = VectorComponents{
                                    base: VectorBase::Spherical,
                                    coords:na::Vector3::new(0.0,0.0,0.0),
                                    };
    for (n,item) in parameters.l.iter().enumerate(){                
    velocity_pulsation.coords += pulsation::v_puls(parameters,
                                                    n,
                                                    sintheta,
                                                    costheta,
                                                    phi,
                                                    velocity_amplitude)?.coords;
    }
    let cartesian_vel = spherical_to_cartesian(
        &velocity_pulsation, theta, phi) ?;
    Ok(project_vector(&cartesian_vel, k)?)
}

///Computes the projected (on the line of sight) rotational velocity. 
///The result will have the same dimensions as V_e. We assume
///a uniform rotation with the rotation axis the positive z-axis. 
pub fn project_vrot(V_e:f64,
                    theta:f64,
                    phi: f64,
                    k: &VectorComponents)->Result<f64,MathErrors>{
    let v_omega = VectorComponents{
                base:VectorBase::Cartesian,
                coords:V_e * na::Vector3::new
                (theta.sin()*phi.sin(),
                theta.sin()*phi.cos(),
                0.0),};

    Ok(project_vector(&v_omega, k)?)
}