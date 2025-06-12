use crate::type_def::{VectorComponents, VectorBase};
use crate::utils::{self,MathErrors};
use nalgebra as na;


/// Given the inclination angle i, this function returns a unit vector k which points in the direction
/// of the observer. The cartesian components are in a reference frame where the z-axis coincides with
/// the rotation axis. 
/// 
/// The inclination angle is in radians. 
pub fn unit_vector_k(inclination_angle: f64) -> VectorComponents {
    VectorComponents { base: VectorBase::Cartesian, 
                       coords: na::Vector3::new(-inclination_angle.sin(),
                                               0.0,
                                               inclination_angle.cos()) }
}

///Transform a vector from spherical to cartesian coordinates 
///Input: vector (vr,vtheta,vphi), theta and phi coordinates
/// 
pub fn spherical_to_cartesian(vector: &VectorComponents,
                              theta: f64,
                              phi: f64)->Result<VectorComponents,MathErrors>{

    match vector.base == VectorBase::Spherical {
        true => {
            Ok(VectorComponents{
                base: VectorBase::Cartesian,
                coords: transformation_matrix(&vector.base,
                                              theta,
                                              phi) * vector.coords,
            })
        }
        false =>{
            Err(MathErrors::DifferentVectorBase)
        }
    }
}

pub fn cartesian_to_spherical(vector: &VectorComponents,
                              theta: f64,
                              phi: f64)->Result<VectorComponents,MathErrors>{

    match vector.base == VectorBase::Cartesian {
        true => {
            Ok(VectorComponents{
                base: VectorBase::Spherical,
                coords: transformation_matrix(&vector.base,
                                              theta,
                                              phi) * vector.coords,
            })
        }
        false =>{
            Err(MathErrors::DifferentVectorBase)
        }
    }
}


//Returns a transformation matrix from spherical->cartesian or cartesian->spherical
fn transformation_matrix(original_vector_base: &VectorBase,
                         theta:f64,
                         phi: f64
                         ) -> na::Matrix3<f64>{

    let sintheta = theta.sin();
    let costheta = theta.cos();
    let sinphi = phi.sin();
    let cosphi = phi.cos();

    let t_matrix = na::Matrix3::new(
                sintheta*costheta, costheta*cosphi, -sinphi,
                sintheta*sinphi, costheta*cosphi, cosphi,
                costheta, -sintheta, 0.0,
            );//spherical to cartesian matrix

    match original_vector_base {
        VectorBase::Spherical => { t_matrix}
        VectorBase::Cartesian => { t_matrix.transpose()}
    }

}