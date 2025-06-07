use nalgebra as na;
use crate::type_def::{VectorComponents,VectorBase,Config};
use super::displacement;
use crate::joris_math::spherical_harmonics::norm_factor::ylmnorm;
///
///This function computes the spherical components of the       
///surface normal vector in a frame where the z-axis coincides  
///with the rotation axis. A surface normal is a vector which
///stand locally perpendicular to the surface and has a length
///equal to the area of the local surface cell.
///WARNING: this function does NOT multiply each component
///with R_0^2 d\theta d\phi. The user must do this him/herself.
///
/// REMARK: If the flag "restriction" is set to 0, then no time dependence
///         is taken into account (like with a static star). But the vector
///         still differs from cell to cell on the surface.
///         If the flag is different from 0, then the full time dependence
///         is taken into account.
///
///
    pub fn surface_normal(
        parameters: &Config,
        theta_rad: f64,
        phi_rad: f64,
        restriction: bool) -> VectorComponents {







        
        
        let mut total_p_ds = na::Vector3::new(0.0,0.0,0.0);
        for (n,_) in parameters.l.iter().enumerate(){
            
            let ampl_radial= parameters.rel_deltar[n] * ylmnorm(parameters.l[n], parameters.m[n]);
            let ampl_tangential= parameters.rel_deltar[n] * parameters.k[n]* ylmnorm(parameters.l[n], parameters.m[n]);
            
            let pulsational_displacement = 
            displacement::displacement(theta_rad,
                                        phi_rad,
                                        parameters, 
                                        n, 
                                        ampl_radial, 
                                        ampl_tangential);
            
            total_p_ds += pulsational_displacement.coords;

        }
        
        //
        VectorComponents{
            base: VectorBase::Spherical,
            coords: total_p_ds,
        }
     }

