use nalgebra as na;
use crate::{type_def::{Config, VectorBase, VectorComponents}, utils::MathErrors};
use super::{displacement,displacement::derivatives};
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
        time_independent: bool) -> Result<VectorComponents,MathErrors> {
        
        let sintheta=theta_rad.sin();

        match time_independent{
        true => {
            Ok(VectorComponents{
                base: VectorBase::Spherical,
                coords: na::Vector3::new(sintheta,0.0,0.0),}
            )
        }
        false =>{
            let costheta = theta_rad.cos();
            
            let mut total_p_ds = na::Vector3::new(0.0,0.0,0.0);
            let mut total_dev1=0.0;
            let mut total_dev2=0.0;
            let mut total_dev3=0.0;
            let mut total_dev4=0.0;
            
            for (n,_) in parameters.l.iter().enumerate(){
                
                let ampl_radial= parameters.rel_deltar[n] 
                                    * ylmnorm(parameters.l[n], parameters.m[n]);
                let ampl_tangential= parameters.rel_deltar[n] 
                                        * parameters.k[n] * 
                                        ylmnorm(parameters.l[n], parameters.m[n]);
                
                let pulsational_displacement = 
                displacement::displacement(parameters,
                                            n,
                                            sintheta,
                                            costheta,
                                            phi_rad,
                                            ampl_radial, 
                                            ampl_tangential) ? ;
                let drdtheta =
                derivatives::d_dtheta::d_dr_rdtheta (parameters,
                                                     n, 
                                                     sintheta,
                                                     costheta,
                                                     phi_rad);
                let drdphi =
                derivatives::d_dphi::d_dr_rdphi(parameters, 
                                                n, 
                                                sintheta, 
                                                costheta, 
                                                phi_rad);
                let dtdtheta = 
                derivatives::d_dtheta::d_dtheta_dtheta(parameters,
                                                       n, 
                                                       sintheta, 
                                                       costheta, 
                                                       phi_rad);
                let dpdphi = 
                derivatives::d_dphi::d_dphi_dphi(parameters,
                                                 n, 
                                                 sintheta,
                                                 costheta,
                                                 phi_rad) ? ;

                total_p_ds += pulsational_displacement.coords;
                total_dev1 += drdtheta;
                total_dev2 += drdphi;
                total_dev3 += dtdtheta;
                total_dev4 += dpdphi;
            }
            let mut normal_vector_coords=
            na::Vector3::new(1.0,
                             total_p_ds[1]-total_dev1,
                             sintheta * total_p_ds[2] - total_p_ds[1]/sintheta);
            let lenght= sintheta * (1.0 + 2.0 * total_p_ds[0]
                             + costheta/sintheta * total_p_ds[1]
                             + total_dev3 + total_dev4);
            
            normal_vector_coords *= lenght;
            //
            Ok(VectorComponents{
                base: VectorBase::Spherical,
                coords: normal_vector_coords,
            })
        }
     }

}