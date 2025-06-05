use nalgebra as na;
use crate::type_def::VectorComponents;
use crate::type_def::VectorBase;
use crate::type_def::Config;
use crate::star_physics::displacement;
///
///
/// REMARK: If the flag "restriction" is set to 0, then no time dependence
///         is taken into account (like with a static star). But the vector
///         still differs from cell to cell on the surface.
///         If the flag is different from 0, then the full time dependence
///         is taken into account.
///
///
impl VectorComponents{
    pub fn surface_normal(
        parameters: Config,
        theta_rad: f64,
        phi_rad: f64,
        restriction: bool) -> VectorComponents {
        ///This function computes the spherical components of the 
        ///surface normal vector in a frame where the z-axis coincides
        ///with the rotation axis. A surface normal is a vector which
        ///stand locally perpendicular to the surface and has a length
        ///equal to the area of the local surface cell.
        ///WARNING: this function does NOT multiply each component
        ///with R_0^2 d\theta d\phi. The user must do this him/herself.
            
        let mut r = 0.0;
        let mut theta = 0.0;
        let mut phi = 0.0;
        
        let mut total_p_ds = na::Vector3::new(0.0,0.0,0.0);
                
        let mut n =1;
        while n<Config::n_modes{
            let pulsational_displacement = 
            VectorComponents::displacement(phase[n],
                                           theta,
                                           phi, 
                                           l[n], 
                                           m[n], 
                                           ampl_radial, 
                                           ampl_tangential);

            n+=1;
        }
        let mut 
        //
        let result = VectorComponents{
            base: VectorBase::Spherical,
            coords: [r,theta,phi]
        };
        return result;
     }
}

