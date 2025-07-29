
use super::{*,na};
use spherical_harmonics::{self,plmcos::plmcos};
use temp_name_lib::joris_math::{spherical_harmonics::norm_factor::ylmnorm};
use std::ops;




/// This enum binds the components of a 3D vector depending on the type of coordinates that are used
/// 
/// as of now there're the following:
/// * 'Spherical' - The component order of this variant is `(r,θ,φ)` 
/// * 'Cartesian' - The component order of this variant is `(x,y,z)`
#[derive(Debug,PartialEq,Clone,Copy)]
pub enum Coordinates{
    Spherical(na::Vector3<f64>),//<-[Ricardo:]  I'll be using the nalgebra crate as it's well suited for linear algebra operation on small fixed size arrays
    Cartesian(na::Vector3<f64>),
}
///This is operation overloading so that we can add vectors component by component on the same basis.
impl ops::Add for Coordinates{
    type Output = Result<Coordinates,MathErrors>;
    fn add(self,other: Coordinates) -> Result<Coordinates,MathErrors>{
        match self {
            Coordinates::Spherical(value) => {
                match other{
                    Coordinates::Spherical(value2)=>{
                        Ok(Coordinates::Spherical(value + value2))
                    }
                    Coordinates::Cartesian(value2)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
            Coordinates::Cartesian(value) =>{
                match other{
                    Coordinates::Cartesian(value2)=>{
                        Ok(Coordinates::Spherical(value + value2))
                    }
                    Coordinates::Spherical(value2)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
        }
    }
}

///This is operation overloading so that we can addasign (self += rhs) vectors component by component having the same basis.
impl ops::AddAssign for Coordinates{
    fn add_assign(&mut self, rhs: Self) {
        let dummy = (*self + rhs).unwrap();
        *self = dummy;       
    }

}

impl Coordinates {
    /// Projects a vector onto another one
    /// ### Arguments:
    /// * `other` - A borrowed reference variant of [Coordinates]
    /// ### Returns:
    /// This method returns A [Result] that has the following values binded to it:
    /// * `Ok(f64)` - where the binded value is the projection of the vector into `other`
    /// * `Err(DiferentVectorBase)` - This returns a Different Vector Base to the error to the caller function.
    pub fn project_vector(&self,other:&Self)->Result<f64,MathErrors>{
        match self{
            Coordinates::Cartesian(value)=>{
                match other {
                    Coordinates::Cartesian(value2)=>{
                        Ok(value.dot(value2))
                    }
                    Coordinates::Spherical(value2)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
            Coordinates::Spherical(value)=>{
                match other {
                    Coordinates::Spherical(value2)=>{
                        Ok(value.dot(value2))
                    }
                    Coordinates::Cartesian(value2)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
            }
        }
    /// This method produces a unit vector k in cartesian coordinates that points towards the observer, 
    /// in this frame of reference, the z-axis coincides with the rotation axis. 
    /// 
    /// ### Arguments: 
    /// * `inclination_angle` - The inclination angle in rads with respect to the observer. 
    /// ### Returns:
    /// * `Coordinate::Cartesian(unit_vector_k)` - The unit vector k in cartesian coordinates.
    pub fn unit_vector_k(inclination_angle:f64)->Self{
        Coordinates::Cartesian(na::Vector3::new(-inclination_angle.sin(),
            0.0,
            inclination_angle.cos() ))  
    }

    /// This method provides a straightforward way to compute the norm of a vector
    /// 
    /// ### Arguments:
    /// 
    /// ### Returns: 
    /// * `f64` - the L2 norm of a vector.
    pub fn vector_length(&self)->f64{
        match self {
            Coordinates::Cartesian(value)=>{value.norm()}
            Coordinates::Spherical(value)=>{value.norm()}
        }
    }

    /// This method obtains the transformation matrix useful for changing basis
    /// 
    /// ### Arguments
    /// * `theta_rad` - The colatitude angle in radians
    /// * `phi_rad` - The azimuthal angle in radians
    /// ### Returns
    /// * `t_a` -The transformation Matrix to change from spherical to cartesian coordinates and vice versa.
    fn transformation_matrix(&self,
        theta_rad:f64,
        phi_rad:f64)-> na::Matrix3<f64>{
            let sintheta = theta_rad.sin();
            let costheta = theta_rad.cos();
            let sinphi = phi_rad.sin();
            let cosphi = phi_rad.cos();


            let t_matrix = na::Matrix3::new(
                sintheta*costheta, costheta*cosphi, -sinphi,
                sintheta*sinphi, costheta*cosphi, cosphi,
                costheta, -sintheta, 0.0,
            );//spherical to cartesian matrix

            match self{
                Coordinates::Spherical(_) => {t_matrix}
                Coordinates::Cartesian(_) => {t_matrix.transpose()}
            }
        }
    /// This method transforms coordinates from spherical to cartesian and vice versa. 
    /// 
    /// ### Arguments:
    /// * `theta_rad` - The colatitude angle in radians.
    /// * `phi_rad` - The azimuthal angle in radians.
    /// 
    /// ### Returns:
    /// * ` ` - a new instance of [Coordinates] with different base.
    pub fn transform(&self,theta_rad:f64,phi_rad:f64)->Coordinates{
        let t_matrix = self.transformation_matrix(theta_rad, phi_rad);
        match self{
            Coordinates::Cartesian(value)=>{
                Coordinates::Spherical(t_matrix * value)}
            Coordinates::Spherical(value)=>{
                Coordinates::Cartesian(t_matrix * value)}
        }
    }
    
    /// This method obtains the r component of the spherical coordinates
    /// 
    /// ### Arguments: 
    /// * `self` - A [Coordinates] instance that should be Spherical variant.
    /// 
    /// ### Returns:
    /// * `Option` - This method returns a `Some(f64)` that has the binded the r component value, in case it was casted by a sperical variant, or a  `None` otherwise.
    pub fn r_component(&self)->Option<f64>{
        match self{
            Coordinates::Spherical(value)=>{
                Some(value[0])
            }
            _=>{None}
        }
    }
}


/// Compute the Lagrangian displacement vector in spherical coordinates
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// * `sintheta` - sine of the colatitude coordinate (theta in rads)
/// * 'costheta' - cosine of the colatitude coordinate (theta in rads)
/// * `phi`   - azimuthal coordinate  in rads
/// * `radial_amplitude`     - amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::joris_math::spherical_harmonics::norm_factors])
/// * `tangential_amplitude` - amplitude in the tangential direction times the normalization factor  'Y_l^m' (see [temp_name_lib::joris_math::spherical_harmonics::norm_factors])
/// 
/// ### Returns:
/// This function can return an [Ok] or [Err] variants of [Result] that will have the following values binded to them:
/// * `Ok(Coordinates::Spherical)` - an Ok  variant that has binded the spherical components of the displacement vector in the`r,θ,φ` order.
/// * `Err(DivisionByZero)` - an Err variant that has binded the error produced if the colatitude  coordinate (theta) is too small.
pub fn ddisplacement(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64)->Result<Coordinates,MathErrors>{
        match sintheta.abs() <= f64::EPSILON.sqrt(){
            true => {Err(MathErrors::DivisionByZero)}
            false => {
                let phase = mode.phase_offset;
                let l = mode.l;
                let m =  mode.m;
                
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

/// This module contains the functions to calculate the derivatives of the lagrangian displacement vector over 
/// the surface of a star using spherical coordinates.
mod displacement_derivatives;


/// This function computes the spherical components of the surface normal vector on a reference frame where the z-axis 
/// coincides with the rotation axis. A surface normal is a vector which stands locally perpendicular to the surfaces and
/// has a length equal to the area of the local surface cell. 
/// 
/// *WARNING* This function does NOT multiply each component with 
/// R_0^2 dθdφ. Users must do this themselves. 
/// 
/// ### Arguments:
/// * `parameters` - The data contained in [PPulstarConfig], here you find the parameters that describe the pulsation modes and the star.
/// * `theta_rad` - The colatitude angle in rads, must not be too small in order to avoid the poles.
/// * `phi_rad` - The azimuthal angle in rads
pub fn ssurface_normal(
parameters: &PPulstarConfig,
theta_rad: f64,
phi_rad: f64
)->Result<Coordinates,MathErrors>{
    let sintheta = theta_rad.sin();
    let costheta = theta_rad.cos();

    let mut total_p_ds = Coordinates::Spherical(na::Vector3::new(0.0,0.0,0.0));
	let mut total_dev1=0.0;
	let mut total_dev2=0.0;
	let mut total_dev3=0.0;
	let mut total_dev4=0.0;
    for mode in parameters.mode_data.iter(){

        let radial_amplitude = ampl_r(mode);
        let tangential_amplitude = ampl_t(mode);
        
        let pulsation_displacement = ddisplacement(
            mode,
            sintheta,
            costheta, 
            phi_rad, 
            radial_amplitude, 
            tangential_amplitude)?;
        
        let drdtheta = displacement_derivatives::d_dr_rdtheta(
            mode,
            sintheta, 
            costheta, 
            phi_rad);
        
        let drdphi = displacement_derivatives::d_dr_rdphi(
            mode,
            sintheta,
            costheta,
            phi_rad);
        
        let dtdtheta = displacement_derivatives::d_dtheta_dtheta(
            mode,
            sintheta,
            costheta,
            phi_rad);
        
        let dpdphi = displacement_derivatives::d_dphi_dphi(
            mode,
            sintheta,
            costheta,
            phi_rad)?;//<- the ? is necessary to pass to the calling function if the colatitude angle theta is too close to the poles.
        
        
        total_p_ds += pulsation_displacement;

        total_dev1 += drdtheta;
        total_dev2 += drdphi;
        total_dev3 += dtdtheta;
        total_dev4 += dpdphi;
    } 
    let mut surface_normal_coords = na::Vector3::new(0.0,0.0,0.0);
    if let Coordinates::Spherical(coords) = total_p_ds{//<- This is an idiomatic way to extract values of an enum variant in Rust
        let r_hat = 1.0;
        let theta_hat = coords[1]-total_dev1;
        let phi_hat = sintheta * coords[2] - total_dev2/sintheta;

        let normal = na::Vector3::new(r_hat,theta_hat,phi_hat);
        let length = sintheta * (1.0 + 2.0 * coords[0]
            + costheta/sintheta * coords[1]
            + total_dev3 + total_dev4);
        
        surface_normal_coords = normal * length;

    }
    Ok(Coordinates::Spherical(surface_normal_coords))

}


///Given the components of both the surface normal vector and the unit vector
///pointing towards the observer, this function gives the cosine of the angle
///between those two vectors. 
/// 
/// The z-axis coincides with the rotation axis.
/// 
/// ### Arguments:
/// * `surface_normal` - a spherical [Coordinates] vector normal to a surface cell that has as length the (normalized) area of the cell. See [ssurface_normal]
/// * `k` - a unit vector pointing towards the observer on a frame of reference where the z-axis coincides with the rotation axis.
/// * `theta_rad` - the colatitude angle on the star in radians. This angle should not be too small. 
/// * `phi_rad` - the azimuthal angle on the star in radians.
/// 
/// ### Returns:
/// * `cos_chi` - a `f64` value that is the cosine of the angle between `surface_normal` and `k`
pub fn cos_chi(surface_normal:&Coordinates,
    k: &Coordinates,theta_rad: f64,phi_rad:f64)->f64{
    match k {
        Coordinates::Cartesian(_)=>{
            let k_spherical = k.transform(theta_rad, phi_rad);
            k_spherical.project_vector(surface_normal).unwrap()
            /(k_spherical.vector_length()*surface_normal.vector_length())
        }
        Coordinates::Spherical(_)=>{
            k.project_vector(surface_normal).unwrap()
            /(k.vector_length()*surface_normal.vector_length())
        }
    }
}


/// This function calculates the amplitude of the relative radial displacement multiplied by the normalization factor `Y_l^m`
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// 
/// ### Returns:
/// * `radial_amplitude` - A `f64` value that contains the amplitude of relative radial displacement (thus without units) caused by the pulsations of a given mode. 
pub fn ampl_r(mode:&PulsationMode)->f64{
    mode.rel_dr * ylmnorm(mode.l, mode.m)
}

/// This function calculates the amplitude of the relative tangential displacement multiplied by the normalization factor `Y_l^m`
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// 
/// ### Returns:
/// * `tangential_amplitude` - A `f64` value that contains the amplitude of relative tangential displacement (thus without units) caused by the pulsations of a given mode. 
pub fn ampl_t(mode:&PulsationMode)->f64{
    mode.rel_dr * ylmnorm(mode.l, mode.m)*mode.k
}