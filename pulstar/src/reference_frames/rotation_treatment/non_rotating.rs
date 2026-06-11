use spherical_harmonics::{plmcos::plmcos};
use spherical_harmonics::{norm_factor::ylmnorm};
use spherical_harmonics::{
    d_plmcos_dtheta::{
        deriv1_plmcos_dtheta as d_plmcos_dtheta,
        deriv2_plmcos_dtheta as d2_plmcos_dtheta}}; 
use crate::PulsationMode;
use crate::reference_frames::{*,na,Coordinates};
use temp_name_lib::utils::MathErrors;




/// Compute the Lagrangian displacement vector in spherical coordinates
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude coordinate (theta in rads)
/// * 'costheta' - cosine of the colatitude coordinate (theta in rads)
/// * `phi`   - azimuthal coordinate  in rads
/// * `radial_amplitude`     - amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::math_module::spherical_harmonics::norm_factors])
/// * `tangential_amplitude` - amplitude in the tangential direction times the normalization factor  'Y_l^m' (see [temp_name_lib::math_module::spherical_harmonics::norm_factors])
/// 
/// ### Returns:
/// This function can return an [Ok] or [Err] variants of [Result] that will have the following values binded to them:
/// * `Ok(Coordinates::Spherical)` - an Ok  variant that has binded the spherical components of the displacement vector in the`r,θ,φ` order.
/// * `Err(DivisionByZero)` - an Err variant that has binded the error produced if the colatitude  coordinate (theta) is too small.
pub fn non_rotating_displacement(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64)->Result<Coordinates,MathErrors>{
        match sintheta.abs() <= f64::EPSILON.sqrt(){
            true => {Err(MathErrors::DivisionByZero)}
            false => {
                let phase = mode.phase;
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


    
///Computes the derivatives of Δr/r0 with respect to θ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to θ
pub fn non_rotating_d_dr_rdtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let r_dr = mode.rel_dr;
    let phase= mode.phase;
    let l = mode.l;
    let m= mode.m;
                            
    r_dr*ylmnorm(l,m)
    * d_plmcos_dtheta(l,m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).cos()
}

///Computes the derivatives of Δθ with respect to θ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the displacement in θ with respect to θ
pub fn non_rotating_d_dtheta_dtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let r_dr = mode.rel_dr;
    let phase= mode.phase;
    let k = mode.k;
    let l = mode.l;
    let m= mode.m;

    r_dr*ylmnorm(l,m)*k
    * d2_plmcos_dtheta(l,m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).cos()
}

///Computes the derivatives of Δr/r0 with respect to φ in the point with spherical
///coordinates θ,φ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// * an `f64` - This value is the derivative of the relative radial displacement with respect to φ 
pub fn non_rotating_d_dr_rdphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let r_dr = mode.rel_dr;
    let phase= mode.phase;
    let l = mode.l;
    let m= mode.m;
    
    r_dr * ylmnorm(l,m) * (-m as f64)
    * plmcos(l, m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).sin()
}

///Computes the derivatives of Δϕ with respect to ϕ in the point with spherical
///coordinates θ,ϕ
/// ### Arguments: 
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `sintheta` - sine of the colatitude angle (theta in rads)
/// * `costheta` - cosine of the colatitude angle (theta in rads)
/// * `phi` - azimuthal coordinate in rads
/// ### Returns:
/// This function returns a [Result] with the following variants:
/// * `Ok(f64)` - Where the binded value is the derivative of the displacement in φ with respect to φ 
/// * `Err(DivisionByZero)` - Where the binded error is returned to the calling function and indicates that the theta value was too small.
pub fn non_rotating_d_dphi_dphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> Result<f64,MathErrors>{

    match sintheta.abs() < MACHINE_PRECISION{  
        false => {
        let r_dr = mode.rel_dr;
        let phase= mode.phase;
        let k= mode.k;
        let l = mode.l;
        let m= mode.m;

        Ok(r_dr * k * ylmnorm(l, m) * (-(m as f64).powi(2))
        * plmcos(l, m.abs() as u16, sintheta, costheta)
        * (phase + (m as f64) * phi).cos()
        /(sintheta.abs().powi(2)) )
        }

        true =>{
        Err(MathErrors::DivisionByZero) //will pass the error in order for the calling function to do something
        }
    }
}