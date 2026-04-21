use spherical_harmonics::{
    {plmcos::plmcos},
    {d_plmcos_dtheta::{deriv1_plmcos_dtheta as d_plmcos_dtheta}} };
use crate::PulsationMode;
use crate::reference_frames::rotation_treatment::non_rotating::{non_rotating_d_dphi_dphi, non_rotating_d_dtheta_dtheta};
use crate::reference_frames::{*,na,Coordinates};
use temp_name_lib::utils::MathErrors;



pub fn amplitude_lp1(radial_amplitude:f64,
    spin_parameter:f64,
    l:f64,
    m:f64,
    k:f64,)->f64{    
    let lp1 = l+1.0;

    radial_amplitude * 2.0 * spin_parameter
    * (lp1 - m)/lp1 * 2.0/(2.0 * l + 1.0)
    * (1.0 - l as f64 * k)
}

pub fn amplitude_lm1(radial_amplitude:f64,
    spin_parameter:f64,
    l:f64,
    m:f64,
    k:f64,)->f64{    

    radial_amplitude * 2.0 * spin_parameter
    * (l + m)/l * 2.0/(2.0*l + 1.0)
    * (1.0 + (l + 1.0)* k)
}

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
pub fn perturbative_displacement(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64,
    spin_parameter:f64)->Result<Coordinates,MathErrors>{
        match sintheta.abs() <= f64::EPSILON.sqrt(){
            true => {Err(MathErrors::DivisionByZero)}
            false => {
                if mode.l == 0 {Err(MathErrors::NotDefinedForRadialOrSectorialPulsations)}
                else{
                let spheroidal_displacement = reference_frames::rotation_treatment::non_rotating::
                non_rotating_displacement(mode, sintheta, costheta, phi, radial_amplitude, tangential_amplitude)?;

                let phase = mode.phase;
                let l = mode.l;
                let lp1 = l + 1;
                let lm1 = l-1;
                let m =  mode.m.abs() as u16;
                let k = mode.k;
                
                let plp1m = plmcos(lp1,m, sintheta, costheta); 
                let dplp1m = d_plmcos_dtheta(lp1, m, sintheta, costheta);
                let plm1m = plmcos(lm1, m, sintheta, costheta); 
                let dplm1m = d_plmcos_dtheta(lm1, m, sintheta, costheta);
                
                //First toroidal term. Taken from Zima 2008 and FAMIAS user manual. 
                let amplitude = amplitude_lp1(radial_amplitude,
                    spin_parameter, l as f64, m as f64, k as f64);

                let delta_r     = 0.0;
                let delta_theta = -amplitude/sintheta * plp1m * (mode.m as f64)
                                    * (phase + 0.5*PI + (mode.m as f64)*phi).sin();
                let delta_phi   = - amplitude * dplp1m 
                                    * (phase + 0.5*PI + (mode.m as f64)*phi).cos();
                
                let first_toroidal_displacement = Coordinates::Spherical(
                    na::Vector3::new(delta_r,delta_theta,delta_phi));
                
                //Second toroidal term. Taken from Zima 2008 and FAMIAS user manual. 
                let amplitude = amplitude_lm1(radial_amplitude,
                    spin_parameter, l as f64, m as f64, k as f64);

                let delta_r     = 0.0;
                let delta_theta = -amplitude/sintheta * plm1m * (mode.m as f64)
                                    * (phase + 0.5*PI + (mode.m as f64)*phi).sin();
                let delta_phi   = - amplitude * dplm1m 
                                    * (phase + 0.5*PI + (mode.m as f64)*phi).cos();
                
                let second_toroidal_displacement = Coordinates::Spherical(
                    na::Vector3::new(delta_r,delta_theta,delta_phi));
                
                Ok(((spheroidal_displacement
                   + first_toroidal_displacement)?
                   + second_toroidal_displacement )? )
                
                }
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
pub fn perturbative_d_dr_rdtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{
    
    rotation_treatment::non_rotating::non_rotating_d_dr_rdphi(mode, sintheta, costheta, phi)
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
pub fn perturbative_d_dtheta_dtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64,
    spin_parameter:f64) -> Result<f64,MathErrors>{
    match sintheta > f64::EPSILON.sqrt(){
        true =>{    
            let spheroidal_part = non_rotating_d_dtheta_dtheta(mode, sintheta, costheta, phi);
        
            // Computations of the toroidal parts.
            let r_dr = ampl_r(mode);
            let phase= mode.phase;
            let k = mode.k;
            let l = mode.l;
            let m= mode.m;
            
            let lp1 = l+1;
            let lm1 = l-1;
        
            let plp1m = plmcos(lp1, m.abs() as u16, sintheta, costheta); 
            let dplp1m = d_plmcos_dtheta(lp1, m.abs() as u16, sintheta, costheta);
            let plm1m = plmcos(lm1, m.abs() as u16, sintheta, costheta); 
            let dplm1m = d_plmcos_dtheta(lm1, m.abs() as u16, sintheta, costheta);
            
            // derivative of the first toroidal part
            let amplitude = amplitude_lp1(r_dr,
                 spin_parameter, l as f64, m as f64, k);
        
            let first_toroidal_part = -amplitude * (mode.m as f64) * (phase + 0.5*PI + (mode.m as f64)*phi).sin()*
                (-costheta/(sintheta.powi(2)) * plp1m + 1.0/sintheta * dplp1m);
        
            // derivative of the second toroidal part
            let amplitude = amplitude_lm1(r_dr, 
                spin_parameter, 
                l as f64, m as f64, k);
            
            let second_toroidal_part = amplitude * (mode.m as f64) * (phase + 0.5*PI + (mode.m as f64)*phi).sin()*
                (-costheta/sintheta.powi(2) * plm1m + 1.0/sintheta * dplm1m);
        
            Ok(spheroidal_part + first_toroidal_part + second_toroidal_part)
        }
            false => { Err(MathErrors::DivisionByZero)}
        }
    
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
pub fn perturbative_d_dr_rdphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{
    
    rotation_treatment::non_rotating::non_rotating_d_dr_rdphi(mode,
         sintheta, costheta, phi)
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
pub fn perturbative_d_dphi_dphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64,
    spin_parameter:f64) -> Result<f64,MathErrors>{

        let spheroidal_part = non_rotating_d_dphi_dphi(mode, sintheta, costheta, phi)?;

        // Computations of the toroidal parts.
        let r_dr = ampl_r(mode);
        let phase= mode.phase;
        let k = mode.k;
        let l = mode.l;
        let m= mode.m;
        
        let lp1 = l+1;
        let lm1 = l-1;
    
        let dplp1m = d_plmcos_dtheta(lp1, m.abs() as u16, sintheta, costheta);
        let dplm1m = d_plmcos_dtheta(lm1, m.abs() as u16, sintheta, costheta);
        
        
        // derivative of the first toroidal part
        let amplitude = amplitude_lp1(r_dr,
             spin_parameter, l as f64, m as f64, k);
    
        let first_toroidal_part = amplitude *dplp1m 
        * (mode.m as f64) * (phase + 0.5*PI + (mode.m as f64)*phi).sin();
        
        // derivative of the second toroidal part
        let amplitude = amplitude_lm1(r_dr,
             spin_parameter, l as f64, m as f64, k);
    
        let second_toroidal_part = amplitude *dplm1m 
        * (mode.m as f64) * (phase + 0.5*PI + (mode.m as f64)*phi).sin();
        
        Ok(spheroidal_part + first_toroidal_part + second_toroidal_part)

}