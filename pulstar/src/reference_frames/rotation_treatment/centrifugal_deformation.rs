use std::ops::{AddAssign,Mul}; 
use crate::PulsationMode;
use crate::reference_frames::rotation_treatment::non_rotating::
{non_rotating_d_dr_rdtheta, non_rotating_d_dtheta_dtheta, 
non_rotating_displacement, non_rotating_d_dr_rdphi, non_rotating_d_dphi_dphi};
use crate::reference_frames::{*,Coordinates};
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
pub fn deformed_displacement(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64)->Result<Coordinates,MathErrors>{

    let displacement_func = 
    |mode:&PulsationMode,
    sintheta:f64,
    costheta: f64,
    phi:f64,
    radial_amplitude: f64,
    tangential_amplitude: f64,|{non_rotating_displacement(mode, sintheta, costheta, phi, radial_amplitude, tangential_amplitude)};

    let deformed_displacement = add_deformations_generic(displacement_func,
    mode, sintheta, costheta, phi, radial_amplitude, tangential_amplitude);

    deformed_displacement
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
pub fn deformed_d_dr_rdtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{



    let disp_function = 
    |mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64, _dummy_parameter1:f64,_dummy_parameter2:f64| {Ok(non_rotating_d_dr_rdtheta(mode,sintheta,costheta,phi))};
    
    if let Ok(return_value) = add_deformations_generic(disp_function,
         mode, sintheta, costheta, phi, 0.0, 0.0){return_value}else{panic!("unable to compute derivative of radial displacement against theta")}
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
pub fn deformed_d_dtheta_dtheta(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{
    
    let disp_function =|mode: &PulsationMode,sintheta: f64,costheta: f64,phi: f64,_dummy1:f64,_dummy2:f64|
    {Ok(non_rotating_d_dtheta_dtheta(mode, sintheta, costheta, phi))};
    
    if let Ok(return_value) = add_deformations_generic(disp_function,
         mode, sintheta, costheta, phi, 0.0, 0.0){return_value}else{panic!("unable to compute derivative of dtheta against theta")}
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
pub fn deformed_d_dr_rdphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> f64{

    let disp_function =|mode: &PulsationMode,sintheta: f64,costheta: f64,phi: f64,_dummy1:f64,_dummy2:f64|
    {Ok(non_rotating_d_dr_rdphi(mode, sintheta, costheta, phi))};
    
    if let Ok(return_value) = add_deformations_generic(disp_function,
         mode, sintheta, costheta, phi, 0.0, 0.0){return_value}else{panic!("unable to compute derivative of dr/r against phi")}
    
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
pub fn deformed_d_dphi_dphi(
    mode: &PulsationMode,
	sintheta: f64,
	costheta: f64,
	phi: f64) -> Result<f64,MathErrors>{

    let disp_function =|mode: &PulsationMode,sintheta: f64,costheta: f64,phi: f64,_dummy1:f64,_dummy2:f64|
    {non_rotating_d_dphi_dphi(mode, sintheta, costheta, phi)};
    
    add_deformations_generic(disp_function,
         mode, sintheta, costheta, phi, 0.0, 0.0)
}



fn compute_index_from_mode(l:u16)->usize{
    if l%2==0{
        (l/2) as usize
    }else{
        ((l-1)/2) as usize
    }
}


pub fn add_deformations_generic<F,T>(
    function_to_eval:F,
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64)
    ->Result<T,MathErrors>
    where
    F:Fn(&PulsationMode,f64,f64,f64,f64,f64)->Result<T,MathErrors>,
    T: AddAssign<T>, 
    f64: Mul<T,Output = T>,
{
    let coeff_expansion = match &mode.rotation_effects{
            RotationRegime::CentrifugalDeformation { coefficient_expansion:vecs }=>{vecs.clone()}
            _=>{return Err(MathErrors::RequestUnrelatedRotationRegime)}
        };
        if coeff_expansion.len()>5{return Err(MathErrors::OrderOfExpansionNotSupported)}
        else{
            // compute the displacement of the main mode of pulsation
            let main_coeff_index=compute_index_from_mode(mode.l);
            let mut first_contribution =  coeff_expansion[main_coeff_index] * function_to_eval(
                mode, 
                sintheta, costheta, phi, radial_amplitude, tangential_amplitude)?;
            
            //add the crossed terms
            for (index,coeff) in coeff_expansion.iter().enumerate(){
                if index != main_coeff_index{
                    let mut temp_pulsation_mode = mode.clone();
                    //extract l from index
                    let odd_or_even = mode.l%2;
                    let temp_l = odd_or_even + 2 * index as u16;
                    if mode.m.abs() as u16 > temp_l{continue};
                    temp_pulsation_mode.l =temp_l;
                    let temp_contribution = *coeff * function_to_eval(&temp_pulsation_mode,
                        sintheta, costheta, phi, radial_amplitude, tangential_amplitude)?;
                    first_contribution += temp_contribution;
                }
            }

        Ok(first_contribution)

    }
}
/*fn add_deformations_coordinates<F>(
    function_to_eval:F,
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64)
    ->Result<Coordinates,MathErrors>
    where
    F:Fn(&PulsationMode,
    f64,f64,f64,f64,f64)->Result<Coordinates,MathErrors>
{
    let coeff_expansion = match &mode.rotation_effects{
            RotationRegime::CentrifugalDeformation { coefficient_expansion:vecs }=>{vecs.clone()}
            _=>{return Err(MathErrors::RequestUnrelatedRotationRegime)}
        };
        if coeff_expansion.len()>5{return Err(MathErrors::OrderOfExpansionNotSupported)}
        else{
            // compute the displacement of the main mode of pulsation
            let main_coeff_index=compute_index_from_mode(mode.l);
            let mut first_contribution = coeff_expansion[main_coeff_index] *function_to_eval(
                mode, 
                sintheta, costheta, phi, radial_amplitude, tangential_amplitude)?;
            
            //add the crossed terms
            for (index,coeff) in coeff_expansion.iter().enumerate(){
                if index != main_coeff_index{
                    let mut temp_pulsation_mode = mode.clone();
                    //extract l from index
                    let odd_or_even = mode.l%2;
                    let temp_l = odd_or_even + 2 * index as u16;
                    if mode.m.abs() as u16 > temp_l{continue};
                    temp_pulsation_mode.l =temp_l;
                    let temp_contribution = coeff_expansion[index]*function_to_eval(&temp_pulsation_mode,
                        sintheta, costheta, phi, radial_amplitude, tangential_amplitude)?;
                    first_contribution += *coeff * temp_contribution;
                }
            }

            Ok(first_contribution)

        }
}*/
