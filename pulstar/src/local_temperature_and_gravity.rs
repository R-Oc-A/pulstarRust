use crate::{PPulstarConfig, PulsationMode};
use crate::reference_frames::{Coordinates,ddisplacement,ampl_r,ampl_t};

/// This function calculates the local temperature and log_g ver a surface cell
/// 
/// ### Arguments:
/// * `parameters` - The data contained in [PPulstarConfig], here you find the parameters that describe the pulsation modes and the star.
/// * `theta_rad` -  The colatitude angle in radians
/// * `phi_rad` - The azimuthal angle in radians
/// * `g0` - The base value of local gravity calculated as 10^(log_g0) on the surface of the star.
/// * `T0` - The base value of the effective temperature  on the surface of the star. 
/// 
/// ### Returns:
/// * `(local_temperature,local_logg)` - a tupple containing the local effective temperature and the local value of log_g
pub fn local_surface_temperature_logg(
    parameters:&PPulstarConfig,
    theta_rad:f64,
    phi_rad:f64,
    g0:f64,
    temperature_0:f64,
    )->(f64,f64){
    let mut local_temperature = 0.0;
    let mut local_g = 0.0;

    let sintheta = theta_rad.sin();
    let costheta = theta_rad.cos();

    for mode in parameters.mode_data.iter(){
        //Check if it's not a trivial case
        if mode.rel_dg != 0.0 || mode.rel_dtemp != 0.0 {
            let radial_amplitude = ampl_r(mode);
            let tangential_amplitude = ampl_t(mode);
            
            if mode.rel_dtemp != 0.0 {
            let ds = local_variable_pulsation_displacement(
            mode, 
            sintheta, 
            costheta, 
            phi_rad, 
            radial_amplitude, 
            tangential_amplitude, 
            mode.phase_rel_dtemp);
            if let Some(ds_r) = ds.r_component(){
            local_temperature += mode.rel_dtemp * ds_r;
            };
            }
        
            if mode.rel_dg !=0.0 {
            let ds = local_variable_pulsation_displacement(
            mode, 
            sintheta, 
            costheta, 
            phi_rad, 
            radial_amplitude, 
            tangential_amplitude, 
            mode.phase_rel_dg);
            if let Some(ds_r) = ds.r_component(){
            local_g += mode.rel_dg * ds_r; 
            };
            }
        }
    }

    local_g += 1.0;
    local_temperature += 1.0;

    local_g *= g0;
    local_temperature *= temperature_0;

    let local_logg = local_g.log10();
    (local_temperature,local_logg)
}

/// This function calculates variations on the pulsation displacement due to diferent phase of some either temperature or log_g
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PPulstarConfig]
/// * `sintheta` - sine of the colatitude coordinate (theta in rads)
/// * 'costheta' - cosine of the colatitude coordinate (theta in rads)
/// * `phi_rad`   - azimuthal coordinate  in rads
/// * `radial_amplitude`     - Amplitude in the radial direction times the normalization factor `Y_l^m`(see [temp_name_lib::joris_math::spherical_harmonics::norm_factors])
/// * `tangential_amplitude` - Amplitude in the tangential direction times the normalization factor  'Y_l^m' (see [temp_name_lib::joris_math::spherical_harmonics::norm_factors])
/// * `dif_phase` - Phase diference in the observed quantity. Could be given by either temperature or `log_g`
/// ### Returns:
/// `Coordinates::Spherical(r,θ,φ)` - [Coordinates] in spherical basis with the pulsation displacement.
fn local_variable_pulsation_displacement(
    mode: &PulsationMode,
    sintheta:f64,
    costheta:f64,
    phi_rad:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64,
    dif_phase:f64)->Coordinates{
    
    //[Ricardo:] There's a shorter version of this, namely 
    //
    // let mut mode_with_dif_phase = *mode.clone();
    // mode_with_dif_phase.phase_offset += dif_phase;
    //
    // However I'm not quite sure I want to implement the Clone and Copy traits on the PulsationMode structure. 
    let mode_with_dif_phase = PulsationMode { 
        l: mode.l,
        m: mode.m,
        rel_dr: mode.rel_dr,
        k: mode.k,
        frequency: mode.frequency,
        phase_offset: mode.phase_offset + dif_phase,//<-- Here is where we add the phase difference
        rel_dtemp:mode.rel_dtemp,
        phase_rel_dtemp: mode.phase_rel_dtemp,
        rel_dg: mode.rel_dg,
        phase_rel_dg: mode.phase_rel_dg};

    ddisplacement(
        &mode_with_dif_phase,
        sintheta,
        costheta, 
        phi_rad, 
        radial_amplitude, 
        tangential_amplitude).unwrap()

}