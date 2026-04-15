use crate::reference_frames::rotation_treatment::tar::TARCollection;
use crate::{PulsationMode, PulstarConfig};
use crate::reference_frames::{Coordinates,displacement,ampl_r,ampl_t};

/// This function calculates the local temperature and log_g ver a surface cell
/// 
/// ### Arguments:
/// * `parameters` - The data contained in [PulstarConfig], here you find the parameters that describe the pulsation modes and the star.
/// * `theta` -  The colatitude angle θ in radians
/// * 'dtheta' - Angular displacement in the colatitude coordinate from one neighbouring [SurfaceCell] to the next, it has the same units as θ.
/// * `phi_rad` - The azimuthal angle in radians
/// * `g0` - The base value of local gravity calculated as 10^(log_g0) on the surface of the star.
/// * `T0` - The base value of the effective temperature  on the surface of the star. 
/// * `tar_functions` - An [Option] enum that has the following variants:
///     - [Some] variant that has binded a reference to a [TARCollection]
///     - [None] in case tar functions are not needed.
/// ### Returns:
/// * `(local_temperature,local_logg)` - a tupple containing the local effective temperature and the local value of log_g
pub fn local_surface_temperature_logg(
    parameters:&PulstarConfig,
    theta:f64,
    dtheta:f64,
    phi:f64,
    g0:f64,
    temperature_0:f64,
    tar_collections:&[Option<TARCollection>],
    )->(f64,f64){
    let mut local_temperature = 0.0;
    let mut local_g = 0.0;

    for (index,mode) in parameters.mode_data.iter().enumerate(){
        //Check if it's not a trivial case
        if mode.rel_dg != 0.0 || mode.rel_dtemp != 0.0 {
            let radial_amplitude = ampl_r(mode);
            let tangential_amplitude = ampl_t(mode);
            
            if mode.rel_dtemp != 0.0 {
            let ds = local_variable_pulsation_displacement(
            mode, 
            theta, 
            dtheta, 
            phi, 
            radial_amplitude, 
            tangential_amplitude, 
            mode.phase_temp,
            &tar_collections[index]);
            if let Some(ds_r) = ds.r_component(){
            local_temperature += mode.rel_dtemp * ds_r;
            };
            }
        
            if mode.rel_dg !=0.0 {
            let ds = local_variable_pulsation_displacement(
            mode, 
            theta, 
            dtheta, 
            phi, 
            radial_amplitude, 
            tangential_amplitude, 
            mode.phase_logg,
            &tar_collections[index]);
            
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
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
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
    theta:f64,
    dtheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64,
    dif_phase:f64,
    tar_functions:&Option<TARCollection>
    )->Coordinates{
    
    //[Ricardo:] There's a shorter version of this, namely 
    let mut mode_with_dif_phase = mode.clone();
     mode_with_dif_phase.phase_offset += dif_phase;//<-- Here is where we ad the phase difference

    displacement(
        &mode_with_dif_phase,
        theta,
        dtheta, 
        phi, 
        radial_amplitude, 
        tangential_amplitude,
        tar_functions).unwrap()

}