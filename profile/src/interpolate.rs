use temp_name_lib::star_physics::local_values::temperature;

use super::{*,intensity::*};

/// Computes by linear interpolation the toal and continuum intensity
/// Input: 
/// Output: (intens_result,cont_result) where intens is the total intensity; cont is the continuum intensity
pub fn interpolate(
    temperature_vec:&[f64],
    log_gravity_vec:&[f64],
    shifted_wavelengths:&[f64],
    flux_vec:&Vec<Vec<f64>>,
    cont_vec:&Vec<Vec<f64>>,
    wavelengths_vec:&Vec<Vec<f64>>,  
    temperature:f64,
    log_g:f64)->(Vec<f64>,Vec<f64>){
    let mut counter:usize = 0;
    let t0 = (temperature-temperature_vec[0])/(temperature_vec[1]-temperature_vec[0]);
    let t1 = (log_g-log_gravity_vec[0])/(log_gravity_vec[1]-log_gravity_vec[0]);

    let mut interpolated_flux = vec![0.0;shifted_wavelengths.len()];
    let mut interpolated_cont= vec![0.0;shifted_wavelengths.len()];
    for (i,temperature) in temperature_vec.iter().enumerate(){
        let delta_temp = 1.0-t0 + (2.0*t0 -1.0) * i as f64;
        for (j,logg) in log_gravity_vec.iter().enumerate(){
            let delta_logg = 1.0-t1 + (2.0 * t1 -1.0) * j as f64;
            for (n,lambda) in shifted_wavelengths.iter().enumerate(){
                let index = search_geq(&wavelengths_vec[counter], *lambda);
                let t2 = (lambda - wavelengths_vec[counter][index-1])/
                (wavelengths_vec[counter][index]-wavelengths_vec[counter][index-1]);
                for k in 0..=1{
                    let delta_lambda = 1.0-t2 + (2.0 * t2 - 1.0) * k as f64;
             
                    interpolated_flux[n] += 
                    interpolate_by_element(delta_temp,
                        delta_logg,
                        delta_lambda,
                        flux_vec[counter][index-1]);
             
                    interpolated_cont[n] += 
                    interpolate_by_element(delta_temp,
                        delta_logg,
                        delta_lambda,
                        cont_vec[counter][index-1]);
                }
            }
            counter += 1;
        }
    }

    (interpolated_flux,interpolated_cont)    
}


fn interpolate_by_element(delta_temp:f64,
delta_logg:f64,
delta_lamb:f64,
value:f64)->f64{
    delta_temp*delta_logg*delta_lamb*value
}