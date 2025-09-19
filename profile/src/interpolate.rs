use super::*;

impl SurfaceCell{
    pub fn interpolate(& self, grids_data: &mut GridsData,global_flux:&mut FluxOfSpectra){
    

    // Counter orders the intensities 
    let mut counter:usize = 0;
    
    grids_data.extract_important_rows(global_flux);
    let t0 = (self.t_eff-grids_data.temperature_vector[0])/
        (grids_data.temperature_vector[2]-grids_data.temperature_vector[0]);
    let t1 = (self.log_g-grids_data.log_g_vector[0])/(grids_data.log_g_vector[1]-grids_data.log_g_vector[0]);
    
    // loop over the four elements of the (temperature, log_gravity) parameter space
    for i in 0..=1{//<- i goes from 0 to 1
        let delta_temp = match i {
            0 => 1.0-t0, 
            1 => t0,
            _ =>{0.0}};//although verbose, computationally this is faster than multiplying floats

        for j in 0..=1{//<- j goes from 0 to 1
            let delta_logg = match j{
                0 => {1.0-t1}
                1 => {t1}
                _=> {0.0}};

            
            for (n,lambda) in global_flux.shifted_wavelength.iter().enumerate(){//<-this is to carry the same process onto every wavelength
                
                let index = grids_data.row_indices[2*n];
                
                let t2 = (lambda - grids_data.grid_wavelengths[index])/
                (grids_data.grid_wavelengths[index+1]-grids_data.grid_wavelengths[index]);
                
                for k in 0..=1{
                    let delta_lamb = match k {
                        0 => {1.0-t2}
                        1 => {t2}
                        _ => {0.0} };
                        
                        global_flux.flux[n] += get_contribution_from_vertex(delta_temp, 
                            delta_logg, 
                            delta_lamb, 
                            grids_data.flux[counter][index]) * self.area;
                        
                        global_flux.continuum[n] += get_contribution_from_vertex(delta_temp,
                            delta_logg,
                            delta_lamb,
                            grids_data.continuum[counter][index]) * self.area;
                }
            }
            counter += 1;//<- counter should from 0 to 3
        }
    }
        
    }
}
/// Computes by linear interpolation (see [trilinear interpolation](https://paulbourke.net/miscellaneous/interpolation/) ) the specific and continuum intensity for a surface cell.
/// 
/// ### Arguments:
/// * `intensity_dfs` - a reference to an instance of [IntensityDataFrames] that contains all of the relevant information parsed from the intensity grid files. 
/// * `relevant_indices` - a [Vec] collection of [usize] that serves as keys to use only the relevant information of the DataFrames.
/// * `shifted_wavelengths` - a vector that contains the doppler shifted observed wavelengths
///  The way the values are ordered is so that  
/// * `shifetd_wavelengths[i]<shifetd_wavelengths[i+1]`
/// 'flux_vec' - a vector of vectors, it contains four intensity flux vectors that were calculated from the relevant intensity grid files
/// 'cont_vec' - a vector of vectors, it contains four continuum flux vectors that were calculated from the relevant intensity grid files
/// * `wavelengths_vec` - a vector of vectors, it contains the four wavelength vectors associated with the intensity fluxes
/// * `temperature` - the temperature over a surface cell of the rasterized star.
/// * `log g` - the log g value over a surface cell of the rasterized star.
///  ### Returns:
/// * `(intens_result,cont_result)` - a tuple that contains two `f64` vectors where intens is the total intensity; cont is the continuum intensity
/// * `intens_result` - is the collection of specific intensity fluxes that correspond to the observed wavelengths
/// * `cont_result` - is the collection of the continuum intensity fluxes that correspond to the observed wavelengths
pub fn interpolate(
    intensity_dfs:&IntensityDataFrames,
    relevant_indices:Vec<usize>,
    shifted_wavelengths:&[f64],
    flux_vec:&Vec<Vec<f64>>,
    cont_vec:&Vec<Vec<f64>>,
    wavelengths_vec:&Vec<Vec<f64>>,  
    temperature:f64,
    log_g:f64)->(Vec<f64>,Vec<f64>){
    
    let mut temperature_vec:Vec<f64>=Vec::new();
    let mut log_gravity_vec:Vec<f64>=Vec::new();
    for index in relevant_indices.iter(){
        temperature_vec.push(intensity_dfs.temperature_vector[*index]);
        log_gravity_vec.push(intensity_dfs.log_g_vector[*index]);
    }
    // Counter orders the intensities 
    let mut counter:usize = 0;
    
    let t0 = (temperature-temperature_vec[0])/(temperature_vec[2]-temperature_vec[0]);
    let t1 = (log_g-log_gravity_vec[0])/(log_gravity_vec[1]-log_gravity_vec[0]);
    
    //Initialize the collecting vectors.
    let mut interpolated_flux = vec![0.0;shifted_wavelengths.len()];
    let mut interpolated_cont= vec![0.0;shifted_wavelengths.len()];
    
    // loop over the four elements of the (temperature, log_gravity) parameter space
    for i in 0..=1{//<- i goes from 0 to 1
        let delta_temp = match i {
            0 => 1.0-t0, 
            1 => t0,
            _ =>{0.0}};//although verbose, computationally this is faster than multiplying floats

        for j in 0..=1{//<- j goes from 0 to 1
            let delta_logg = match j{
                0 => {1.0-t1}
                1 => {t1}
                _=> {0.0}};

            for (n,lambda) in shifted_wavelengths.iter().enumerate(){//<-this is to carry the same process onto every wavelength

                let index = search_geq(&wavelengths_vec[counter], *lambda);
                let t2 = (lambda - wavelengths_vec[counter][index-1])/
                (wavelengths_vec[counter][index]-wavelengths_vec[counter][index-1]);
                
                for k in 0..=1{
                    let delta_lambda = match k {
                        0 => {1.0-t2}
                        1 => {t2}
                        _ => {0.0} };
             
                    interpolated_flux[n] += 
                    get_contribution_from_vertex(delta_temp,
                        delta_logg,
                        delta_lambda,
                        flux_vec[counter][index-1]);
             
                    interpolated_cont[n] += 
                    get_contribution_from_vertex(delta_temp,
                        delta_logg,
                        delta_lambda,
                        cont_vec[counter][index-1]);
                }
            }
            counter += 1;//<- counter should from 0 to 3
        }
    }

    (interpolated_flux,interpolated_cont)    
}

/// The trilinear interpolation estimates a quantity that depends on position on a point inside of a cube where the
/// values are known on each of the vertices. This function returns the contribution given by the value on a particular vertex.
/// ### Arguments:
/// * `delta_temp` - This is the distance from the observed temperature value to the temperature value of a specific grid.
/// * `delta_logg` - This is the distance from the observed log gravity value to the log gravity value of a specific grid.
/// * `delta_lamb` - This is the distance from the observed wavelength value to the wavelength value inside a specific grid.
/// * `value` - This is the value on a specific vertex of the cube defined on the parameter space (temperature,log g, wavelength).
/// ### Returns:
/// * `f64` - This will hold the contribution given by a vertex. 
fn get_contribution_from_vertex(delta_temp:f64,
delta_logg:f64,
delta_lamb:f64,
value:f64)->f64{
    delta_temp*delta_logg*delta_lamb*value
}  