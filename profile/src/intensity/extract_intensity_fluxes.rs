use super::*;
use parse_intensity_grids::{read_intensity_grid_file,filter1_if_contains_wavelenghts, extract_relevant_wavelengths};

/// This function extracts the intensity
/// 
/// ### Arguments: 
/// * `name` - file name of the intensity grid file 
/// * `shifted_wavelengths` - collection of doppler shifted wavelenghts from where intensity shall be queried.
/// * `coschi` - projection of the surface cell normal onto the observer's unit position vector;
/// ### Returns: 
/// `(Flux, Continuum, wavelength_vec)` - a tupple containing three vectors of type `f64` where
/// `Flux` - Is the intensity flux
/// `Continuum` - continuus intensity flux
/// `wavelength_vec` - is the wavelength values where the intensities are computed.
pub fn get_flux_continuum(grid_filename:String,
    shifted_wavelengths:&[f64],
    coschi:f64,
    start_computation:&std::time::Instant)->Option< (Vec<f64> , Vec<f64>, Vec<f64>)>{
    if let Ok(lf) = read_intensity_grid_file(&grid_filename){

        // if the number of flux points to be calculated is really small, it's useful to filter out some of the wavelengths from the data frame so that the queries are performed faster.
        let lf_in_range = match filter1_if_contains_wavelenghts(shifted_wavelengths){
            Some(relevant_range) => {lf.filter(relevant_range)}
            None => {lf}
        };
        let df_unfiltered = lf_in_range.collect().unwrap();
        let df_filtered = extract_relevant_wavelengths(shifted_wavelengths, df_unfiltered);
        let lf_relevant = df_filtered.lazy();

        /*lf;match shifted_wavelengths.len() < (N_FLUX_POINTS/100) as usize{
            true => {
                lf
                .filter(filter1_if_contains_wavelenght(&shifted_wavelengths,
                0.2).unwrap())}
            false => {lf}
        };*/

        // Here the lazy frame is transformed into a data frame that contains three columns: the intensity fluxes (specific, continuous) and the relevant wavelengths.
        let df_flux = create_lf_wavelength_flux_continuum(lf_relevant,coschi).collect().unwrap();

        let lbd_from_df = extract_column_as_vectorf64("wavelength", &df_flux);
        let flux_from_df = extract_column_as_vectorf64("flux", &df_flux);
        let continuum_from_df = extract_column_as_vectorf64("continuum", &df_flux);

        Some(( flux_from_df,
            continuum_from_df,
            lbd_from_df))
    }else{
    None
    }
}

/// This function receives a lazy frame containing the coefficients of the limb
/// darkening law
/// $I_{\lambda}(\mu) =a_0 \sum_{k=1}^{3} a_k \left( 1 - \mu^k  \right)$
/// and returns another lazy frame with the wavelengths and their intensity fluxes.
///
/// ### Arguments:
/// * `lf` - A lazy frame read from the intensity grid file
/// * `coschi` - a `f64` value that is the projection of the unit surface normal onto the unit vector towards the observer.
/// ### Returns:
/// * `LazyFrame` - The column headers of this lazy frame are given by 
/// * `|wavelength|flux|continuum|`
 fn create_lf_wavelength_flux_continuum(lf:LazyFrame, coschi:f64)->LazyFrame{
    let mu = coschi.sqrt();

    let bcoef = 1.0 - mu; 
    let ccoef = 1.0 - coschi;
    let dcoef = 1.0-mu.powi(3);

    let i_lambda = col("a") 
        + col("b") * lit(bcoef) 
        + col("c") * lit(ccoef)
        + col("d") * lit(dcoef);

    let i_cont = col("ac") 
        + col("bc") * lit(bcoef) 
        + col("cc") * lit(ccoef)
        + col("dc") * lit(dcoef);
    
    lf.select([
        col("wavelength"),i_lambda.alias("flux"),
        i_cont.alias("continuum")
    ])
}
/// This function selects the four points closer to
/// (temperature,Log_g)
/// on the parameter space whose coordinates have an associated intensity grid file
/// 
/// ### Arguments:
/// * `grid_id_lf` - A polars LazyFrame that contains the information of the intensity grids.
/// * `temperature`` - a `f64` value for the temperature on a surface cell
/// * `log g`- a `f64` value for the logarithm of gravity on a surface cell
/// ### Returns:
/// * `(temperature_vec, log_g_vec, Filenames_vec)` - This is a tuple that contains three vectors where
/// * `temperature vec` - is a vector of `f64` that contains 4 ordered temperatures of the coordinates in the parameter space
/// * `log_g_vec` - is a vector of `f64 that contains 4 ordered values of log_g of the coordinates in the parameter space
/// * `Filenames_vec` - is a vector of `String` that contains the paths for the 4 intensity grid files.
pub fn get_temp_logg_filenames(grid_id_lf:LazyFrame,
    temperature:f64,
    log_gravity:f64)->(Vec<f64>,Vec<f64>,Vec<String>){
    let lf_relevant = define_parameter_space::get_rectangles_in_parameter_space(grid_id_lf, temperature, log_gravity);
    let df_relevant = lf_relevant.collect().unwrap();

    let temperature_vec = extract_column_as_vectorf64("temperature", &df_relevant);
    let log_g_vec = extract_column_as_vectorf64("log_gravity", &df_relevant);
    let name_of_files = extract_column_as_vector_string("file name", &df_relevant);

    (temperature_vec,log_g_vec , name_of_files)
}
// 1) Get the maximum doppler shift value, temperature, and log_g
// 2) Open and filter all of the relevant grid files
// 3) rewrite the methods herebelow so that queries and column operations are carried on lazyframes related to this df
// 4) do some benchmark
// 5) If nothing works, have a meeting with Joris about How to improve speed. Might be necessary to look into parallelization. 

