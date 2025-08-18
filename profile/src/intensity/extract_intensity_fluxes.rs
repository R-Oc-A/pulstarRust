use super::*;

/// This is a function that extracts only the relevant wavelengths as [Vec<f64>]. 
/// 
/// ### Arguments: 
/// * `shifted_wavelengts` - A reference to a [Vec<f64>] that contains the doppler shifted observed wavelengths.
/// * `grid_dataframe` - A [DataFrame] that contains the wavelengths from the grid file. This DataFrame has exluded all the wavelengths outside of range. 
/// ### Returns: 
/// * `relevant_df` - A [DataFrame] that only contains the wavelengths that'll be used for the [interpolate] step.
pub fn extract_relevant_wavelengths(
    shifted_wavelengths:&[f64],
    grid_dataframe:&DataFrame)->DataFrame{
    
    // Extract df's collumns into vectors; Wavelengths are ordered from lower to greater so we'll profit from that.
        let df_wavelengths=extract_column_as_vectorf64("wavelength", grid_dataframe);
        let grid_precision = (df_wavelengths[1]-df_wavelengths[0]).abs()+5.0e-3;
        let wavelength_precision = (shifted_wavelengths[1]-shifted_wavelengths[0]).abs();
    let relevant_df = match wavelength_precision>grid_precision{// If the requested wavelength spacing is bellow the grid precision, then all of the grids wavelengths are relevant ones.
        true => {
            let df_a = extract_column_as_vectorf64("a", grid_dataframe);
	        let df_b = extract_column_as_vectorf64("b", grid_dataframe);
	        let df_c = extract_column_as_vectorf64("c", grid_dataframe);
	        let df_d = extract_column_as_vectorf64("d", grid_dataframe);
	        let df_ac = extract_column_as_vectorf64("ac", grid_dataframe);
	        let df_bc = extract_column_as_vectorf64("bc", grid_dataframe);
	        let df_cc = extract_column_as_vectorf64("cc", grid_dataframe);
	        let df_dc = extract_column_as_vectorf64("dc", grid_dataframe);
		    // Get relevant indices of the df_wavelength vector into an index vector
	        let indices = extract_important_indices(shifted_wavelengths, &df_wavelengths);
    	    // Construct new vectors using only the relevant indices
	        let wavelengths = construct_new_vec_using_indices(df_wavelengths, &indices);
	        let a = construct_new_vec_using_indices(df_a, &indices);
	        let b = construct_new_vec_using_indices(df_b, &indices);
	        let c = construct_new_vec_using_indices(df_c, &indices);
	        let d = construct_new_vec_using_indices(df_d, &indices);
	        let ac = construct_new_vec_using_indices(df_ac, &indices);
	        let bc = construct_new_vec_using_indices(df_bc, &indices);
	        let cc = construct_new_vec_using_indices(df_cc, &indices);
	        let dc = construct_new_vec_using_indices(df_dc, &indices);
	
	        // Construct a data frame using the relevant indices
		    df![
		        "wavelength"=>wavelengths,
		        "a" => a,
		        "b" => b,
		        "c" => c,
		        "d" => d,
		        "ac" => ac,
		        "bc" => bc,
		        "cc" => cc,
		        "dc" => dc
		    ].unwrap()
	        }
	    false => {
	            grid_dataframe.clone()
	        }
	    };
	    relevant_df
}

/// This function is useful to get the relevant indices out of the wavelength vector. 
///
/// ### Arguments:
/// * `shifted_wavelengts` - A reference to a [Vec<f64>] that contains the doppler shifted observed wavelengths.
/// * `df_wavelengts` - A reference to a [Vec<f64>] that contains the grid's wavelengths extracted from the grid [DataFrame].
/// ### Returns:
/// * `indices` - A [Vec<usize>] that contains the indices of the df_wavelengths that'll be used for the [interpolate] step.
fn extract_important_indices(
    shifted_wavelengths:&[f64],
    df_wavelengths: &[f64]
)->Vec<usize>{
    let mut indices:Vec<usize> = Vec::new();
    for wavelength in shifted_wavelengths.iter(){
        let index = search_geq(df_wavelengths, *wavelength);
        //extract last index
        /*if let Some(last_index)= indices.last(){
            if df_wavelengths[*last_index]!=df_wavelengths[index-1]{
                indices.push(index-1)
            }
        }else{
            indices.push(index-1);
        }*/
        indices.push(index-1);
        indices.push(index);
    }
    indices
}

/// This function filters a vector using a collection of indices
/// 
/// ### Arguments: 
/// * `original_vector` - A [Vec<f64>] that contains the unfiltered data from the grid files
/// * `index_collection` - A [ &[usize] ] that contains the relevant indices.
/// ### Returns:
/// * `filtered_vector` - A [Vec<f64>] that contains the filtered data from the grid files. 
fn construct_new_vec_using_indices(
    original_vector:Vec<f64>,
    index_collection:&[usize]
)->Vec<f64>{
    let mut filtered_vector:Vec<f64>= Vec::new();
    for index in index_collection.iter(){
        filtered_vector.push(original_vector[*index]);
    }
    filtered_vector
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
    let dcoef = 1.0 - mu.powi(3);

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
pub fn get_flux_continuum(
    grid_dataframe: &DataFrame,
    shifted_wavelengths:&[f64],
    coschi:f64
    )->Option< (Vec<f64> , Vec<f64>, Vec<f64>)>{

    let df = extract_relevant_wavelengths(shifted_wavelengths, grid_dataframe);


    // Here the lazy frame is transformed into a data frame that contains three columns: the intensity fluxes (specific, continuous) and the relevant wavelengths.
    let df_flux = create_lf_wavelength_flux_continuum(df.lazy(),coschi).collect().unwrap();

    let lbd_from_df = extract_column_as_vectorf64("wavelength", &df_flux);
    let flux_from_df = extract_column_as_vectorf64("flux", &df_flux);
    let continuum_from_df = extract_column_as_vectorf64("continuum", &df_flux);

    Some(( flux_from_df,
        continuum_from_df,
        lbd_from_df))
}

