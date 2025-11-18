use super::*;

/// This function constructs a polars expression [Expr] that filters out wavelengths that are greater than or less than
/// observed (requested) ones. This reduces memory consumption and computation time.
/// 
/// ### Arguments: 
/// * `shifted_wavelengths` - The observed wavelengths
/// * `threshold` - a f64 value that specifies how close wavelengths need to be.
/// ### Returns:
/// *`Option<Expr>` -  where `Expr` is polars expression  that filters out unrelevant wavelengths from the lazyframe of an intensity grid file.
/// 
pub fn filter1_if_contains_wavelenghts(
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64)->Option<Expr>{

    
    let epsilon = 0.01;

    let min_wavelength = minval_rel_dopplershift 
        * wavelengths.iter()
        .fold(wavelengths[0],
            |accumulator,wavelength_val| 
            if*wavelength_val < accumulator {*wavelength_val} else {accumulator}
        )-epsilon;
    let max_wavelength = maxval_rel_dopplershift 
        * wavelengths.iter()
        .fold(wavelengths[0],
            |accumulator,wavelength_val| 
            if*wavelength_val > accumulator {*wavelength_val} else {accumulator}
        )+epsilon;    

    let filter_lower_expr = col("wavelength").gt(lit(min_wavelength));
    let filter_greater_expr = col("wavelength").lt(lit(max_wavelength));

    let combined_filter_exp = filter_lower_expr.or(filter_greater_expr);
    Some(combined_filter_exp)
}

// There might be necessary to apply another filter of the kind 
// Include wavelength only if observed_wavelength*min_rel_dopplershift - epsilon < wavelength <observed_wavelength * max_rel_dopplershift + epsilon
pub fn filter2_sift_wavelengths(
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64)->Option<Expr> {

    let epsilon = 1.0e-3;    
    let mut combined_expresion: Option<Expr> = None;
    
    for wavelength in wavelengths.iter(){
        let lb_wavelength = col("wavelength").gt(lit(wavelength*minval_rel_dopplershift - epsilon));
        let ub_wavelength = col("wavelength").lt(lit(wavelength*maxval_rel_dopplershift + epsilon));
        
        let current_mask = lb_wavelength.and(ub_wavelength);
        combined_expresion = match combined_expresion{
            Some(expression) => {Some(expression.or(current_mask))}
            None => {Some(current_mask)}
        }
    }    
    combined_expresion
}



///  This function materializes all of the filtered intensity data frames that will be used throughout the full program. It 
///  also fills an instance of the [IntensityDataFrames] 
/// 
/// ### Arguments:
/// * `grids_db` - a [DataFrame] that contains the name of the relevant gridfiles and their associated temperature and log_g values.
/// * `wavelengths` - a &[f64] vector (slice) that contains the observed wavelengths
/// * `max_rel_dopplershift` - a [f64] value that contains the maximum relative dopplershift
/// * `min_rel_dopplershift` - a [f64] value that contians the minimum relative dopplershift
/// 
/// ### Returns:
/// * an Instance of [IntensityDataFrames] where the temperature, log_g are ordered and it contains the [DataFrame]s of the intensity grid files
/// 
pub fn parse_relevant_intensity_grids(
    grids_db: DataFrame,
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64,
)->IntensityDataFrames{
    
    let temperature_vector = extract_column_as_vectorf64("temperature", &grids_db);
    let log_g_vector = extract_column_as_vectorf64("log_gravity", &grids_db);
    let filenames = extract_column_as_vector_string("file name", &grids_db);

    let is_low_resolution = (wavelengths[1]-wavelengths[0])>5.0e-3;//if the requested wavelengths are separated by more than 0.005 nm
    let mut intensity_dfs:Vec<DataFrame> = Vec::new();
    for name in filenames{
        let lf = read_intensity_grid_file(&name).expect("Unable to read intensity grid file");
        println!("loading {} grid file",name);

        let temp_df =lf.filter(
            filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift).unwrap()
            ).collect().expect("unable to produce dataframe using the intensity grid files");
        if is_low_resolution {
            intensity_dfs.push(
                sifted_wavelengths_dataframe(wavelengths,
                    temp_df.lazy(),
                    maxval_rel_dopplershift,
                    minval_rel_dopplershift).expect("unable to produce dataframe using the intensity grid files")
                );
        }
        else {
            intensity_dfs.push(temp_df);
        }
    }
    IntensityDataFrames{
        temperature_vector:temperature_vector,
        log_g_vector: log_g_vector,
        intensity_dfs:intensity_dfs,
    }
}


/// This function returns an intensity dataframe with only relevant wavelengths
/// The main purpose of this function is when you have a  low_resolution wavelength spectra
/// with bigger number of waves than 500. 
/// 
/// [polars] has a difficult time dealing with the combined expresion of more than 500 conditionals so the strategy implemented here
/// is to sepparate an intensity dataframe into chunks and sift the relevant wavelengths by chunks. 
/// 
/// ### Arguments:
/// * `wavelengths` - a [Vec<f64>] collection that contains the requested wavelengths.
/// * `intensity_df` - An intensity [DataFrame] about to be filtered
/// ### Returns: 
/// This function returns a [PolarsResult] with the following variants
/// * `Ok(DataFrame)` -where the binded data frame contains the sifted wavelengths
fn sifted_wavelengths_dataframe(
    wavelengths:&[f64],
    temp_lf: LazyFrame,
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64
)->PolarsResult<DataFrame>{
    //separate wavelengths into chunks of 500
    let chunk_size = 500usize;

    //initiallize collecting dataframe
    let mut collecting_df = df!("empty"=>[0.0]).unwrap();

    //loop over chunks
    for (iteration,chunk_of_wavelengths) in wavelengths.chunks(chunk_size).enumerate(){

        let partially_sifted_lf = temp_lf.clone().filter(filter2_sift_wavelengths(
            chunk_of_wavelengths,
            maxval_rel_dopplershift,
            minval_rel_dopplershift).unwrap());

        let new_df = 
        if iteration == 0{
            partially_sifted_lf.collect()?
        }
        else{
            let old_lf = collecting_df.lazy();
            concat([old_lf,partially_sifted_lf],
            UnionArgs::default())?.collect()?
        };

        collecting_df = new_df;        
    }
    Ok(collecting_df)
    //initialize the collecting dataframe
    //There's an iterator function called .chunks(usize) I shall look into it
    //loop over chunks
        //create data frame with the filetered wavelenghts from chunk
        //append with collecting data frame
            //create a lazyframe with thefiltered wavelengths dataframe 
            // create a lazyframe witht the collecting dataframe
    //
    //
}
