use super::*;

/// This structure has a collection of the [DataFrame]s extracted from the intensity grid file. This [DataFrame]s have been filtered
/// to include only the usefull wavelengths
/// 
pub struct IntensityDataFrames{
    /// Collection of the associated temperature values for the intensity dataframe
    temperature_vector: Vec<f64>,
    /// Collection of the associated log gravity values for the intensity dataframe
    log_g_vector: Vec<f64>,
    /// Collection of the dataframes extracted from the grid files.
    intensity_dfs: Vec<DataFrame>,
}

/// This function is used to create a polars LazyFrame out of the intensity grid file in order to perform 
/// column wise operations faster. 
/// The Kurukz intensity files contain the coefficients $a_k$ useful to calculate the intensity via the limb darkening law
/// given by 
/// $ I_{\lambda}(\mu) =a_0 \sum_{k=1}^{3} a_k \left( 1 - \mu^k  \right) $
/// where $\mu$ is the $cos(\xi)$
///  
/// ### Arguments:
/// * `path` - a reference to a string slice that contains the relative path to the intensity grid file
/// ### Returns:
/// * `PolarsResult<LazyFrame>` - where the lazy frame has as headers
/// * `|wavelength|a|b|c|d|ac|bc|cc|dc|`
// This function returns a lazy frame with the apropiate column names of the intensity grid file
// the column names are 
// "wavelength","a","b","c","d","ac","bc","cc","dc"
// Input: path of the intensity grid file as a str
 pub fn read_intensity_grid_file(path:&str) -> PolarsResult<LazyFrame> {
    //let file = File::open(&path)?;

    let schema =  Schema::from_iter(vec![
		Field::new("wavelength".into(), DataType::Float64),
		Field::new("a".into(), DataType::Float64),
		Field::new("b".into(), DataType::Float64),
		Field::new("c".into(), DataType::Float64),
		Field::new("d".into(), DataType::Float64),
		Field::new("ac".into(), DataType::Float64),
		Field::new("bc".into(), DataType::Float64),
		Field::new("cc".into(), DataType::Float64),
		Field::new("dc".into(), DataType::Float64)
    ]); 
   
    let lf= LazyCsvReader::new(path)
    .with_separator(b' ')
    .with_has_header(false)
    .with_schema(Some(Arc::new(schema))).finish()?;
    
    Ok(lf)

}

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

    let epsilon = 0.01;    
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
    let filenames = extract_column_as_vector_string("file names", &grids_db);

    let mut intensity_dfs:Vec<DataFrame> = Vec::new();
    for name in filenames{
        let lf = read_intensity_grid_file(&name).expect("Unable to read intensity grid file");
        intensity_dfs.push(lf.filter(
            filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift).unwrap().and(
            filter2_sift_wavelengths(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift).unwrap()
            )).collect().expect("unable to produce dataframe using the intensity grid files"));
    }
    IntensityDataFrames{
        temperature_vector:temperature_vector,
        log_g_vector: log_g_vector,
        intensity_dfs:intensity_dfs,
    }
}
