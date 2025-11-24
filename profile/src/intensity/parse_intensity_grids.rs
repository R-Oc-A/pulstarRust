use super::*;

pub mod joris_grids;

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

    let filter_lower_expr = col("wavelengths").gt(lit(min_wavelength));
    let filter_greater_expr = col("wavelengths").lt(lit(max_wavelength));

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
        let lb_wavelength = col("wavelengths").gt(lit(wavelength*minval_rel_dopplershift - epsilon));
        let ub_wavelength = col("wavelengths").lt(lit(wavelength*maxval_rel_dopplershift + epsilon));
        
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
pub fn filter_wavelength_range(
    grids_lf: LazyFrame,
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64,
)->LazyFrame{
    
    //Nadya's grids are in Angstroms while Joris's are in nm. To check if the requested wavelengths are low resolution i.e. dλ,1e-3nm, I need to specify that or use the same wavelenght units.
    let is_low_resolution=false;
    let combined_expresion=match is_low_resolution{
        true => {filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift).or(
            filter2_sift_wavelengths(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift)
        )}
        false => {filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift)}
    };

    match combined_expresion{
        Some(expresion)=>{grids_lf.filter(expresion)}
        None=>{panic!("unable to produce dataframe using the intensity grid files")}
        
    }
}


//add interpolating profile test for each fractional coordinate.
