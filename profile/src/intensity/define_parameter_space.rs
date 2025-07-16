use polars::prelude::*;

/// This function helps find out the 4 (temperature,log_gravity) grid files
/// that can approximate the log_g and temperature on a cell of the rasterized star.
/// 
/// ### Arguments: 
/// 
/// `lf` - Is a lazy frame that contains a database of the intensity grid files. Its column headers are `|temperature|log_gravity|file name|`
/// 
/// `temperature` - a `f64` value that holds the temperature over a specific surface cell. 
/// 
/// `log_gravity` -  a `f64` value that holds the log gravity over a specific surface cell.
/// 
/// ### Returns:
/// 
/// `LazyFrame` - a LazyFrame that has  only the four relevant rows for a specific surface cell to perform interpolation. 
///  
pub fn get_rectangles_in_parameter_space(lf:LazyFrame,
    temperature:f64,
    log_gravity:f64)->LazyFrame{

    //find the intensity grid file that has the closest lower temperature and lower log g to the requested values
    let lower_temp_lower_grav_lf = 
    get_row(true,true,
            &lf, temperature, log_gravity);
    
    //find the intensity grid file that has the closest greater temperature and lower log g to the requested values
    let upper_temp_lower_grav_lf = 
    get_row(false,true,
            &lf, temperature, log_gravity);

    //find the intensity grid file that has the closest lower temperature and greater log g to the requested values
    let lower_temp_upper_grav_lf = 
    get_row(true,false,
            &lf, temperature, log_gravity);

    //find the intensity grid file that has the closest greater temperature and greater log g to the requested values
    let upper_temp_upper_grav_lf = 
    get_row(false,false,
            &lf, temperature, log_gravity);

    concat([lower_temp_lower_grav_lf,
        lower_temp_upper_grav_lf,
        upper_temp_lower_grav_lf,
        upper_temp_upper_grav_lf],UnionArgs::default()).unwrap()
}


/// This function looks for the grid file that has the closest temperature and log g to requested values.
/// 
/// ### Arguments:
/// 
/// `is_lower_temp_bnd` - a boolean that indicates whether it's searching for the closest lower of greater temperature.
/// 
/// `is_lower_grav_bnd` - a boolean that indicates whether it's searching for the closest lower of greater gravity.
/// 
/// `lf` - a polars LazyFrame that functions as a database for the intensity grid files.
/// 
/// `temperature` - a `f64` value that holds the reference temperature (in Kelvin)
/// 
/// `log_gravity` - a `f64` value that holds the log gravity reference.
fn get_row(is_lower_temp_bnd:bool,
    is_lower_grav_bnd:bool,
    lf:&LazyFrame,
    temperature:f64,
    log_gravity:f64)->LazyFrame{
    let mut combined_filter_mask = col("temperature");

    match is_lower_temp_bnd{
        true =>{ combined_filter_mask = combined_filter_mask.lt_eq(lit(temperature))}
        false => {combined_filter_mask = combined_filter_mask.gt_eq(lit(temperature))}
    };
    combined_filter_mask = match is_lower_grav_bnd{
        true => {combined_filter_mask.and(col("log_gravity").lt_eq(lit(log_gravity)))}
        false => {combined_filter_mask.and(col("log_gravity").gt_eq(lit(log_gravity)))}
    };

    lf.clone()
    //get the appropriate "quadrant" in the parameter space
    .filter(combined_filter_mask)
    //sort the obtained quadrant
    .sort( ["temperature","log_gravity"],
    SortMultipleOptions::new()
    .with_order_descending_multi([is_lower_temp_bnd,is_lower_grav_bnd]))
    //get the first element of the quadrant
    .first()
}