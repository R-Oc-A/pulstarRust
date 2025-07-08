use polars::prelude::*;

// This function helps find out the 4 (temperature,log_gravity) grid files
// that can approximate the log_g and temperature on a cell of the raster sphere.
// Input: Is a lazy frame that contains as columns all names of the 
// intensity grids, the temperature associated with the grid, and the log gravity associated with the grid.
// Output: A lazy frame that filters only the four relevant rows.
pub fn get_rectangles_in_parameter_space(lf:LazyFrame,
    temperature:f64,
    log_gravity:f64)->LazyFrame{

    let lower_temp_lower_grav_lf = 
    get_row(true,true,
            &lf, temperature, log_gravity);
    
    let upper_temp_lower_grav_lf = 
    get_row(false,true,
            &lf, temperature, log_gravity);

    let lower_temp_upper_grav_lf = 
    get_row(true,false,
            &lf, temperature, log_gravity);

    let upper_temp_upper_grav_lf = 
    get_row(false,false,
            &lf, temperature, log_gravity);

    concat([lower_temp_lower_grav_lf,
        lower_temp_upper_grav_lf,
        upper_temp_lower_grav_lf,
        upper_temp_upper_grav_lf],UnionArgs::default()).unwrap()
}

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
    .filter(combined_filter_mask)
    .sort( ["temperature","log_gravity"],
    SortMultipleOptions::new()
    .with_order_descending_multi([is_lower_temp_bnd,is_lower_grav_bnd]))
    .first()
}