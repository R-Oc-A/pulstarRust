use polars::prelude::*;

pub fn get_rectangles_lazy(lf:LazyFrame,
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

    let combined_filter_mask = col("temperature").lt_eq(lit(temperature))
        .and(col("log_gravity").lt_eq(lit(log_gravity)));
    let lower_temp_lower_grav_lf = lf.clone()
    .filter(combined_filter_mask)
    .sort(["temperature","log_gravity"], 
    SortMultipleOptions::new()
    .with_order_descending_multi([true,true]))
    .first();

    let combined_filter_mask = col("temperature").lt_eq(lit(temperature))
        .and(col("log_gravity").gt_eq(lit(log_gravity)));
    let lower_temp_upper_grav_lf = lf.clone()
    .filter(combined_filter_mask)
    .sort(["temperature","log_gravity"], 
    SortMultipleOptions::new()
    .with_order_descending_multi([true,false]))
    .first();
    
    let combined_filter_mask = col("temperature").gt_eq(lit(temperature))
        .and(col("log_gravity").lt_eq(lit(log_gravity)));
    let upper_temp_lower_grav_lf = lf.clone()
    .filter(combined_filter_mask)
    .sort(["temperature","log_gravity"], 
    SortMultipleOptions::new()
    .with_order_descending_multi([false,true]))
    .first();

    let combined_filter_mask = col("temperature").gt_eq(lit(temperature))
        .and(col("log_gravity").gt_eq(lit(log_gravity)));
    let upper_temp_upper_grav_lf = lf.clone()
    .filter(combined_filter_mask)
    .sort(["temperature","log_gravity"], 
    SortMultipleOptions::new()
    .with_order_descending_multi([false,false]))
    .first();

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
    .with_order_descending([is_lower_temp_bnd,is_lower_grav_bnd]))
    .first();
}