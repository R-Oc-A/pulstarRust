use polars::prelude::*;
use temp_name_lib::type_def::{CLIGHT, N_FLUX_POINTS};
use std::fs::File;
use std::sync::Arc;
pub struct config{
    pub lambda_0:f64,
    pub lambda_f:f64,
    pub delta_lbd:f64,
    pub v_max: f64,
    pub n_phases:u32,
}


pub struct IntensityGrids{
    temperatures:Vec<f64>,
    log_g:Vec<f64>,
    filenames:Vec<String>,
}

impl IntensityGrids{
    pub fn create_dataframe_sorted(self)->DataFrame{
        let filenames_str:Vec<&str> = self.filenames.iter().map(|s| s.as_str()).collect();

        let temperature_series = Series::new("temperature".into(),self.temperatures);
        let log_gravity_serie=Series::new("log_gravity".into(),self.log_g);
        let file_name_serie = Series::new("file name".into(),filenames_str);
        
        let df = DataFrame::new(vec![temperature_series.into(),
            log_gravity_serie.into(),
            file_name_serie.into()]).unwrap()
            .sort(["temperature","log_gravity"], 
            SortMultipleOptions::default()
            .with_order_descending_multi([false,false]));

       df.unwrap()
    }

    pub fn sort_intensitygrid(self)->IntensityGrids{
        let df = self.create_dataframe_sorted();
        let temperature_series = df.column("temperature").unwrap();
        let log_gravity_series = df.column("log_gravity").unwrap();
        let file_name_series = df.column("file name").unwrap();

        let temperatures:Vec<f64> = temperature_series.f64().unwrap().into_iter().flatten().collect();
        let log_g:Vec<f64> = log_gravity_series.f64().unwrap().into_iter().flatten().collect();
        let filenames_str:Vec<&str> = file_name_series.str().unwrap().into_iter().flatten().collect();

        let filenames:Vec<String> = filenames_str.iter().map(|s| s.to_string()).collect();

        IntensityGrids { temperatures, log_g, filenames}
    }
}

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

pub fn get_rectangles_lazy(lf:LazyFrame,
    temperature:f64,
    log_gravity:f64)->LazyFrame{

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

//function that appends doppler shift
pub fn append_doppler_shift(df:LazyFrame)// I take ownership of the data frame since I will produce a new one and want the old one to be dropped after appending
->LazyFrame{
    let doppler_shift = lit(1.0) - col("velocity")/lit(CLIGHT)*lit(1.0e3);
    let result =
    df.clone().lazy().
    with_column(doppler_shift.clone().alias("doppler shift"));
    result
}

//this function receives a lazy frame read from the .txt intensity grid files. It's purpose is to create a lf containing
// lambda column (only with relevant wavelengths), Flux_column created using the limb_darkening law and the continuum column.
pub fn create_lf_wavelength_flux_continuum(lf:LazyFrame, coschi:f64)->LazyFrame{
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
        col("wavelength"),i_lambda.alias("flux"),i_cont.alias("continuum")
    ])
}
//Estoy Aqui.


//get flux, continuum for each temperature gravity as an array
pub fn get_flux_continuum(profile_config:config,
    grid_filename:String,
    wavelengths:&[f64],
    doppler_shift:f64,
    coschi:f64)->Option< Vec<Vec<f64>> >{
    if let Ok(lf) = read_intensity_grid_file(&grid_filename){
        let shifted_wavelengths = get_doppler_shifted_wavelengths(doppler_shift,
             wavelengths);
        
        let lf_relevant = match wavelengths.len() < (N_FLUX_POINTS/10) as usize{
            true => {
                lf
                .filter(filter_if_contains_wavelenght(&shifted_wavelengths,
                0.2).unwrap())}
            false => {lf}
        };

        let lf_flux = create_lf_wavelength_flux_continuum(lf_relevant,coschi).collect().unwrap();
        
        let lbd_from_df = extract_column_as_vector("wavelength", &lf_flux);
        let flux_from_df = extract_column_as_vector("flux", &lf_flux);
        let continuum_from_df = extract_column_as_vector("continuum", &lf_flux);

        Some(vec![lbd_from_df,
            flux_from_df,
            continuum_from_df])
    }else{
    None
    }
}
//----------------------------------------
//EStoy aqui
fn extract_column_as_vector(column_name: &str,df:&DataFrame)->Vec<f64>{
        let column = df.column(column_name).unwrap();
        let vec_result:Vec<f64> = column.f64().unwrap().into_iter().flatten().collect();
        vec_result
}

//function that returns the flux for each cell
pub fn return_flux_for_cell_thetaphi(
    coschi:f64,
    temperature:f64,
    log_gravity:f64,
    doppler_shift:f64,
    area:f64,
    wavelengths:&[f64],
    )->Vec<Vec<f64>>{
	//function that selects 4 grids temperature +_ times gravity+_ ..and 
    //a function that creates lazyframes for the intensity grids
    //filter values that are outside of wavelength range
    //construct a lambda array dopplershifted
    //look for a function that select lambdas that are close to the values of the lambda array
    let shifted_wavelengths = get_doppler_shifted_wavelengths(doppler_shift, wavelengths);
    let filter_exp = filter_if_contains_wavelenght(&shifted_wavelengths, 
        0.1); 
    
    //function that appends Ic and I to each data frame
	//function that linearly interpolates Ic and I from the 8 grids
    //

    vec![vec![0.0],vec![0.0]]
}

fn get_doppler_shifted_wavelengths(doppler_shift:f64,
        wavelengths:&[f64])->Vec<f64>{
            let mut shifted:Vec<f64>=Vec::new();
            for lambda in wavelengths.iter(){
                shifted.push(*lambda);
            }
            shifted
}

fn filter_if_contains_wavelenght(shifted_wavelengths:&[f64],threshold:f64)->Option<Expr>{

    let mut combined_filter_exp: Option<Expr> = None;

    for wavelength in shifted_wavelengths.iter(){
        let lower_bound = lit(*wavelength-threshold);
        let upper_bound = lit(*wavelength + threshold);
        
        let current_filter_exp = col("wave_length")
        .gt_eq(lower_bound)
        .and(col("wave_length").lt_eq(upper_bound));

        combined_filter_exp = match combined_filter_exp{
            Some(expr) => {Some(expr.or(current_filter_exp))}
            None => {Some(current_filter_exp)}
        }
    }   
    combined_filter_exp
}


/// Computes by linear interpolation the toal and continuum intensity
/// Input: 
/// Output: (intens_result,cont_result) where intens is the total intensity; cont is the continuum intensity
/// 
pub fn interpolate(
    intensity: IntensityGrids,
    t_eff_value:f64,
    log_g_value:f64,
    cos_chi:f64,
    lambda:f64,
    )->[f64;2]{

        //let 
    //let grid_index = search_geq(&intensity.grids.t_eff, t_eff_value);

    //opens relevant grid files only
    //filters useless lambda values
    

    [0.0,0.0]
}



fn search_geq(vector:&[f64],key:f64)-> usize {
    //non empty vector
    let first_element =vector.get(0);
    if let Some(first_value) = first_element{
        if key <= *first_value {0usize}
        else{
            let size = vector.len();
            let mut top = size - 1;
            let mut bottom = 0;
            let mut middle = bottom + (top - bottom) / 2;
            if key > vector[top]{panic!("key value is too large")}

            while bottom<top{
                if vector[middle] <key{
                    bottom = middle + 1;
                } else {
                    top = middle;
                }
                middle = bottom + (top - bottom)/2;
            }
            top
        }
    }
    else{ panic!("empty vector!")}
}