use polars::prelude::*;
use temp_name_lib::type_def::CLIGHT;
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
    pub grids:Vec<GridFeatures>,//<--effective temperature,log gravity;name_of_grid_file
    pub lambda:Vec<f64>,//<-relevant wavelengths of grids

    pub t_eff_useful_index: usize,
    pub log_g_useful_index: usize,
    pub lambda_useful_index: usize,
}

pub struct GridFeatures{
    pub t_eff:f64,
    pub log_g:f64,
    pub file_name:String,
}

impl GridFeatures{
    pub fn new(t_eff_value:f64,log_g_value:f64,file_name:&str)->Self{
        GridFeatures { t_eff: t_eff_value,
             log_g: log_g_value,
             file_name: file_name.to_string() }
    }

}


pub fn read_intensity_grid_file(path:&str) -> PolarsResult<LazyFrame> {
    //let file = File::open(&path)?;

    let schema =  Schema::from_iter(vec![
		Field::new("wave_length".into(), DataType::Float64),
		Field::new("a".into(), DataType::Float64),
		Field::new("b".into(), DataType::Float64),
		Field::new("c".into(), DataType::Float64),
		Field::new("d".into(), DataType::Float64),
		Field::new("ac".into(), DataType::Float64),
		Field::new("bc".into(), DataType::Float64),
		Field::new("cc".into(), DataType::Float64),
		Field::new("dc".into(), DataType::Float64)
    ]); 
   
    /*let parse_options = CsvParseOptions{
        separator: b' ',
        ..Default::default()
    };*/

/*    let read_options = CsvReadOptions{
        has_header:false,
        schema: Some(Arc::new(schema)),
        parse_options: Arc::new(parse_options),
        ..Default::default()     
    };*/
    let lf= LazyCsvReader::new(path)
    .with_separator(b' ')
    .with_has_header(false)
    .with_schema(Some(Arc::new(schema))).finish()?;
    
/*    let df = CsvReader::new(file)
        .with_options (read_options)
        .finish()?;
*/
    Ok(lf)

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

    vec![vec![0.0]]
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