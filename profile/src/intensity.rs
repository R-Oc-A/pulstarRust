use polars::prelude::*;
use super::*;
use temp_name_lib::type_def::N_FLUX_POINTS;
mod define_parameter_space;
use std::{fs::File, io::ErrorKind};

impl ProfileConfig{
    /// This function checks within the collection of intensity grid files provided by the toml file
    /// to see if all of the grid files are loaded into the directory.
    /// If they are not, it should return a ErrorKind::NotFound
    pub fn intensity_grids_are_loaded(&self)->Result<(),std::io::Error>{
        for grid in self.intensity_grids.iter(){
            grid.is_there_a_file(&self.path_to_grids)?;
        }
        Ok(())
    }

    /// This function unwrapps the collection of the grid file id. 
    /// Producing an ordered vector with all of the temperatures
    /// an ordered vector with all of the log g
    /// an ordered vector with all of the file names
    fn unwrap_grids(&self)->UnwrappedGrids{

        let mut temperature_vector:Vec<f64>=Vec::new();
        let mut logg_vector:Vec<f64>=Vec::new();
        let mut filename_vector:Vec<String>=Vec::new();
        for grid in self.intensity_grids.iter(){
            temperature_vector.push(grid.temperature);
            logg_vector.push(grid.log_gravity);
            filename_vector.push(
                format!("{}{}",
                self.path_to_grids,
                grid.filename.clone())
            );
        }
        UnwrappedGrids { temperatures: temperature_vector, 
            log_g: logg_vector, 
            filenames: filename_vector}
    }

}
impl IntensityGrid{
    /// This function asses whether the purported intensity grid file is stored in a given directory.
    /// ### Argument:
    /// * path - a string that indicates the relative path to the directory of the intensity grid files.
    /// ### Returns:
    /// * Ok(File) if it's able to locate the intensity file and is able to open it.
    /// * Error - otherwise.
    fn is_there_a_file(&self,path:&str)->Result<File,std::io::Error>{
        let file_name = self.filename.clone();
        let full_name = format!("{}{}",String::from(path), file_name);
        File::open(&full_name)
    }


}

/// This struct is useful to create a the database that handles all of the intensity grid files.
/// The data base produced using this structure will be easier to query from.
struct UnwrappedGrids{
    pub temperatures:Vec<f64>,
    pub log_g:Vec<f64>,
    pub filenames:Vec<String>,
}

impl IntensityGrids{
    //This function creates a data frame with the columns "temperature","log_gravity","file name"
    // that are Float64, Float64, &[str] in that order. The data frame is ordered from lower to greater both in temperature and gravity. 
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

        let temperatures = extract_column_as_vectorf64("temperature", &df);
        let log_g = extract_column_as_vectorf64("log_gravity", &df);
        let filenames = extract_column_as_vector_string("file name", &df);
        IntensityGrids { temperatures, log_g, filenames}
    }
}


//----------------------------------------
// Functions for parsing intensity grids
//---------------------------------------- 

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


// Returns the flux, continuum for each temperature,gravity as a f64 vector
// Input: name -> file name of the intensity grid file 
//        shifted_wavelengths -> collection of doppler shifted wavelenghts from where intensity shall be queried.
//        coschi -> projection of the surface cell normal onto the observer's unitary position vector;
pub fn get_flux_continuum(grid_filename:String,
    shifted_wavelengths:&[f64],
    coschi:f64)->Option< (Vec<f64> , Vec<f64>, Vec<f64>)>{
    if let Ok(lf) = read_intensity_grid_file(&grid_filename){
        
        let lf_relevant = match shifted_wavelengths.len() < (N_FLUX_POINTS/10) as usize{
            true => {
                lf
                .filter(filter_if_contains_wavelenght(&shifted_wavelengths,
                0.2).unwrap())}
            false => {lf}
        };

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


// This function returns a lazy frame with the apropiate column names of the intensity grid file
// the column names are 
// "wavelength","a","b","c","d","ac","bc","cc","dc"
// Input: path of the intensity grid file as a str
 fn read_intensity_grid_file(path:&str) -> PolarsResult<LazyFrame> {
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

// this function's purpose is to create a lazy frame containing
// lambda column (only with relevant wavelengths), Flux column created using the limb_darkening law and the continuum column.
// Input: lf is a lazy frame created with the read_intensity_grid_file function.
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



// This function applies a doppler shift to the observed wavelengths
pub fn get_doppler_shifted_wavelengths(doppler_shift:f64,
        wavelengths:&[f64])->Vec<f64>{
            let mut shifted:Vec<f64>=Vec::new();
            for lambda in wavelengths.iter(){
                shifted.push(*lambda * doppler_shift);
            }
            shifted
}

//this function filters
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

