use polars::prelude::*;
use super::*;
use temp_name_lib::type_def::N_FLUX_POINTS;
mod define_parameter_space;
use std::{fmt::LowerExp, fs::File};

/// This module contains the methods and functions to process the intensity grids. 
/// It has two parts.
/// 
/// The first part consist of methods to produce a DataFrame that stores the file name of the relevant
/// intensity grid files
/// 
/// The second part consist of methods to obtain a LazyFrame from the contents of the intensity grid file and compute 
/// the Energy fluxes. 
/// 




/// The implementations for the profile config structure defined in this module are used to 
/// extract a DataFrame that works as a database for the intensity grid files. 
/// 
/// With this approach searching for the appropriate intensity grid files are carried out faster. 
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

    /// This function unwrapps the collection of the "grid file identification data", 
    /// producing
    /// * an ordered vector with all of the temperatures,
    /// * an ordered vector with all of the log g,
    /// * an ordered vector with all of the file names.
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
    /// This function creates a Data frame out of the intensity grid collection stored
    /// in the profile_config structure. 
    /// 
    /// The columns have the following headers
    /// 
    /// `temperature`| `log_gravity`| `file name`
    pub fn get_intensity_grids_dataframe(&self)->DataFrame{
        self.unwrap_grids().create_dataframe_sorted()
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

impl UnwrappedGrids{
    /// This function creates a data frame with the columns "temperature","log_gravity","file name"
    /// that are Float64, Float64, &[str] in that order. 
    /// 
    /// The data frame is ordered from lower to greater both in temperature and gravity. 
    fn create_dataframe_sorted(self)->DataFrame{
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
}


//----------------------------------------
//------parsing intensity grids-----------
//---------------------------------------- 

/// This function selects the four points closer to
/// (temperature,Log_g)
/// on the parameter space whose coordinates have an associated intensity grid file
/// 
/// ### Arguments:
/// * `grid_id_lf` - A polars LazyFrame that contains the information of the intensity grids.
/// * `temperature`` - a `f64` value for the temperature on a surface cell
/// * `log g`- a `f64` value for the logarithm of gravity on a surface cell
/// ### Returns:
/// * `(temperature_vec, log_g_vec, Filenames_vec)` - This is a tuple that contains three vectors where
/// * `temperature vec` - is a vector of `f64` that contains 4 ordered temperatures of the coordinates in the parameter space
/// * `log_g_vec` - is a vector of `f64 that contains 4 ordered values of log_g of the coordinates in the parameter space
/// * `Filenames_vec` - is a vector of `String` that contains the paths for the 4 intensity grid files.
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
pub fn get_flux_continuum(grid_filename:String,
    shifted_wavelengths:&[f64],
    coschi:f64,
    start_computation:&std::time::Instant)->Option< (Vec<f64> , Vec<f64>, Vec<f64>)>{
    if let Ok(lf) = read_intensity_grid_file(&grid_filename){

        // if the number of flux points to be calculated is really small, it's useful to filter out some of the wavelengths from the data frame so that the queries are performed faster.
        let lf_in_range = match filter_if_contains_wavelenghts(shifted_wavelengths){
            Some(relevant_range) => {lf.filter(relevant_range)}
            None => {lf}
        };
        let mut start = start_computation.elapsed();
        let mut duration = start_computation.elapsed() - start;
        println!("\n read intensity grid file and created lazyframe in {:?}", duration);
        start =start_computation.elapsed();
        let df_unfiltered = lf_in_range.collect().unwrap();
        duration = start_computation.elapsed()-start;
        println!(" created dataframe in {:?}", duration);
        
        start = start_computation.elapsed();
        let df_filtered = extract_relevant_wavelengths(shifted_wavelengths, df_unfiltered);
        duration = start_computation.elapsed()-start;
        println!("filtered dataframe in {:?}", duration);

        let lf_relevant = df_filtered.lazy();

        /*lf;match shifted_wavelengths.len() < (N_FLUX_POINTS/100) as usize{
            true => {
                lf
                .filter(filter_if_contains_wavelenght(&shifted_wavelengths,
                0.2).unwrap())}
            false => {lf}
        };*/

        // Here the lazy frame is transformed into a data frame that contains three columns: the intensity fluxes (specific, continuous) and the relevant wavelengths.
        start = start_computation.elapsed();        
        let df_flux = create_lf_wavelength_flux_continuum(lf_relevant,coschi).collect().unwrap();
        duration = start_computation.elapsed()-start;
        println!("final dataframe in {:?}", duration);

        start = start_computation.elapsed();
        let lbd_from_df = extract_column_as_vectorf64("wavelength", &df_flux);
        let flux_from_df = extract_column_as_vectorf64("flux", &df_flux);
        let continuum_from_df = extract_column_as_vectorf64("continuum", &df_flux);
        duration = start_computation.elapsed()-start;
        println!("gotten some vectors {:?}\n", duration);

        Some(( flux_from_df,
            continuum_from_df,
            lbd_from_df))
    }else{
    None
    }
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

/// This function constructs a polars expression [Expr] that filters out wavelengths that are not close to the 
/// observed (requested) ones. This reduces memory consumption and computation time.
/// 
/// ### Arguments: 
/// * `shifted_wavelengths` - The observed wavelengths
/// * `threshold` - a f64 value that specifies how close wavelengths need to be.
/// ### Returns:
/// *`Option<Expr>` -  where `Expr` is polars expression  that filters out unrelevant wavelengths from the lazyframe of an intensity grid file.
/// 
fn filter_if_contains_wavelenghts(shifted_wavelengths:&[f64])->Option<Expr>{

    let max_wavelength = shifted_wavelengths.last().unwrap().clone();
    let min_wavelength = shifted_wavelengths.get(0).unwrap().clone();

    let filter_lower_expr = col("wavelength").gt(lit(min_wavelength));
    let filter_greater_expr = col("wavelength").lt(lit(max_wavelength));

    let combined_filter_exp = filter_lower_expr.or(filter_greater_expr);
    Some(combined_filter_exp)
}

/// This is a function that extracts only the relevant wavelengths. 
/// 
/// ### Arguments: 
/// * `shifted_wavelengts` - A reference to a [Vec<f64>] that contains the doppler shifted observed wavelengths.
/// * `grid_dataframe` - A [DataFrame] that contains the wavelengths from the grid file. This DataFrame has exluded all the wavelengths outside of range. 
/// ### Returns: 
/// * `relevant_df` - A [DataFrame] that only contains the wavelengths that'll be used for the [interpolate] step.
fn extract_relevant_wavelengths(
    shifted_wavelengths:&[f64],
    grid_dataframe:DataFrame)->DataFrame{
    
    let df = &grid_dataframe;
    // Extract df's collumns into vectors; Wavelengths are ordered from lower to greater so we'll profit from that.
        let df_wavelengths=extract_column_as_vectorf64("wavelength", &grid_dataframe);
        let df_a = extract_column_as_vectorf64("a", &grid_dataframe);
        let df_b = extract_column_as_vectorf64("b", &grid_dataframe);
        let df_c = extract_column_as_vectorf64("c", &grid_dataframe);
        let df_d = extract_column_as_vectorf64("d", &grid_dataframe);
        let df_ac = extract_column_as_vectorf64("ac", &grid_dataframe);
        let df_bc = extract_column_as_vectorf64("bc", &grid_dataframe);
        let df_cc = extract_column_as_vectorf64("cc", &grid_dataframe);
        let df_dc = extract_column_as_vectorf64("dc", &grid_dataframe);
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
    let relevant_df= df![
        "wavelength"=>wavelengths,
        "a" => a,
        "b" => b,
        "c" => c,
        "d" => d,
        "ac" => ac,
        "bc" => bc,
        "cc" => cc,
        "dc" => dc
    ].unwrap();
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