use polars::{error::ErrString, prelude::*};
use serde::Deserialize;
use temp_name_lib::type_def::{CLIGHT,N_FLUX_POINTS};//Velocity of light in m/s
use ndarray;

use std::fs;

/// The profile program reads a the quantities calculated over a Rasterized star on selected time points
/// and reurns a parquet file that contains the time series of the mean intensity flux produced 
/// for relevant observed wavelengths.


// I'll translate everything to Nadyalike grids. It's easier to work with only one type of spectral grids. 

mod intensity;
pub mod utils;

pub mod regresor_template;

/// This structure holds the data to construct the synthetic normalized flux.
#[derive(Clone)]
pub struct FluxOfSpectra{
    /// [Vec<f64>] containing the current phase of the pulsation cycle.
    pub time: Vec<f64>,
    /// [Vec<f64>] containing the current observed wavelengths.  
    pub wavelengths: Vec<f64>,
    /// [Vec<f64>] the doppler shifted wavelengths due to the total velocity.
    pub shifted_wavelength:Vec<f64>,
    /// [Vec<f64>] cointaining the flux of the observed wavelengths.  
    pub flux: Vec<f64>,
    /// [Vec<f64>] containing the flux of the continuum expectra (i.e. blackbody radiation) if the requested wavelenghts. 
    pub continuum:Vec<f64>,
}


/// This structure holds the local quantities over a surface element of the star. 
pub struct SurfaceCell{
    ///Effective temperature.
    t_eff: f64,
    ////Log value of the surface gravity.
    log_g: f64,
    ///Normalized surface area.
    area: f64,
    ///Cosine of the angle between the surface cell normal and the line of sight. 
    pub coschi: f64,
    ///Relative Doppler shift of a wavelength.
    rel_dlamb: f64,
    ///Total velocity. It is not 
    v_tot: f64,
}

///Contains the relevant information to perform interpolation and get specific intensity values. 
#[derive(Clone)]
pub struct SpectralGrid {
        /// Effective temperature of the plane parallel atmosphere.
        t_eff:[f64;2],
        /// Logarithm of the surface gravity of a plane parallel atmosphere. 
        log_g:[f64;2],
        /// Specific intensity and continuum intensity values dependant of the wavelength and χ.
        grid_values:ndarray::Array3<f64>,
        /// Array containing the wavelengths. 
        wavelengths: Vec<f64>,
        /// µ=sqrt(cos(χ)),
        /// where χ is the angle of the normal of a parallel atmosphere plane with respect to the unit vector in direction of the observer.
        mu_values:[f64;7],

        /// Important indices 
        row_indices: Vec<usize>
}

impl FluxOfSpectra{
    /// This function creates a new instance of [FluxOfSpectra] where all of its members contain values of only 0.0.
    /// This function should be used to construct a mutable instance at the beginning of the profile program. 
    pub fn new(profile_input: &ProfileConfig)->FluxOfSpectra{
        let wavelengths = profile_input.wavelength_range.get_wavelength_vector();
        let shifted_wavelengths = wavelengths.clone();
        let time = vec![0.0;wavelengths.len()];
        let flux = vec![0.0;wavelengths.len()];
        let continuum = vec![0.0;wavelengths.len()];

        FluxOfSpectra { time: time,
			wavelengths: wavelengths,
			shifted_wavelength: shifted_wavelengths,
			flux: flux,
			continuum: continuum }
    }

    /// This function sets the specific intensity flux and continuum specific intensity as 0.0, it also stores the new phase of pulsation of the calculation. 
    pub fn restart(&mut self, time_point:f64){
         self.time.fill(time_point);
         self.flux.fill(0.0);
         self.continuum.fill(0.0);
    }

    /// This function fills the `shifted_wavelength` member of [FluxOfSpectra] by multiplying the wavelength vector  with the relative doppler shift stored in a [SurfaceCell]. 
    pub fn get_doppler_shifted_wavelengths(&mut self,cell:&SurfaceCell){
        for (n,wavelength) in self.wavelengths.iter().enumerate(){
            self.shifted_wavelength[n] = wavelength * cell.rel_dlamb;
        }
    }
}


// This are the necessary parameters to run the profile program.
/// The initialization of the program requires to specify the wavelength range over which the intensity fluxes will be computed.
/// 
/// This structure also captures the maximum pulsation velocity, which is used to look within the intensity grids, and 
/// a collection of the intensity grid files that will be used for the interpolation. And the path to the intensity grids which are meant to be collected on the same folder. 
#[derive(Deserialize,Debug,PartialEq)]
pub struct ProfileConfig{
    /// This is the requested range of wavelengths.
    pub wavelength_range:WavelengthRange,
    /// This is the path to the directory containing the intensity grids
    pub path_to_grids: String,
    /// This is a [Vec] collection of [IntensityGrid]s. 
    pub intensity_grids:Vec<IntensityGrid>,
}
/// The wave length range is defined in nanometers.
/// The start should be bigger than the end and the step should be reasonable enough
#[derive(Deserialize,Debug,PartialEq)]
pub struct WavelengthRange{
    pub start: f64,
    pub end: f64,
    pub step: f64,
}
/// The intensity grids are characterized by
/// 
/// - the file name stored as a string,
/// - the temperature in Kelvin
/// - the logarithm of the surface gravity
#[derive(Deserialize,Debug,PartialEq)]
pub enum IntensityGrid{
    Joris{
        temperature: f64,
        log_gravity: f64,
        filename: String,
    },
    /// Nadya's intensity grids also enable parameter study in metalicity
    /// They provide specific intensity values dependant on µ.
    Nadya{
        temperature: f64,
        log_gravity: f64,
        metalicity: f64,
        filename:String,
    }
}

//maybe I'll use this struct CommonGridId{ temp, logg, fname}


impl ProfileConfig {
    /// This function is used to fill the parameters required for the profile program to run out of the toml configuration file.
    /// #### Arguments:
    /// * `path_to_file` - this is a string that indicates the path to the `profile_input.toml` file
    /// #### Returns:
    // * new instance of the profile config structure.
    pub fn read_from_toml(path_to_file:&str)->Self{
        let contents = match fs::read_to_string(path_to_file){
            Ok(c)=>c,
             Err(_) => { panic!("Could not read file {}",path_to_file)}
            };
        let params: ProfileConfig = match toml::from_str(&contents){
            Ok(d) => d,
            Err(e) => {println!("Unable to load data from {}",path_to_file);
                panic!("error {}",e)}
        }; 

        match params.intensity_grids_are_loaded(){
            Ok(()) => { params }
            Err(e) => { panic!("Intensity files are not properly loaded into the directory, '{}' not found.
 Either load the file into the directory or check the toml file to see if there was a misspelling",e)}
        }
    }

}

impl WavelengthRange{
    /// This method returns the wavelength vector out of  the range specified on the toml file
    /// It also checks if the step size is reasonable enough
    pub fn get_wavelength_vector(&self)->Vec<f64>{
        if self.end < self.start {panic!("Wave length range is ill defined. Please correct the toml file")}
        let capacity = ((self.end-self.start)/self.step).floor() as usize + 1usize;
        if capacity >= N_FLUX_POINTS as usize {panic!("Error, too many flux points requested.")}
	    
        let mut wavelength:Vec<f64> = Vec::with_capacity(capacity);
	    wavelength.push(self.start);
	    for i in 0..=capacity {//<- inclussive loop so wavelength[capacity]==λ_f.
	            wavelength.push( self.start + self.step * (i as f64) );
        }
        wavelength
    }
}


/// This function that inserts the relative Doppler wavelength shift to the pulstar lazyframe.
/// For a surface cell, the pulsation + rotation speed creates a doppler shift for all of the wavelengths given by 
/// delta_lambda = lambda * (1 - total_velocity/speed of light)
/// this last factor is inserted into the data frame so that it reduces the frequency of times it needs to be calculated
/// ### Arguments: 
/// * `lf` - a [LazyFrame] created out of the rasterized_star.parquet file.
/// ### Returns:
/// * a [LazyFrame] with the same headers as the `rasterized_star.parquet` with an extra column |relative shift| that holds the value
/// (1 - v/c)
pub fn insert_col_relative_dlambda(lf:LazyFrame)// I take ownership of the data frame since I will produce a new one and want the old one to be dropped after appending
->LazyFrame{
    let doppler_shift = lit(1.0) - col("velocity")/lit(CLIGHT)*lit(1.0e3);
    lf.clone().lazy().
    with_column(doppler_shift.clone().alias("relative shift"))
}

/// This function takes a polars data frame and returns all of the values from a given column that holds f64 values. 
/// ### Arguments: 
/// * `column_name` - a string slice that holds the name of a column. The column should hold f64 values.
/// * `df`- a polars DataFrame
/// ### Returns:
/// * `Vec<f64>` - a vector that contains all of the values on the column.
fn extract_column_as_vectorf64(column_name: &str,df:&DataFrame)->Vec<f64>{
    let column = df.column(column_name).unwrap();
    column.f64().unwrap().into_iter().flatten().collect()
}

/// This function takes a polars data frame and returns all of the values from a given column that holds string values. 
/// ### Arguments: 
/// * `column_name` - a string slice that holds the name of a column. The column should hold String values.
/// * `df`- a polars DataFrame
/// ### Returns:
/// * `Vec<String>` - a vector that contains all of the values on the column.
fn extract_column_as_vector_string (column_name: &str, df:&DataFrame)->Vec<String>{
    let column = df.column(column_name).unwrap();
    let vec_str:Vec<&str> = column.str().unwrap().into_iter().flatten().collect();

    let vecc:Vec<String> = vec_str.iter().map(|s| s.to_string()).collect();

    vecc
}


impl SurfaceCell {
    /// This function creates a vector collection of [SurfaceCell]s from the pulstar's output. 
    /// ### Arguments: 
    /// `star` - A [DataFrame] that contains the pulstar's output, which is a rasterized information of the parameters of the star over its surface. 
    /// ### Retruns: 
    /// - A [Vec] of [SurfaceCell]s that hold all of the local values. 
    pub fn extract_cells_from_df(star: DataFrame)-> Vec<Self>{

        let rel_dlamb_vector = extract_column_as_vectorf64("relative shift", & star);
        let area_vector = extract_column_as_vectorf64("area", & star);
        let coschi_vector = extract_column_as_vectorf64("coschi", & star);
        let temperature_vector = extract_column_as_vectorf64("temperature", & star);
        let log_g_vector = extract_column_as_vectorf64("log gravity", & star);
        let velocity_vector = extract_column_as_vectorf64("velocity", & star);

        let mut cells: Vec<SurfaceCell> = Vec::new();
        for index in 0..area_vector.len(){
            cells.push(
                SurfaceCell { t_eff: temperature_vector[index],
				log_g: log_g_vector[index],
				area: area_vector[index],
				coschi: coschi_vector[index],
				v_tot: velocity_vector[index],
				rel_dlamb: rel_dlamb_vector[index] }
            )
        }
        cells
    }
}

/// This function applies binary search on an ordered vector to find the index where a reference value
/// would be inserted. 
/// 
/// This is useful on this program because it gives the first and last values of the vector that are less than and greather than respectively.
/// ### Arguments: 
/// * `vector` - a borrowed vector of `f64` 
/// * `key` - a reference value.
/// ### Returns:
/// * `index` - a `usize` value witht the property that  `vector[index-1]<key<vector[index]`
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

/// This function returns the maximum or minimum value of a column of [f64] from a [DataFrame]. This function is adviced to be used seldomly (as in outside of loops).
/// 
/// ### Arguments:
/// * `column_name` -  The name of the column of [f64] values from which an extreme value should be extracted.
/// * `lf` - a [LazyFrame] from which the value will be extracted
/// * `is max` - a [bool] that indicates wheter to return `max`(true) or `min`(false).
/// ### Returns:
/// This function returns a [PolarsResult] with the following variants:
/// * `Ok([f64])` - where the binded value is the extremal requested one.
/// * Err(E) - returns an error to the calling function.
pub fn extremal_val_from_col(
    column_name: &str,
    lf: LazyFrame,
    is_max: bool )->PolarsResult<f64>
    {
        let df_maxval = match is_max{
            true => {lf.select([col(column_name)]).max().collect()?}
            false => {lf.select([col(column_name)]).min().collect()?}
        };
        let max_val = extract_column_as_vectorf64(column_name, &df_maxval);
        match max_val.get(0){
            Some(value) => { Ok(*value)}
            None => {Err(PolarsError::InvalidOperation
                (ErrString::new_static("Couldn't extract extreme value. The dataframe might be corrupted")
            ))}
        }
    }


