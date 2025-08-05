
use polars::{error::ErrString, prelude::*};
use super::*;
mod define_parameter_space;
use std::fs::File;

/// This module contains the methods and functions to process the intensity grids. 
/// It has three parts.
/// 
/// The first part consist of methods to produce a [DataFrame] that stores the file name of the relevant
/// intensity grid files
/// 
/// The second part consist of methods to obtain a [DataFrame]s from the contents of the intensity grid files 
/// 
/// The third part has the methods used to filter and compute intensity and continuum energy fluxes out the intensity [DataFrame]s
/// 
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

    /// This function creates an [Expr] that filters only the relevant grid files
    /// 
    /// ### Arguments:
    /// * `temperatures` - a reference to a vector that contains all of the local temperatures for each surface cell
    /// * `log_g` - a reference to a vector that contains all of the local values for log_g
    /// * `delta_temp` - temperature spacing between grid files
    /// * `delta_logg` - logg spacing between grid files
    /// ### Returns:
    /// This function returns a [Option] with the following variants:
    /// * `Ok(Expr)` - where Expr defines the relevant temperatures and log g in the parameter space.
    /// * 'None' - when there's an error on the computation of the local values. 
    fn filter_grids(temperatures:&[f64],log_g:&[f64],delta_temp:f64,delta_logg:f64)->Option<Expr>{
        
        let max_temp = temperatures.iter().fold(0.0, |acc,x| if acc<*x {*x}else{acc});
        let max_logg = log_g.iter().fold(0.0, |acc,x| if acc<*x {*x}else{acc});
        let min_temp = temperatures.iter().fold(0.0, |acc,x| if acc>*x {*x}else{acc});
        let min_logg = log_g.iter().fold(0.0, |acc,x| if acc>*x {*x}else{acc});

        if max_temp*max_logg*min_temp*min_logg == 0.0{ //<--- If there's some maximum or minimum missing, filter nothing.
            None
        }else{

            let upper_bound_temp = col("temperature").lt_eq(lit(max_temp+delta_temp));
            let upper_bound_logg= col("log_gravity").lt_eq(lit(max_logg+delta_logg));
            let lower_bound_temp = col("temperature").gt_eq(lit(min_temp-delta_temp));
            let lower_bound_logg= col("log_gravity").gt_eq(lit(min_logg-delta_logg));
        
            let uu_expr = upper_bound_temp.clone().and(upper_bound_logg.clone());
            let ul_expr = upper_bound_temp.and(lower_bound_logg.clone());
            let lu_expr = lower_bound_temp.clone().and(upper_bound_logg);
            let ll_expr = lower_bound_temp.and(lower_bound_logg);
        
            let combined_expr = uu_expr.and(lu_expr.and(ul_expr.and(ll_expr)));
        
            Some(combined_expr)
        }
        
    }
    /// This function creates a data frame with the columns "temperature","log_gravity","file name"
    /// that are Float64, Float64, &[str] in that order. 
    /// 
    /// The data frame is ordered from lower to greater both in temperature and gravity. 
    /// Only the relevant grid files from the computation will be loaded
    fn load_gridfiles_into_dataframe(self,temperatures:&[f64],log_g:&[f64])->PolarsResult<DataFrame>{
        let db_temp = self.temperatures.clone();
        let db_log_g = self.log_g.clone();
        let delta_temp = (db_temp.last().unwrap() - db_temp.get(0).unwrap())/(db_temp.len() as f64);
        let delta_logg = (db_log_g.last().unwrap() - db_log_g.get(0).unwrap())/(db_log_g.len() as f64);
        
        let df = self.create_dataframe_sorted();

        let lf = df.lazy();
        
        if let Some(param_space_expression) = Self::filter_grids(temperatures, log_g, delta_temp, delta_logg){
            lf.filter(param_space_expression).collect()
        }else { Err(PolarsError::InvalidOperation(ErrString::new_static("The loaded grid files are not adequate to carry on the computation."))) }
    }
}


//----------------------------------------
//------parsing intensity grids-----------
//---------------------------------------- 
pub mod parse_intensity_grids;

//--------------------------------------------------
//-----Extracting Intensity and continuum fluxes----
//--------------------------------------------------
pub mod extract_intensity_fluxes;