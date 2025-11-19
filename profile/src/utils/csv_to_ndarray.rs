use crate::IntensityGrid;
use crate::ProfileConfig;
use crate::SpectralGrid;
use crate::extract_column_as_vectorf64;
use crate::intensity::parse_intensity_grids::{
    filter_wavelength_range,
    joris_grids::convert_joris_grid_to_regular_grid};
use polars::prelude::*;
use ndarray::{Array2,Array3};

/// This module is used to obtain a [SpectralGrid] out of the Csv files provided by Nadya and Joris. 



impl IntensityGrid {
    /// This function is used to provide an schema to the csv files read from disk
    /// ### Returns:
    /// * - This function returns a [Vec] collection of the [Field]s used on [polars] to indicate the headers of the columns;
    fn get_schema(& self)->Vec<Field>{
        match self{
            // In Joris grids the first column is the wavelenght array,
            // the other 4 columns are the coefficients of the fourth order limb darkening law used to compute the specific intensity and the last four are for the continuum intensity.
            IntensityGrid::Joris { temperature:_, log_gravity:_,filename:_ }=>{
                vec![
            		Field::new("wavelengths".into(), DataType::Float64),
            		Field::new("a".into(), DataType::Float64),
            		Field::new("b".into(), DataType::Float64),
            		Field::new("c".into(), DataType::Float64),
            		Field::new("d".into(), DataType::Float64),
            		Field::new("ac".into(), DataType::Float64),
            		Field::new("bc".into(), DataType::Float64),
            		Field::new("cc".into(), DataType::Float64),
            		Field::new("dc".into(), DataType::Float64)
                ]
            }
            // In  Nadya's grids the first column is the wavelenght array,
            // the other 7 columns are the specific intensity values for  a given mu, and the last seven are for the continuum intensity for these same mu values.
            IntensityGrid::Nadya { temperature:_, log_gravity:_, metalicity:_, filename:_ }=>{
                vec![
                    Field::new("wavelengths".into(), DataType::Float64),
                    Field::new("mu1_s".into(),DataType::Float64),//specific intensity for mu=0.9636
                    Field::new("mu2_s".into(),DataType::Float64),//specific intensity for mu=0.8864
                    Field::new("mu3_s".into(),DataType::Float64),//specific intensity for mu=0.8018
                    Field::new("mu4_s".into(),DataType::Float64),//specific intensity for mu=0.7071
                    Field::new("mu5_s".into(),DataType::Float64),//specific intensity for mu=0.5976
                    Field::new("mu6_s".into(),DataType::Float64),//specific intensity for mu=0.4629
                    Field::new("mu7_s".into(),DataType::Float64),//specific intensity for mu=0.2673
                    Field::new("mu1_c".into(),DataType::Float64),//continuum intensity for mu=0.9636
                    Field::new("mu2_c".into(),DataType::Float64),//continuum intensity for mu=0.8864
                    Field::new("mu3_c".into(),DataType::Float64),//continuum intensity for mu=0.8018
                    Field::new("mu4_c".into(),DataType::Float64),//continuum intensity for mu=0.7071
                    Field::new("mu5_c".into(),DataType::Float64),//continuum intensity for mu=0.5976
                    Field::new("mu6_c".into(),DataType::Float64),//continuum intensity for mu=0.4629
                    Field::new("mu7_c".into(),DataType::Float64),//continuum intensity for mu=0.2673
                ]
            }
        }


    }


     /// This function is used to create a polars LazyFrame out of the intensity grid file in order to perform 
     /// column wise operations faster. 
     /// ### Arguments:
     /// * `path` - a reference to a string slice that contains the relative path to the intensity grid file
     /// ### Returns:
     /// * `PolarsResult<LazyFrame>` - where the lazy frame has as headers
     /// * `|wavelength|a|b|c|d|ac|bc|cc|dc|`- In case we're dealing with Joris intensity grids.
     /// * `|wavelength|I_s(mu1)|I_s(mu2)|I_s(mu3)|I_s(mu4)|I_s(mu5)|I_s(mu6)|I_s(mu7)|I_c(mu1)|I_c(mu2)|I_c(mu3)|I_c(mu4)|I_c(mu5)|I_c(mu6)|I_c(mu7)|`- In case they're Nadya's intensity grids.
     fn read_intensity_grid_file(& self,path_to_grid: &str) -> PolarsResult<LazyFrame> {
        let filename = match self {
            Self::Joris { temperature:_, log_gravity:_, filename:path_to_grid}=>{path_to_grid} 
            Self::Nadya { temperature:_, log_gravity:_, metalicity:_, filename: path_to_grid }=>{path_to_grid}
        };
        let path = format!("{}{}",path_to_grid,filename);
        let schema = Schema::from_iter(self.get_schema());
        let lf= LazyCsvReader::new(path)
        .with_separator(b' ')
        .with_has_header(false)
        .with_schema(Some(Arc::new(schema))).finish()?;
        Ok(lf)
    }


    /// This function extracts all of the csv data binded to a [IntensityGrid] into a [Vec<f64>]
    fn extract_grid_into_vector(& self, 
        wavelengths:&[f64],
        maxval_rel_dopplershift:f64,
        minval_rel_dopplershift:f64,
        path_to_grid: &str)->(usize,Vec<f64>){

        let filtered_lf = filter_wavelength_range(
            self.read_intensity_grid_file(path_to_grid).unwrap(),
            wavelengths,
            maxval_rel_dopplershift,
            minval_rel_dopplershift);
        let grid_df = match self{
                Self::Nadya { temperature:_, log_gravity:_, metalicity:_, filename:_ }=>{filtered_lf.collect().unwrap()}
                Self::Joris { temperature:_, log_gravity:_, filename:_ }=>{convert_joris_grid_to_regular_grid(filtered_lf.clone()).collect().unwrap()}
                _=>{panic!("this grid type {:?} hasn't been coded!",self)}
        };
        let mut collection:Vec<f64>= Vec::new();
            //I'm going to order the mu values from lower to greater.
        collection.append(&mut extract_column_as_vectorf64("mu7_s", &grid_df));
        let nwavelengths = collection.len();
        collection.append(&mut extract_column_as_vectorf64("mu6_s", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu5_s", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu4_s", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu3_s", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu2_s", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu1_s", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu7_c", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu6_c", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu5_c", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu4_c", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu3_c", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu2_c", &grid_df));
        collection.append(&mut extract_column_as_vectorf64("mu1_c", &grid_df));
            
        (nwavelengths,collection)
    }

    ///This function extracts the data contained in the csv file binded to a [IntensityGrid] into a 2D array
    /// ### Returns:
    /// * - A [Array2<f64>] that has the same tabular structure as the csv file of the specific intensity grid. 
    fn extract_grid_into_array2(& self,
        wavelengths:&[f64],
        maxval_rel_dopplershift:f64,
        minval_rel_dopplershift:f64,
        path_to_grid: &str
        )->Array2<f64>{
        let ncols= 14usize;
        let (nrows, collection) = self.extract_grid_into_vector(wavelengths,maxval_rel_dopplershift,minval_rel_dopplershift,path_to_grid);
        //the shape of the matrix is inverted because polars stores data column wise, and ndarray does it row wise.
        let shape = (ncols,nrows);
        //Then the array obtained used ndarray's from_shape_vec method is the transpose of the one we want
        let transposed = ndarray::Array2::from_shape_vec(shape, collection).unwrap();
        //So the final array must have the axes reversed.
        transposed.reversed_axes()
    }
}

impl ProfileConfig{
    /// This function is used to obtain a [SpectralGrid] data structure from the appropriately loaded specific intensity grids.
    /// ### Arguments:
    /// * This is an implementation on the [ProfileConfig] data structure that contains the user's inputs. 
    /// ### Returns:
    /// * This implementation returns a [SpectralGrid] that contains the domain on the parameter space as well as the tabular data of the intensity grids. 
    pub fn init_spectral_grid_from_csv(&self,
        maxval_rel_dopplershift:f64,
        minval_rel_dopplershift:f64)->SpectralGrid{

        let intensity_grids = &self.intensity_grids;
        
        //nested array2 of the data in the grid files. 
        let obs_wavelengths= self.wavelength_range.get_wavelength_vector();
        let wavelengths = self.extract_wavelength_array_from_grid(
            &obs_wavelengths,
            maxval_rel_dopplershift,
            minval_rel_dopplershift
        );
        let mut t_eff:[f64;2] = [0.0;2];
        let mut log_g:[f64;2] = [0.0;2];

        let nrows=wavelengths.len();
        let ncols= 14usize;
        
        let mut nested:Vec<Array2<f64>>=Vec::new();

        for (n,grid) in intensity_grids.iter().enumerate(){
            nested.push(grid.extract_grid_into_array2(&wavelengths,
            maxval_rel_dopplershift,
            minval_rel_dopplershift,
        &self.path_to_grids));
            if n>3{panic!("Exceeded the number of spectral grids.")}
        }
        
        let flat:Vec<f64> = nested.iter().flatten().cloned().collect();
        let shape = (4usize,nrows,ncols);        
        let array3 = Array3::from_shape_vec(shape, flat).expect("There was an error loading the spectral grids.");

        //Fill last members of coordinate space
        match self.intensity_grids[3]{
            IntensityGrid::Joris { temperature , log_gravity , filename:_ }=>{
                t_eff[1]=temperature;
                log_g[1]=log_gravity;
            }
            IntensityGrid::Nadya { temperature , log_gravity , metalicity:_ , filename:_}=>{
                t_eff[1]=temperature;
                log_g[1]=log_gravity;
            }
        }

        //Fill first members of coordinate space and create spectral grids;
        let mu_values=[0.2673,0.4629,0.5976,0.7071,0.8018,0.8864,0.9636];
        match self.intensity_grids[0]{
            IntensityGrid::Joris { temperature , log_gravity , filename:_ }=>{
                t_eff[0]=temperature;
                log_g[0]=log_gravity;
            }
            IntensityGrid::Nadya { temperature , log_gravity , metalicity:_ , filename:_  }=>{
                t_eff[0]=temperature;
                log_g[0]=log_gravity; 
            }
        }
        let row_indices = vec![0usize;2*wavelengths.len()];
        SpectralGrid{ t_eff:t_eff, log_g:log_g, grid_values: array3, wavelengths:wavelengths, mu_values:mu_values, row_indices:row_indices}
    }

    ///This function is used to obtain the wavelengths of the specific intensity grids. 
    fn extract_wavelength_array_from_grid(&self,
        obs_wavelengths:&[f64],
        maxval_rel_dopplershift:f64,
        minval_rel_dopplershift:f64)->Vec<f64>{
        let df = filter_wavelength_range(self.intensity_grids[0].read_intensity_grid_file(&self.path_to_grids).unwrap(),
                    obs_wavelengths,
                    maxval_rel_dopplershift,
                    minval_rel_dopplershift
                    )
                    .collect()
                    .unwrap();

        extract_column_as_vectorf64("wavelengths", &df)
    }
}

