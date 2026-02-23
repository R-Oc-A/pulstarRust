use core::time;

use crate::utils::IntensityFlux;

use super::*;
use temp_name_lib::interpolation::ParameterSpaceHypercube;

pub fn parsing_star(path_to_star:&str)->(LazyFrame,Vec<f64>){
   //---------------------------------------- 
   //----Parsing rasterized_star.parquet-----
   //----------------------------------------
   // Obtain the lazy frame of the parquet file, Obtain the time points, obtain the theta points
    let rasterized_star_path = path_to_star;
    let lf = LazyFrame::scan_parquet(rasterized_star_path, Default::default()).unwrap();
    //get vector of time_points
    
    let tf = lf.clone().select([col("time").unique(),]).collect().unwrap();
    let extract_time_series = tf.column("time").unwrap();
    let time_points:Vec<f64> = extract_time_series.f64().unwrap().into_iter().flatten().collect();
    (lf.clone(),time_points)
}

pub fn loading_intensity_grids(star_lf:LazyFrame,
profile_config:& ProfileConfig)->(
SpectralGrid,//SpectralGrid
ParameterSpaceHypercube,//hypercube3d
ParameterSpaceHypercube,//hypercube4d
){
    let max_vel = extremal_val_from_col(
        "velocity",
         star_lf.clone(),
          true).unwrap();
    let min_vel = extremal_val_from_col(
        "velocity",
         star_lf.clone(),
          false).unwrap();
    
    let maxval_rel_dopplershift =  1.0+max_vel/CLIGHT*1.0e3;
    let minval_rel_dopplershift = 1.0+min_vel/CLIGHT*1.0e3;
    println!("min relative dopplershift is {}",minval_rel_dopplershift);
    println!("max relative dopplershift is {}", maxval_rel_dopplershift);

    println!("creating the spectral grids data structures from csv files...or neural network regresor");
    let spectral_grids = profile_config.init_spectral_grid_from_csv(maxval_rel_dopplershift, minval_rel_dopplershift);
    println!("allocating memory for hypercube in the parameter space");
    let hypercube4d= spectral_grids.new_hypercube(4usize);
    let hypercube3d= spectral_grids.new_hypercube(3usize);

    (spectral_grids,hypercube3d,hypercube4d)
}

impl FluxOfSpectra {
    pub fn integrate(& mut self,
        star_lf:LazyFrame,
        pulsation_phase:f64,
        spectral_grid:& mut SpectralGrid,
        hypercube3d:& mut ParameterSpaceHypercube,
        hypercube4d:& mut ParameterSpaceHypercube){
        let expr = col("time").eq(lit(pulsation_phase));
        let sphere_frame = star_lf.clone().filter(expr);
        //--------------------------------------------------
        //----Collect fluxes over the whole star------------
        //--------------------------------------------------
    
        // Filter if surface cell is visible.
        let expr = col("coschi").gt(lit(0.08));//.and(col("coschi").lt(lit(0.9285)));
        let visible_lf =sphere_frame.filter(expr);
            
        // Append relative doppler wavelength shift 
        let observed_sphere_df = insert_col_relative_dlambda(visible_lf).collect().unwrap();
    
        // Obtain the relevant quantities to compute the flux on each cell of the surface of the rasterized star
        // |--> relative doppler wavelength shift
        // |--> normalized area of each cell projected onto the unit vector of directed towards the observer
        // |--> coschi is projection of the unit vector normal to the cell surface towards the observer.
        // |--> temperature over the surface cell
        // |--> log gravity value over the surface cell
        let surface_cells = SurfaceCell::extract_cells_from_df(observed_sphere_df);
    
        // Integrate specific intensity.        
        self.restart(pulsation_phase);
        for cell in surface_cells.iter(){
            self.get_doppler_shifted_wavelengths(cell);
            match cell.coschi>0.9285{
                true => {self.collect_flux_from_cell(cell,  spectral_grid, hypercube3d)}
                false => {self.collect_flux_from_cell(cell,  spectral_grid, hypercube4d)}
            }
        }
    
    }

    /// So far I've only coded the version to write into a parquet file. 
    pub fn write_output(&mut self,time_point:u16)->PolarsResult<()>{
        utils::write_into_parquet(time_point + 1, self.clone())
    }
}

pub fn profile_main(toml_string:&str,star_df:DataFrame)->DataFrame{
   //---------------------------------------- 
   //------Parsing profile_input.toml--------
   //----------------------------------------
   // |--> Check that the toml file exists
   // |--> Check if the Profile_input.toml is well written.
   // |--> Check if the Intensity Grid files exist.
   // |--> Initialize the profile parameters.
    //let profile_config = ProfileConfig::read_from_toml(toml_string);
    let profile_config:ProfileConfig=match toml::from_str(toml_string){
        Ok(config)=>{config}
        _=>{panic!("error parsing toml for profile config")}
    };
    let mut fluxes = FluxOfSpectra::new(&profile_config);
    let mut intensity_collection = IntensityFlux::new();
   //---------------------------------------- 
   //----Parsing rasterized_star.parquet-----
   //----------------------------------------
    let lf = star_df.lazy();
    let tf = lf.clone().select([col("time").unique(),]).collect().unwrap();
    let extract_time_series = tf.column("time").unwrap();
    let time_points:Vec<f64> = extract_time_series.f64().unwrap().into_iter().flatten().collect();
   // Obtain the lazy frame of the parquet file, Obtain the time points, obtain the theta points
   let (
        mut spectral_grid,
        mut hypercube3d,
        mut hypercube4d,
    )= loading_intensity_grids(lf.clone(), & profile_config);
    //----------------------------------------------------------------
    //-------------- Collect fluxes for each time point  -------------
    //----------------------------------------------------------------

    //time loop    
    let mut last_timepoint= 0u16;
    for (time_point_number,pulsation_phase) in time_points.iter().enumerate() {
        fluxes.integrate(
            lf.clone(),
            *pulsation_phase,
            & mut spectral_grid,
            & mut hypercube3d,
            & mut hypercube4d);
        println!("done computing flux");

        println!("finished collecting fluxes {}",pulsation_phase);
        //fluxes.write_output(time_point_number as u16).expect(&format!("Unable to write parquet file for {} time point",*pulsation_phase));
        intensity_collection = intensity_collection.append_fluxes(fluxes.clone());
        last_timepoint=time_point_number as u16;
    }
    intensity_collection.data_frame
    //if let Ok(_)= intensity_collection.write_output(last_timepoint){
    //println!("finished computation for a star's pulsation")}
    //else{panic!("unable to write parquetfile")};
}
