use std::env;
use polars::prelude::*;
//use profile::intensity::IntensityGrids;
use profile::{utils::{write_into_parquet}, *};
use temp_name_lib::type_def::CLIGHT;
use std::time::Instant;

fn main() {

   // First receive the arguments for running the program. 
   // Important to check:
   // There's a correct number of arguments. 
   // 1) The profile file
   // 2) The profile_input.toml file
   // 3) The rasterized star parquet file.
   let env_args:Vec<String> = env::args().collect(); 
   if env_args.len() < 3usize {
    println!("Not enough arguments.");
    panic!("USAGE: profile -- profile_input.toml rasterized_star.parquet");
   }


   //--------------------------------------------------
   //---------Program Start!---------------------------
   //--------------------------------------------------
    let start_computing_time = Instant::now();   



   //---------------------------------------- 
   //------Parsing profile_input.toml--------
   //----------------------------------------
   // |--> Check that the toml file exists
   // |--> Check if the Profile_input.toml is well written.
   // |--> Check if the Intensity Grid files exist.
   // |--> Initialize the profile parameters.
   let profile_config_path = env_args[1].clone();
   let profile_config = ProfileConfig::read_from_toml(&profile_config_path);
   
    let mut fluxes = FluxOfSpectra::new(&profile_config);

   //|-->Initialize the wavelength array that tracks the intensity flux profile.
    //let wavelength = profile_config.wavelength_range.get_wavelength_vector();

   //---------------------------------------- 
   //----Parsing rasterized_star.parquet-----
   //----------------------------------------
   // Obtain the lazy frame of the parquet file, Obtain the time points, obtain the theta points
    let rasterized_star_path = env_args[2].clone();
    let lf = LazyFrame::scan_parquet(rasterized_star_path, Default::default()).unwrap();
    //get vector of time_points
    
    let tf = lf.clone().select([col("time").unique(),]).collect().unwrap();
    let extract_time_series = tf.column("time").unwrap();
    let time_points:Vec<f64> = extract_time_series.f64().unwrap().into_iter().flatten().collect();
    println!("Finish getting time_points,");
    println!("here they are:");
    println!("{:#?}",time_points);
    println!("time_elapsed is {:?} seconds \n",start_computing_time.elapsed());

    //--------------------------------------- 
    //----Loading intensity flux grids-------
    //---------------------------------------

    let max_vel = extremal_val_from_col(
        "velocity",
         lf.clone(),
          true).unwrap();
    let min_vel = extremal_val_from_col(
        "velocity",
         lf.clone(),
          false).unwrap();
    
    let maxval_rel_dopplershift =  1.0+max_vel/CLIGHT*1.0e3;
    let minval_rel_dopplershift = 1.0+min_vel/CLIGHT*1.0e3;
    println!("min relative dopplershift is {}",minval_rel_dopplershift);
    println!("max relative dopplershift is {}", maxval_rel_dopplershift);

    println!("creating the spectral grids data structures from csv files...or neural network regresor");
    let mut spectral_grids = profile_config.init_spectral_grid_from_csv(maxval_rel_dopplershift, minval_rel_dopplershift);
    println!("allocating memory for hypercube in the parameter space");
    let mut hypercube= spectral_grids.new_hypercube();
    //----------------------------------------------------------------
    //-------------- Collect fluxes for each time point  -------------
    //----------------------------------------------------------------


    //time loop    

    for (time_point_number,phase) in time_points.iter().enumerate() {
        
        let expr = col("time").eq(lit(*phase));
        let sphere_frame = lf.clone().filter(expr);
        //--------------------------------------------------
        //----Collect fluxes over the whole star------------
        //--------------------------------------------------

        // Filter if surface cell is visible.
        let expr = col("coschi").gt(lit(0.0));
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

        println!("finished collecting a the star pulsation profile for the timestep {}",phase);
        println!("time_elapsed is {:?} seconds",start_computing_time.elapsed());

        // Integrate specific intensity.        
        fluxes.restart(*phase);
        for cell in surface_cells.iter(){
            fluxes.collect_flux_from_cell(cell, & mut spectral_grids,&mut hypercube);
        }
        println!("done computing flux");

        println!("finished collecting fluxes {}",phase);
        println!("time_elapsed is {:?} seconds",start_computing_time.elapsed());
        
        write_into_parquet(time_point_number as u16 + 1, fluxes.clone()).expect(&format!("Unable to write parquet file for {} time point",*phase));
    }
    println!("finished computation for a star's pulsation");
    println!("Total computation time is {:#?}",start_computing_time.elapsed());
}
