use std::env;
use polars::prelude::*;
//use profile::intensity::IntensityGrids;
use profile::*;
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
   //|-->Initialize the wavelength array that tracks the intensity flux profile.
    let wavelength = profile_config.wavelength_range.get_wavelength_vector();

    // // |-->Create a data frame with all of the intensity grids information
    //let grid_id_df = profile_config.get_intensity_grids_dataframe();
    // |--> Create a lazy frame that will be cloned inside the loop so that querys can be chained from here.
    //let grid_id_lf_original = grid_id_df.lazy();


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
    println!("time_elapsed is {:?} seconds",start_computing_time.elapsed());
    // Initialize the arrays that contain all the output information. This will be saved at the end of the computation on a parquet file.
    let mut all_fluxes:Vec<f64> = Vec::new(); // Intensity fluxes.
    let mut all_cont: Vec<f64> = Vec::new(); // Continuum fluxes.
    let mut all_times: Vec<f64> = Vec::new(); // time point where the intensity and continuum are calculated.
    let mut all_wavelengths: Vec<f64> = Vec::new(); // Wavelengths whose fluxes are computed.

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
    let maxval_rel_dopplershift =  (1.0+max_vel)/CLIGHT*1.0e3;
    let minval_rel_dopplershift = (1.0+min_vel)/CLIGHT*1.0e3;
    
    let intensity_dfs = profile_config.get_filtered_intensity_dataframes(
        &wavelength,
        maxval_rel_dopplershift, 
        minval_rel_dopplershift);
 
    //----------------------------------------------------------------
    //-------------- Collect fluxes for each time point  -------------
    //----------------------------------------------------------------
    //time loop    
    for phase in time_points.iter() {
        
        let expr = col("time").eq(lit(*phase));
        let sphere_frame = lf.clone().filter(expr);
        //--------------------------------------------------
        //----Collect fluxes over the whole star------------
        //--------------------------------------------------
        
        // Define the fluxes over the whole star.
        let capacity = wavelength.len();
        let mut flux = vec![0.0;capacity];
        let mut cont = vec![0.0;capacity];
        let mut time_vec = vec![*phase;capacity];

        println!("finished collecting a the star pulsation profile for the timestep {}",phase);
        println!("time_elapsed is {:?} seconds",start_computing_time.elapsed());

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

        let doppler_shift_series = observed_sphere_df.column("relative shift").unwrap();
        let area_series = observed_sphere_df.column("area").unwrap();
        let coschi_series = observed_sphere_df.column("coschi").unwrap();
        let temperature_series = observed_sphere_df.column("temperature").unwrap();
        let log_gravity_series = observed_sphere_df.column("log gravity").unwrap();
        let doppler_shift:Vec<f64> = doppler_shift_series.f64().unwrap().into_iter().flatten().collect();
        let area:Vec<f64> = area_series.f64().unwrap().into_iter().flatten().collect();
        let coschi:Vec<f64> = coschi_series.f64().unwrap().into_iter().flatten().collect();
        let temperature:Vec<f64> = temperature_series.f64().unwrap().into_iter().flatten().collect();
        let log_gravity:Vec<f64> = log_gravity_series.f64().unwrap().into_iter().flatten().collect();
        // phi loop
        for i in 0..coschi.len(){
            // Calculate the fluxes over the surface cell
            let fluxcont =return_flux_for_cell_thetaphi(
                coschi[i],
                temperature[i],
                log_gravity[i],
                doppler_shift[i],
                area[i],
                &wavelength,
                &intensity_dfs,
                );
                // Store the fluxes over the cell on vectors
                let flux_thetaphi = fluxcont.0;
                let cont_thetaphi = fluxcont.1;
                // Collect fluxes onto the global vectors
                for (index,flux_item) in flux_thetaphi.into_iter().enumerate(){
                    flux[index] += flux_item;
                    cont[index] += cont_thetaphi[index]; 
                }
                
            }
        
        all_fluxes.append(& mut flux);
        all_cont.append(& mut cont);
        all_times.append(& mut time_vec);
        all_wavelengths.extend(wavelength.iter().copied());
        
    }//end time loop
    println!("finished computation for a star's pulsation");
    
    //Write parquet file
    let output_file_name= String::from("wavelengths.parquet");
    utils::write_into_parquet(&output_file_name,
         &all_wavelengths, 
         &all_fluxes, 
         &all_cont, 
         &all_times);
    
    
}
