use std::env;
use polars::prelude::*;
use profile::intensity::IntensityGrids;
use temp_name_lib::type_def::N_FLUX_POINTS;
use profile::*;

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


   //---------------------------------------- 
   //------Parsing profile_input.toml--------
   //----------------------------------------
   // |--> Check that the toml file exists
   // |--> Check if the Profile_input.toml is well written.
   // |--> Check if the Intensity Grid files exist.
   // |--> Initialize the profile parameters.
   // let profile_config_path = eng_args[1];
   let profile_config = ProfileConfig::read_from_toml(&profile_config_path);

   let profile_config=Config{
    lambda_0:412.0,
    lambda_f:413.0,
    delta_lbd:0.01,
    v_max:1.0e2,
   };
   //|-->Initialize the wavelength array that tracks the intensity flux profile.
   //let wavelength = profile_config.get_wavelength();
   let capacity = ((profile_config.lambda_f - profile_config.lambda_0)
        /profile_config.delta_lbd).floor() as usize + 1usize ;
    if capacity >= N_FLUX_POINTS as usize {panic!("Error, too many flux points requested.")}
    let mut wavelength:Vec<f64> = Vec::with_capacity(capacity);
    wavelength.push(profile_config.lambda_0);
    for i in 0..=capacity {//<- inclussive loop so wavelength[capacity]==λ_f.
            wavelength.push(profile_config.lambda_0 
                + profile_config.delta_lbd * i as f64 );
    }

    // |-->Create a data frame with all of the intensity grids information
    let grid_id_df = IntensityGrids{
        temperatures:vec![20000.0,20000.0,21000.0,21000.0],
        log_g: vec![3.5,3.5,4.0,4.0],
        filenames: vec![String::from("name1"),String::from("name2"),String::from("name3"),String::from("name4")]
    }.create_dataframe_sorted();
    // |--> Create a lazy frame that will be cloned inside the loop so that querys can be chained from here.
    let grid_id_lf_original = grid_id_df.lazy();


   //---------------------------------------- 
   //----Parsing rasterized_star.parquet-----
   //----------------------------------------
   // Obtain the lazy frame of the parquet file, Obtain the time points, obtain the theta points
    let path= String::from("pulstar_output.parquet");
    let lf = LazyFrame::scan_parquet(path, Default::default()).unwrap();
    //get vector of time_points
    let tf = lf.clone().select([col("time").unique(),]).collect().unwrap();
    let extract_time_series = tf.column("time").unwrap();
    let time_points:Vec<f64> = extract_time_series.f64().unwrap().into_iter().flatten().collect();

    // Initialize the arrays that contain all the output information. This will be saved at the end of the computation on a parquet file.
    let mut all_fluxes:Vec<f64> = Vec::new(); // Intensity fluxes.
    let mut all_cont: Vec<f64> = Vec::new(); // Continuum fluxes.
    let mut all_times: Vec<f64> = Vec::new(); // time point where the intensity and continuum are calculated.
    let mut all_wavelengths: Vec<f64> = Vec::new(); // Wavelengths whose fluxes are computed.


    //----------------------------------------------------------------
    //-------------- Collect fluxes for each time point  -------------
    //----------------------------------------------------------------
    //time loop    
    for phase in time_points.iter() {
        
        let expr = col("time").eq(lit(*phase));
        let theta_frame = lf.clone().filter(expr);
        //--------------------------------------------------
        //----Collect fluxes over the whole star------------
        //--------------------------------------------------
        
        // Define the fluxes over the whole star.
        let mut flux = vec![0.0;capacity];
        let mut cont = vec![0.0;capacity];
        let mut time_vec = vec![*phase;capacity];

        //get vector of theta_points
        let theta_df = theta_frame.clone().select([col("theta").unique()]).collect().unwrap();
        let theta_steps_serie = theta_df.column("theta").unwrap();
        let theta_steps:Vec<f64> = theta_steps_serie.f64().unwrap().into_iter().flatten().collect();

        //theta loop
        for theta_step in theta_steps.iter(){

            // create a lazy frame for a latitudinal strip of the rasterized star (constant theta over a sphere)
            let expr: Expr = col("theta").eq(lit(*theta_step));
            let phi_frame=theta_frame.clone().filter(expr);
            
            // Filter if surface cell is visible.
            let expr = col("coschi").gt(lit(0.0));
            let pf_visible =phi_frame.filter(expr);
            
            // Append relative doppler wavelength shift 
            let pf = append_doppler_shift(pf_visible).collect().unwrap();

            // Obtain the relevant quantities to compute the flux on each cell of the surface of the rasterized star
            // |--> relative doppler wavelength shift
            // |--> normalized area of each cell projected onto the unit vector of directed towards the observer
            // |--> coschi is projection of the unit vector normal to the cell surface towards the observer.
            // |--> temperature over the surface cell
            // |--> log gravity value over the surface cell

            let doppler_shift_series = pf.column("doppler shift").unwrap();
            let area_series = pf.column("area").unwrap();
            let coschi_series = pf.column("coschi").unwrap();
            let temperature_series = pf.column("temperature").unwrap();
            let log_gravity_series = pf.column("log gravity").unwrap();
            let doppler_shift:Vec<f64> = doppler_shift_series.f64().unwrap().into_iter().flatten().collect();
            let area:Vec<f64> = area_series.f64().unwrap().into_iter().flatten().collect();
            let coschi:Vec<f64> = coschi_series.f64().unwrap().into_iter().flatten().collect();
            let temperature:Vec<f64> = temperature_series.f64().unwrap().into_iter().flatten().collect();
            let log_gravity:Vec<f64> = log_gravity_series.f64().unwrap().into_iter().flatten().collect();
            
            // phi loop
            for i in 0..coschi.len(){
                // Clone the lazy frame as the program will search on the Intensity grids data base
                let grid_id_lf = grid_id_lf_original.clone();
                // Calculate the fluxes over the surface cell
                let fluxcont =return_flux_for_cell_thetaphi(
                    coschi[i],
                    temperature[i],
                    log_gravity[i],
                    doppler_shift[i],
                    area[i],
                    &wavelength,
                    grid_id_lf);
                // Store the fluxes over the cell on vectors
                let flux_thetaphi = fluxcont.0;
                let cont_thetaphi = fluxcont.1;
                // Collect fluxes onto the global vectors
                for (index,flux_item) in flux_thetaphi.into_iter().enumerate(){
                    flux[index] += flux_item;
                    cont[index] += cont_thetaphi[index]; 
                }

            }//end phi loop
        }//end theta loop

        all_fluxes.append(& mut flux);
        all_cont.append(& mut cont);
        all_times.append(& mut time_vec);
        all_wavelengths.extend(wavelength.iter().copied());
    }//end time loop
    
    //Write parquet file
    let output_file_name= String::from("wavelengths.parquet");
    utils::write_into_parquet(&output_file_name,
         &all_wavelengths, 
         &all_fluxes, 
         &all_cont, 
         &all_times);

    
}
