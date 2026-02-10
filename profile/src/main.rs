use std::env;
use profile::*;
use std::time::Instant;
use profile::profile_mkr::*;

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

   //---------------------------------------- 
   //----Parsing rasterized_star.parquet-----
   //----------------------------------------
   // Obtain the lazy frame of the parquet file, Obtain the time points, obtain the theta points
   let (lf,time_points)=parsing_star(&env_args[2].clone());
   let (
        mut spectral_grid,
        mut hypercube3d,
        mut hypercube4d,
    )= loading_intensity_grids(lf.clone(), & profile_config);
    //----------------------------------------------------------------
    //-------------- Collect fluxes for each time point  -------------
    //----------------------------------------------------------------

    //time loop    
    for (time_point_number,pulsation_phase) in time_points.iter().enumerate() {
        fluxes.integrate(
            lf.clone(),
            *pulsation_phase,
            & mut spectral_grid,
            & mut hypercube3d,
            & mut hypercube4d);
        println!("done computing flux");

        println!("finished collecting fluxes {}",pulsation_phase);
        println!("time_elapsed is {:?} seconds",start_computing_time.elapsed());
        
        fluxes.write_output(time_point_number as u16).expect(&format!("Unable to write parquet file for {} time point",*pulsation_phase));

    }
    println!("finished computation for a star's pulsation");
    println!("Total computation time is {:#?}",start_computing_time.elapsed());
}
