use crate::{reference_frames::{Coordinates}, 
            utils::{print_info::{ print_report}, 
                    write_grid_data::write_output_to_parquet},
             PulstarConfig,};
use std::{time::Instant};
use crate::{AdvanceInTime,ParsingFromToml};

/// This function is used to get the output of the pulstar code
pub fn pulstar_main(path:&str){
    println!("--------------------");
    println!("|PULSTARust launched|");
    println!("--------------------");
    
    let now = Instant::now();
    //----------------------------------------
    //----------Read input file---------------
    //----------------------------------------

    //let path = String::from("pulstar_input.toml");
    let mut pulse_config = PulstarConfig::read_from_toml(path);
    let time_points = pulse_config.get_time_points(); 

    
    //----------------------------------------
    //---Initialize some useful parameters.---
    //----------------------------------------
    let mut star = pulse_config.rasterize_star();

    //--The components of a unit vector pointing towards the observer
    let k = Coordinates::unit_vector_k(
        pulse_config.star_data.inclination_angle.to_radians());

    //---------------------------------------- 
    //----------Start of loop-----------------
    //---------------------------------------- 

    for (n,time_stamp) in time_points.iter().enumerate(){
        println!("\n +-- Computing surface data for time point number {} with time stamp {:.3}.", n,*time_stamp);

        //--Compute the time phases
        pulse_config.advance_in_time(*time_stamp);
        star.advance_in_time(*time_stamp);
        //--Initialize the minimum and maximum arrays
        
        //--Computes effective temperature, log gravity, radial component of total velocity, cosχ, etc. on all surface cells avoiding the poles.  
        star.compute_local_quantities(&pulse_config, &k);
        
        //--Save the data of the current phase.
        write_output_to_parquet(&star, n as u16 +1).unwrap();
    }//end for time loop
    
    // Prints some values of the run
    print_report(&now, &pulse_config, time_points.len());
    
    println!("----------------------");
    println!("|PULSTARust Finished |");
    println!("----------------------");

}