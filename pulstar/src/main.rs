use temp_name_lib::type_def::{CYCLI2RAD};
use pulstar::{local_pulsation_velocity::observed_pulsation_velocity, local_temperature_and_gravity::local_surface_temperature_logg, reference_frames::{self, surface_normal, Coordinates}, utils::{print_info::print_report, write_grid_data::write_output_to_parquet}, PulstarConfig, RasterizedStar};
use pulstar::utils::write_grid_data::{RasterizedStarOutput};
use std::{env, f64::consts::PI, time::Instant};
use pulstar::{AdvanceInTime,ParsingFromToml};
fn main() {

    // program start!
    let env_args: Vec<String> = env::args().collect();
    
    //having the right number of arguments
    if env_args.len() != 2usize {
        panic!("USAGE: pulstar <parameter file>");
    }

    let path = &env_args[1];
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
        
        //--Start the loop over the whole stellar surface. Avoiding the poles
        star.compute_local_quantities(&pulse_config, &k);

        write_output_to_parquet(&star, n as u16 +1).unwrap();
    }//end for time loop
    println!("\n PULSTARust finished...")    
}