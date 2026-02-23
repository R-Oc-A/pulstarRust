use crate::{PulstarConfig, reference_frames::Coordinates, utils::{print_info::print_report, 
                    write_grid_data::{output_to_parquet, write_output}}};
use std::{time::Instant};
use crate::{AdvanceInTime,ParsingFromToml};
use polars::prelude::*;

/// This function is used to get the output of the pulstar code
pub fn pulstar_main(path:&str)->Option<DataFrame>{
    println!("--------------------");
    println!("|PULSTARust launched|");
    println!("--------------------");
    
    let now = Instant::now();
    //----------------------------------------
    //----------Read input file---------------
    //----------------------------------------

    //let path = String::from("pulstar_input.toml");
    let mut pulse_config = PulstarConfig::read_from_toml(path);
    /*let mut pulse_config:PulstarConfig = match toml::from_str(path){
        Ok(config)=>{config}
        _=>{panic!("error loading pulstar config toml file")}
    };
    */
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
    
    let mut collection_df:Option<DataFrame> = None;
    for (n,time_stamp) in time_points.iter().enumerate(){
        println!("\n +-- Computing surface data for time point number {} with time stamp {:.3}.", n,*time_stamp);

        //--Compute the time phases
        pulse_config.advance_in_time(*time_stamp);
        star.advance_in_time(*time_stamp);
        //--Initialize the minimum and maximum arrays
        
        //--Computes effective temperature, log gravity, radial component of total velocity, cosχ, etc. on all surface cells avoiding the poles.  
        star.compute_local_quantities(&pulse_config, &k);
        
        //--Save the data of the current phase.
        collection_df=Some(write_output(&star,collection_df).unwrap());
    }//end for time loop
    
    
    println!("----------------------");
    println!("|PULSTARust Finished |");
    println!("----------------------");

    collection_df

}