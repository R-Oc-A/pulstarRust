use temp_name_lib::joris_math::geometry::surface_normal::surface_normal;
use pulstar::utils::{self,parse_file};
use nalgebra as na;
use std::{env,fs};
fn main() {

    // program start!

    let env_args: Vec<String> = env::args().collect();
    
    //having the right number of arguments
    if ( env_args.len() != 3usize){
        panic!("USAGE: pulstar <parameter file> <time points file>");
    }
    let parameter_file = &env_args[1];
    let time_points_file = &env_args[2];

    println!("--------------------");
    println!("|PULSTARust launced|");
    println!("--------------------");

    //Read input file
    let build_pulse_config= parse_file::parse_from_file(parameter_file);

    let time_points = parse_file::parse_timepoints(&time_points_file,&build_pulse_config);
    //asserting 


    //initialize param

    //start of loop
    /*let mut theta_index:usize= 1;
    'loop_on_theta: loop{
        let mut phi_index:usize=1;
        'loop_on_phi : loop{
            phi_index += 1;
            if phi_index > 359 {break 'loop_on_phi};

            let theta_rad = theta as f64 * pulstar::DEG2RAD;
            let phi_rad = phi as f64 * pulstar::DEG2RAD;

            let surf_normal = surface_normal(&parameters,
                                            theta_rad,
                                            phi_rad,
                                            false).unwrap_or_else
                                            (|error| "problem with theta too small {error:?}");

            let coschi = cos_chi
                        (&surf_normal,
                       cartesian_to_sphere(k,
                                          theta_rad,
                                          phi_rad));
        }
        theta_index += 1;
        if theta_index > 179 {break 'loop_on_theta};
    }
    */
    println!("Hello, world!");
}
