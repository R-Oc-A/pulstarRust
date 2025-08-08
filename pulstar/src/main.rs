use temp_name_lib::type_def::{CYCLI2RAD};
use pulstar::{local_pulsation_velocity::observed_pulsation_velocity, local_temperature_and_gravity::local_surface_temperature_logg, reference_frames::{self, surface_normal, Coordinates}, utils::{print_info::print_report, write_grid_data::write_output_to_parquet}, PulstarConfig};
use pulstar::utils::write_grid_data::{RasterizedStarOutput, collect_output};
use std::{env, f64::consts::PI, time::Instant};
use pulstar::AdvanceInTime;
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

    //let path = String::from("pulstar_input.toml");
    let mut pulse_config = PulstarConfig::read_from_toml(path);
    let time_points = pulse_config.get_time_points(); 
    let mesh_structure = pulse_config.get_mesh_structure();
    let theta_step = mesh_structure.0;
    let phi_step = mesh_structure.1;
    println!("This is a cute little test for printing ");
    println!("the param config {:#?}",pulse_config);
    //----------------------------------------
    //----------Read input file---------------
    //----------------------------------------
    
    //----------------------------------------
    //---Initialize some useful parameters.---
    //----------------------------------------
    let mut freqrad:Vec<f64> = Vec::new();
    let mut period:Vec<f64> = Vec::new();
    let mut k_theory:Vec<f64> = Vec::new();
    let mut max_veloc:Vec<f64> = Vec::new();
    let mut min_veloc:Vec<f64> = Vec::new();
    let mut max_temp:Vec<f64> = Vec::new();
    let mut min_temp:Vec<f64> = Vec::new();
    let mut max_logg:Vec<f64> = Vec::new();
    let mut min_logg:Vec<f64> = Vec::new();
    //iterate in an idiomatic way

    //----------------------------------------
    //-----Constructing the pulsation values--
    //----------------------------------------
    for (n,mode) in pulse_config.mode_data.iter().enumerate(){
        // Frequency in rad/s
        freqrad.push(mode.frequency*CYCLI2RAD);
        // Period of pulsation in hours
        period.push(2.0*PI/3.6e3/freqrad[n]);
        // Theoretical zero order K value, mass & radius are in solar units
        if mode.l != 0u16 {
            k_theory.push(74.437 
                * pulse_config.star_data.mass 
                /pulse_config.star_data.radius.powi(3)//r^3
                /mode.frequency.powi(2));//freq^2
        }
        else{
            k_theory.push(0.0);
        }
    }

    //--The inclination angle in radians
    let incl_rad = pulse_config.star_data.inclination_angle.to_radians();
    //--Equilibrium log(g_0) (gravity g_0 is in cgs units)
    //--Mass & radius are in solar units
    let log_g0 = 4.438 + pulse_config.star_data.mass.log10()
            - 2.0 * pulse_config.star_data.radius.log10();
    let g0 = 10.0_f64.powf(log_g0);
    //--The components of a unit vector pointing towards the observer
    let k = Coordinates::unit_vector_k(incl_rad);
    //--Reset the maximum length of velocity and relative displacement vector
    
    //---------------------------------------- 
    //----------Start of loop-----------------
    //---------------------------------------- 

    for (n,time_stamp) in time_points.iter().enumerate(){
        println!("\n +-- Computing surface data for time point number {} with time stamp {:.3}.", n,*time_stamp);

        let mut star_output = RasterizedStarOutput::new_sphere(
            (180.0/theta_step) as usize, 
            (360.0/phi_step) as usize, 
            time_points.len());
        //--Compute the time phases
        pulse_config.advance_in_time(*time_stamp);
        //--Initialize the minimum and maximum arrays
        max_veloc.push(-1.0e30);
        min_veloc.push(1.0e30);
        max_temp.push(-1.0e30);
        min_temp.push(1.0e30);
        max_logg.push(-1.0e30);
        min_logg.push(1.0e30);
        
        //--Start the loop over the whole stellar surface. Avoid the poles.
        let mut theta:f64 = 1.0;
        while  theta < 180.0 {
            let theta_rad = theta.to_radians();
            let mut phi:f64 = 1.0;
         while phi <= 360.0{
                let phi_rad= phi.to_radians();

                let mut local_veloc =0.0;
                let mut local_temp =0.0;
                let mut local_logg = 0.0;
                let mut local_area = 0.0;

                let s_normal = surface_normal(&pulse_config,
                     theta_rad, phi_rad).unwrap();
                let k_spherical = k.transform(theta_rad, phi_rad);
                let cos_chi = reference_frames::cos_chi(&s_normal, &k_spherical, theta_rad, phi_rad);

                if cos_chi > 0.0{ //If the shifted mass element is visible

                    local_veloc = observed_pulsation_velocity(&pulse_config, theta_rad, phi_rad, &k).unwrap();
                    let local_values = local_surface_temperature_logg(
                        &pulse_config,
                        theta_rad, 
                        phi_rad, 
                        g0, 
                        pulse_config.star_data.effective_temperature);
                    local_temp = local_values.0;
                    local_logg = local_values.1;
                    local_area = s_normal.project_vector(&k_spherical).unwrap();
                    
                    if max_veloc[n] < local_veloc {max_veloc[n]=local_veloc}
                    if max_temp[n] < local_temp {max_temp[n]=local_temp}
                    if max_logg[n] < local_logg {max_logg[n]=local_logg}
                    if min_veloc[n] > local_veloc {min_veloc[n]=local_veloc}
                    if min_temp[n] > local_temp {min_temp[n]=local_temp}
                    if min_logg[n] > local_logg {min_logg[n]=local_logg}
                    
                }//end if coschi > 0

                let coschiout = |cos_chi: f64|{if cos_chi>0.0 {cos_chi} else { 0.0}};
                
                
                collect_output(&mut star_output,
                    theta_rad,
                    phi_rad,
					*time_stamp,
					local_veloc,
					local_temp,
					local_logg,
					coschiout(cos_chi),
					local_area);

                phi += phi_step ;
            }//end phi loop
            theta += theta_step;
        }//end theta loop
        write_output_to_parquet(star_output, n as u16 +1).unwrap();
    }//end for time loop
    //match write_output_to_parquet(&pulse_config, Star_Output){
    //    Ok(_)=>{println!("\n Returned nice and well to main function")}
    //    Err(e)=>{println!("Couldnt create data frame and panicqued in main function with error {}",e)}
    //};

    print_report(&now,//<---This gives the time of the computation
            &pulse_config ,
            &k_theory,
            &min_veloc, 
            &max_veloc, 
            &min_temp, 
            &max_temp, 
            &min_logg, 
            &max_logg, 
            &freqrad,
            &period, 
            log_g0);
}
