use temp_name_lib::type_def::{PHI_STEP, THETA_STEP,CYCLI2RAD,RADIUSSUN};
use pulstar::{local_pulsation_velocity::observed_pulsation_velocity, local_temperature_and_gravity::local_surface_temperature_logg, reference_frames::{self, ssurface_normal, Coordinates}, utils::{parse_file, print_info::print_report, write_grid_data::write_output_to_parquet}, PPulstarConfig};
use pulstar::utils::write_grid_data::{RasterizedStarOutput, collect_output};

use std::{env, f64::consts::PI, time::Instant};
fn main() {

    // program start!

    let env_args: Vec<String> = env::args().collect();
    
    //having the right number of arguments
    if env_args.len() != 3usize {
        panic!("USAGE: pulstar <parameter file> <time points file>");
    }

        

    let parameter_file = &env_args[1];
    let time_points_file = &env_args[2];

    println!("--------------------");
    println!("|PULSTARust launched|");
    println!("--------------------");
    
    let now = Instant::now();

    let ppath = String::from("pulstar_input.toml");
    let ppulse_config = PPulstarConfig::read_from_toml(&ppath);
    let ttime_points = ppulse_config.get_time_points(); 
    let mesh_structure = ppulse_config.get_mesh_structure();
    let theta_step = mesh_structure.0;
    let phi_step = mesh_structure.1;
    println!("This is a cute little test for printing ");
    println!("the param config {:#?}",ppulse_config);
    //----------------------------------------
    //----------Read input file---------------
    //----------------------------------------
    let pulse_config= parse_file::parse_from_file(parameter_file);

    let time_points = parse_file::parse_timepoints(&time_points_file,&pulse_config);
    
    //----------------------------------------
    //---Initialize some useful parameters.---
    //----------------------------------------
    let mut freqrad:Vec<f64> = Vec::new();
    let mut period:Vec<f64> = Vec::new();
    let mut k_theory:Vec<f64> = Vec::new();
    let mut vampl:Vec<f64> = Vec::new();
    let mut t_phase_rad:Vec<f64> = Vec::new();
    let mut g_phase_rad:Vec<f64> = Vec::new();
    let mut max_veloc:Vec<f64> = Vec::new();
    let mut min_veloc:Vec<f64> = Vec::new();
    let mut max_temp:Vec<f64> = Vec::new();
    let mut min_temp:Vec<f64> = Vec::new();
    let mut max_logg:Vec<f64> = Vec::new();
    let mut min_logg:Vec<f64> = Vec::new();
    let mut maxvel_length =0.0;
    let mut maxrel_length =0.0;
    //iterate in an idiomatic way

    //----------------------------------------
    //-----Constructing the pulsation values--
    //----------------------------------------
    for (n,mode) in ppulse_config.mode_data.iter().enumerate(){
        // Frequency in rad/s
        freqrad.push(mode.frequency*CYCLI2RAD);
        // Period of pulsation in hours
        period.push(2.0*PI/3.6e3/freqrad[n]);
        // Theoretical zero order K value, mass & radius are in solar units
        if mode.l != 0u16 {
            k_theory.push(74.437 
                * ppulse_config.star_data.mass 
                /ppulse_config.star_data.radius.powi(3)//r^3
                /mode.frequency.powi(2));//freq^2
        }
        else{
            k_theory.push(0.0);
        }
        
        //--Convert the phase difference of the effective temperature from degrees to radians
        t_phase_rad.push(mode.phase_rel_dtemp.to_radians());
        //--Convert the phase difference of the effective gravity from degrees to radians
        g_phase_rad.push(mode.phase_rel_dg.to_radians());
    }

    //--The inclination angle in radians
    let incl_rad = ppulse_config.star_data.inclination_angle.to_radians();
    //--Equilibrium log(g_0) (gravity g_0 is in cgs units)
    //--Mass & radius are in solar units
    let log_g0 = 4.438 + ppulse_config.star_data.mass.log10()
            - 2.0 * ppulse_config.star_data.radius.log10();
    let g0 = 10.0_f64.powf(log_g0);
    //--The components of a unit vector pointing towards the observer
    let kk = Coordinates::unit_vector_k(incl_rad);
    //--Reset the maximum length of velocity and relative displacement vector
    
    //---------------------------------------- 
    //----------Start of loop-----------------
    //---------------------------------------- 

    for (n,time_stamp) in ttime_points.iter().enumerate(){
        println!("\n +-- Computing surface data for time point number {} with time stamp {:.3}.", n,*time_stamp);

        let mut star_output = RasterizedStarOutput::new_sphere(
            (180.0/theta_step) as usize, 
            (360.0/phi_step) as usize, 
            ttime_points.len());
        //Collecting the phase_offset of all the modes. 
        let mut phase_offset:Vec<f64> = Vec::new();
        for mode in ppulse_config.mode_data.iter(){
            phase_offset.push(2.0 * PI * (mode.phase_offset + 
                mode.frequency * *time_stamp) )//<- [Ricardo]: Here I'm dereferencing to obtain the f64 value pointed by time stamp
                //[Ricardo]: I'm seriously considering iterating time_stamps by into_iter() so I don't have to  dereference.
        }
        
        //--Initialize the minimum and maximum arrays
        max_veloc.push(-1.0e30);
        min_veloc.push(1.0e30);
        max_temp.push(-1.0e30);
        min_temp.push(1.0e30);
        max_logg.push(-1.0e30);
        min_logg.push(1.0e30);
    
        //--Start the loop over the whole stellar surface. Avoid the poles.
        let mut theta:u16 = 1;
        'theta_loop: while  theta < 180 {
            let theta_rad = (theta as f64).to_radians();
            let mut phi:u16 = 1;
            'phi_loop: while phi <= 360{
                let phi_rad= (phi as f64).to_radians();

                let mut local_veloc =0.0;
                let mut local_temp =0.0;
                let mut local_logg = 0.0;
                let mut local_area = 0.0;

                let s_normal = ssurface_normal(&ppulse_config,
                     theta_rad, phi_rad).unwrap();
                let k_spherical = kk.transform(theta_rad, phi_rad);
                let cos_chi = reference_frames::cos_chi(&s_normal, &k_spherical, theta_rad, phi_rad);

                if cos_chi > 0.0{ //If the shifted mass element is visible

                    local_veloc = observed_pulsation_velocity(&ppulse_config, theta_rad, phi_rad, &kk).unwrap();
                    let local_values = local_surface_temperature_logg(
                        &ppulse_config,
                        theta_rad, 
                        phi_rad, 
                        g0, 
                        ppulse_config.star_data.effective_temperature);
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
                else{
                    //[Ricardo]: already nullified local values, so do nothing...
                }//end else coschi > 0   

                //writing data into pulstar.out files


                let mut vel_length =0.0;
                let mut rel_lenght = 0.0;                
                /*if pulse_config.print_amplitude {
                    let sintheta = theta_rad.sin();
                    let costheta = theta_rad.cos();
                    for (index,rel_deltar) in pulse_config.mode_config.rel_deltar.iter().enumerate(){
                        vel_length += v_puls(&pulse_config.mode_config,
                            index, 
                            sintheta, 
                            costheta, 
                            phi_rad,
                            &vampl).unwrap().coords.norm();
                        
                        let ampl_radial= rel_deltar 
                                        * ylmnorm(pulse_config.mode_config.l[index], pulse_config.mode_config.m[index]);
                        let ampl_tangential= *rel_deltar
                                        * pulse_config.mode_config.k[index] * 
                                        ylmnorm(pulse_config.mode_config.l[index], pulse_config.mode_config.m[index]);
                        
                        rel_lenght += displacement(&pulse_config.mode_config,
                            index,
                            sintheta,
                            costheta,
                            phi_rad,
                            ampl_radial,
                            ampl_tangential).unwrap().coords.norm();
                    }
                    if rel_lenght > maxrel_length { maxrel_length=rel_lenght}
                    if vel_length > maxvel_length { maxvel_length=vel_length}    
                }
                */

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

                phi += PHI_STEP;
            }//end phi loop
            theta += THETA_STEP;
        }//end theta loop
        write_output_to_parquet(star_output, n as u16 +1).unwrap();
    }//end for time loop
    //match write_output_to_parquet(&pulse_config, Star_Output){
    //    Ok(_)=>{println!("\n Returned nice and well to main function")}
    //    Err(e)=>{println!("Couldnt create data frame and panicqued in main function with error {}",e)}
    //};

    print_report(&now,//<---This gives the time of the computation
            &ppulse_config ,
            &k_theory,
            &min_veloc, 
            &max_veloc, 
            &min_temp, 
            &max_temp, 
            &min_logg, 
            &max_logg, 
            &freqrad,
            &period, 
            maxvel_length, 
            maxrel_length,
            log_g0);
}
