use temp_name_lib::{
    joris_math::{
        geometry::{
            displacement::displacement, 
            projections::{cos_chi, project_vector}, 
            surface_normal::surface_normal
            }, 
        ref_frame_convrs::{cartesian_to_spherical, unit_vector_k},
        spherical_harmonics::norm_factor::ylmnorm
        }, 
    star_physics::{
        local_values::local_value, 
        velocity::{
            pulsation::v_puls, 
            velocity_projection::{project_vpuls, project_vrot}
            }
        }, type_def::{PHI_STEP, THETA_STEP}
    };
use pulstar::{utils::{parse_file, print_info::print_report, write_grid_data::write_output_to_parquet},
            DEG2RAD, PI, RADIUSSUN};
use pulstar::utils::write_grid_data::{RasterStarOutput, collect_output};

use std::{env,
        time::Instant};
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
    
    let now = Instant::now();

    
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
    for (n,l_val) in pulse_config.mode_config.l.iter().enumerate(){
        //--Frequency in rad/s--
        freqrad.push(pulse_config.freqcycli[n]*pulstar::CYCLI2RAD);
        //--the period of pulsation in hours
        period.push(2.0 * pulstar::PI/3.6e3/freqrad[n]);
        //--Theoretical zero order K value, mass & radius are in solar units.
        if *l_val != 0u16 {
            k_theory.push(74.437 
                * pulse_config.star_config.mass 
                /pulse_config.star_config.radius.powi(3)//r^3
                /pulse_config.freqcycli[n].powi(2));//freq^2
        }
        else {
            k_theory.push(0.0);
        }
        //--Amplitude of the pulse velocity (/1000.0 because from m-> km)
        vampl.push(pulse_config.star_config.radius
            * freqrad[n] * RADIUSSUN * 1.0e-3
            * pulse_config.mode_config.rel_deltar[n]
        );
        //--Convert the phase difference of the effective temperature from degrees to radians
        t_phase_rad.push(pulse_config.temperature_config[n].phasedif * DEG2RAD);
        //--Convert the phase difference of the effective gravity from degrees to radians
        g_phase_rad.push(pulse_config.gravity_config[n].phasedif * DEG2RAD);
    } 
    //--The inclination angle in radians
    let incl_rad = pulse_config.star_config.inclination_angle as f64 * DEG2RAD;
    //--Equilibrium log(g_0) (gravity g_0 is in cgs units)
    //--Mass & radius are in solar units
    let log_g0 = 4.438 + (pulse_config.star_config.mass).log10()
            - 2.0 * (pulse_config.star_config.radius).log10();
    //--The components of a unit vector pointing towards the observer
    let k = unit_vector_k(incl_rad);//cartesian
    //--Reset the maximum length of velocity and relative displacement vector
    let mut maxrellength =0.0;
    let mut maxvelolength = 0.0;
    
    //---------------------------------------- 
    //----------Start of loop-----------------
    //---------------------------------------- 
    //idiomatic loop
    let mut Star_Output= RasterStarOutput::new_sphere(
        180/THETA_STEP as usize, 
        360/PHI_STEP as usize, 
        pulse_config.time_pts_nmbr as usize);

    for (n,time_stamp) in time_points.iter().enumerate(){
        println!("\n +-- Computing surface data for time point number {} with time stamp {}.", n,*time_stamp);

        //idiomatic loop
        let mut phase_offset:Vec<f64> = Vec::new();
        for (m,phase_val) in pulse_config.mode_config.phase.iter().enumerate(){
            phase_offset.push(2.0 * (*phase_val //here I'm dereferencing to obtain the phase val
                + pulse_config.freqcycli[m] * *time_stamp)//here I'm dereferencing to obtain time value
                * PI );
        }//end for loop n_modes

        // initialize files to be created.[Ricardo]: This might not be necessary on Rust.
        // initialize format conventions. [Ricardo]: This might not be necessary on Rust.

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
            let theta_rad = theta as f64 * DEG2RAD;
            let mut phi:u16 = 1;
            'phi_loop: while phi <= 360{
                let phi_rad= phi as f64 * DEG2RAD;

                let mut local_veloc =0.0;
                let mut local_temp =0.0;
                let mut local_logg = 0.0;
                let mut local_area = 0.0;

                let s_normal = surface_normal(&pulse_config.mode_config,
                     theta_rad, 
                     phi_rad,
                    pulse_config.is_time_dependent).unwrap();
                let k_spherical = cartesian_to_spherical(&k, 
                    theta_rad, phi_rad).unwrap();
                let cos_chi = cos_chi(&s_normal,
                     &k_spherical).unwrap();
                if cos_chi > 0.0{ //If the shifted mass element is visible
                    match pulse_config.suppress_pulse{
                        true => {
                            local_veloc = project_vrot(pulse_config.star_config.rotation_velocity,
                                                    theta_rad,
                                                    phi_rad, &k).unwrap();
                        }
                        false => {
                            local_veloc = project_vpuls(&pulse_config.mode_config,
                                 theta_rad,
                                 phi_rad,
                                 &k,
                                &vampl).unwrap()
                                + project_vrot(pulse_config.star_config.rotation_velocity,
                                    theta_rad,
                                    phi_rad, 
                                    &k).unwrap();
                        }
                    }
                    local_temp = local_value(&pulse_config.temperature_config,
                        &pulse_config.mode_config,
                        theta_rad,
                        phi_rad,
                        pulse_config.star_config.effective_temperature).unwrap();
                    local_logg= local_value(&pulse_config.gravity_config,
                        &pulse_config.mode_config,
                        theta_rad,
                        phi_rad,
                        log_g0).unwrap();
                    local_area = project_vector(&s_normal,
                        &k_spherical).unwrap();

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
                if pulse_config.print_amplitude {
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
                let coschiout = |cos_chi: f64|{if cos_chi>0.0 {cos_chi} else { 0.0}};
                collect_output(&mut Star_Output,
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


        
    }//end for time loop
    match write_output_to_parquet(&pulse_config, Star_Output){
        Ok(_)=>{println!("\n Returned nice and well to main function")}
        Err(e)=>{println!("Couldnt create data frame and panicqued in main function with error {}",e)}
    };

    print_report(&now,//<---This gives the time of the computation
            &pulse_config,
            &vampl,
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
