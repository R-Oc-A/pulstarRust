use temp_name_lib::{joris_math::{geometry::{projections::{cos_chi, project_vector}, surface_normal::surface_normal}, ref_frame_convrs::{cartesian_to_spherical, unit_vector_k}}, star_physics::{local_values::local_value, velocity::velocity_projection::{project_vpuls, project_vrot}}};
use pulstar::{utils::{self,parse_file,PulstarConfig}, DEG2RAD, PI, RADIUSSUN};
use nalgebra as na;
use core::time;
use std::{env};
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
                }//end if coschi > 0
                else{
                    //[Ricardo]: already nullified local values, so do nothing...
                }//end else coschi > 0   
                phi += 1;
            }//end phi loop
            theta += 1;
        }//end theta loop
    }//end for loop time points
    /*
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
