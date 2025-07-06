use std::time::{Instant};
use crate::{utils::PulstarConfig, };
use temp_name_lib::{
    joris_math::spherical_harmonics::norm_factor::ylmnorm,
    type_def::{THETA_STEP,PHI_STEP}
    };
pub fn print_report(now:&Instant,
    parameters: &PulstarConfig,
    vampl:&[f64],
    k_theory:&[f64],
    min_vel:&[f64],
    max_vel:&[f64],
    min_t:&[f64],
    max_t:&[f64],
    min_logg:&[f64],
    max_logg:&[f64],
    freqrad:&[f64],
    period:&[f64],
    maxvellength:f64,
    maxrellength:f64,
    log_g0:f64){
    println!("+---------------------------------------------+");
    println!("PULSTARRust: REPORT - TIME: {} s", now.elapsed().as_secs());
    println!("+---------------------------------------------+\n");
    


    println!("PULSATION PARAMETERS:");
    println!("+---+-------+------------+--------------+------------+---------------+");
    println!("| # | (l,m) | freq (c/d) | freq (rad/s) | period (h) | ampl.(xi_r/r) |");
    println!("+---+-------+------------+--------------+------------+---------------+");

    for (index,l) in parameters.mode_config.l.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",l,parameters.mode_config.m[index]);
        print!("|  {:8.5}  ",parameters.freqcycli[index]);
        print!("|   {:8.5}   ",freqrad[index]);
        print!("|  {:8.5}  ",period[index]);
        print!("|    {:8.5}    ",parameters.mode_config.rel_deltar[index]);
    }
    println!("\n+---+-------+------------+--------------+------------+---------------+");
    println!(  "| # | (l,m) |  Vp (km/s) |   K (user)   | K (theory) | Y_l^m norm    |");
    println!(  "+---+-------+------------+--------------+------------+---------------+");
    for (index,l) in parameters.mode_config.l.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",l,parameters.mode_config.m[index]);
        print!("|  {:8.5}  ",vampl[index]);
        print!("|   {:8.5}   ",parameters.mode_config.k[index]);
        print!("|  {:8.5}  ",k_theory[index]);
        print!("|    {:8.5}    ", ylmnorm(*l,parameters.mode_config.m[index])); 
    }
    println!("\n+---+-------+------------+--------------+------------+---------------+");
    println!(  "| # | (l,m) | T_e factor |T_e phase dif |  g factor  | g phase dif   |");
    println!(  "+---+-------+------------+--------------+------------+---------------+");
    for (index,l) in parameters.mode_config.l.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",l,parameters.mode_config.m[index]);
        print!("|  {:8.3e}  ",(parameters.temperature_config[index].ampl).to_radians());
        print!("|   {:8.3e}   ",(parameters.temperature_config[index].phasedif).to_radians());
        print!("|  {:8.3e}  ",parameters.gravity_config[index].ampl);
        print!("|    {:8.3e}    ", (parameters.gravity_config[index].phasedif).to_radians()); 
    }

    println!("\n+---+-------+--------------+");
    println!(  "| # | (l,m) | phase offset |");
    println!(  "+---+-------+--------------+");
    for (index,l) in parameters.mode_config.l.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",l,parameters.mode_config.m[index]);
        print!("|   {:8.3e}  \n",parameters.mode_config.phase[index]);
    }

    print!("- Ve: {:8.5} km/s ",parameters.star_config.rotation_velocity);
    print!(" Vsini: {:8.5} km/s ",parameters.star_config.rotation_velocity 
        * ( (parameters.star_config.inclination_angle as f64).to_radians()).sin());
    println!(" Inclination angle: {} degrees", parameters.star_config.inclination_angle);

    println!("\nVISIBLE SURFACE DATA");
    println!("+-------+-----------------+-----------------+------------+-----------+--------------+-------------+");
    println!("| phase | Min. Proj. Vtot | Max. Proj. Vtot |  Min. T    |  Max. T   | Min. log(g)  | Max. log(g) |");
    println!("+-------+-----------------+-----------------+------------+-----------+--------------+-------------+");

    for (n,min_v) in min_vel.iter().enumerate(){
        print!("|  {:0>3}  ",n);
        print!("|     {:8.3e}    ",*min_v);
        print!("|     {:8.3e}    ",max_vel[n]);
        print!("|  {:8.3e} ",min_t[n]);
        print!("|  {:8.3e} ",max_t[n]);
        print!("|  {:8.3e}   ",min_logg[n]);
        print!("|  {:8.3e}   |\n",max_logg[n]);
    }

    let maxvel = max_vel.iter().fold(max_vel[0], |prev,curr| prev.max(*curr));
    let minvel = min_vel.iter().fold(min_vel[0], |prev,curr| prev.min(*curr));
    let maxt = max_t.iter().fold(max_t[0], |prev,curr| prev.max(*curr));
    let mint = min_t.iter().fold(min_t[0], |prev,curr| prev.min(*curr));
    let maxlogg = max_logg.iter().fold(max_logg[0], |prev,curr| prev.max(*curr));
    let minlogg = min_logg.iter().fold(min_logg[0], |prev,curr| prev.min(*curr));
    
    println!("EXTREMA:");
    println!("Total veloc.:  {:8.3} -> {:8.3} (diff. = {:8.3} )",minvel,maxvel,maxvel-minvel);
    println!("T_eff :  {:8.3} -> {:8.3} (diff. = {:8.3} )",mint,maxt,maxt-mint);
    println!("- log(g) :  {:8.3} -> {:8.3} (diff. = {:8.3} )",minlogg,maxlogg,maxlogg-minlogg);

    if parameters.print_amplitude {
        println!("\n AMPLITUDE GUARDING:");
        println!(" - Maximum length of puls. velocity vector: {:8.4}", maxvellength);
        println!(" - Maximum length of rel. displacement vecotr: {:8.4} \n", maxrellength);
    }

    println!("STAR PARAMETERS");
    println!(" - Mass/Mass_sun: {:8.4}",parameters.star_config.mass);
    println!(" - Radius/Radius_sun: {:8.4}", parameters.star_config.radius);
    println!(" - T_eff (Kelvin): {:8.4}", parameters.star_config.effective_temperature);
    println!(" - Log(g_0): {:8.4}\n",log_g0);

    println!("RESOLUTION");
    println!(" - Delta theta: {}",THETA_STEP);
    println!(" - Delta phi: {}\n",PHI_STEP);

    println!("COMPUTATIONAL RESTRICTIONS");
    println!(" - Surface normal time variation taken into account: ");
    match parameters.is_time_dependent{
        true => { println!( "   yes\n")}
        false => {println!("   no\n")}
    }
    println!(" -  Artificial suppression of pulsational velocity field: ");
    match parameters.suppress_pulse{
        true=> {println!("   yes\n")}
        false=> {println!("   no\n")}
    }
}