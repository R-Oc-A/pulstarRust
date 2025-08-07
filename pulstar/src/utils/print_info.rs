use std::time::{Instant};
use crate::PulstarConfig;
use temp_name_lib::{
    math_module::spherical_harmonics::norm_factor::ylmnorm,
    };
pub fn print_report(now:&Instant,
    parameters: &PulstarConfig,
    k_theory:&[f64],
    min_vel:&[f64],
    max_vel:&[f64],
    min_t:&[f64],
    max_t:&[f64],
    min_logg:&[f64],
    max_logg:&[f64],
    freqrad:&[f64],
    period:&[f64],
    log_g0:f64){
    let mesh = parameters.get_mesh_structure();
    let theta_step = mesh.0;
    let phi_step = mesh.1;
    let vampl = parameters.get_velocity_amplitudes();
    println!("+---------------------------------------------+");
    println!("PULSTARRust: REPORT - TIME: {} s", now.elapsed().as_secs());
    println!("+---------------------------------------------+\n");
    


    println!("PULSATION PARAMETERS:");
    println!("+---+-------+------------+--------------+------------+---------------+");
    println!("| # | (l,m) | freq (c/d) | freq (rad/s) | period (h) | ampl.(xi_r/r) |");
    println!("+---+-------+------------+--------------+------------+---------------+");

    for (index,mode) in parameters.mode_data.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",mode.l,mode.m);
        print!("|  {:8.5}  ",mode.frequency);
        print!("|   {:8.5}   ",freqrad[index]);
        print!("|  {:8.5}  ",period[index]);
        print!("|    {:8.5}    ",mode.rel_dr);
    }
    println!("\n+---+-------+------------+--------------+------------+---------------+");
    println!(  "| # | (l,m) |  Vp (km/s) |   K (user)   | K (theory) | Y_l^m norm    |");
    println!(  "+---+-------+------------+--------------+------------+---------------+");
    for (index,mode) in parameters.mode_data.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",mode.l,mode.m);
        print!("|  {:8.5}  ",vampl[index]);
        print!("|   {:8.5}   ",mode.k);
        print!("|  {:8.5}  ",k_theory[index]);
        print!("|    {:8.5}    ", ylmnorm(mode.l,mode.m)); 
    }
    println!("\n+---+-------+------------+--------------+------------+---------------+");
    println!(  "| # | (l,m) | T_e factor |T_e phase dif |  g factor  | g phase dif   |");
    println!(  "+---+-------+------------+--------------+------------+---------------+");
    for (index,mode) in parameters.mode_data.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",mode.l,mode.m);
        print!("|  {:8.3e}  ",mode.rel_dtemp);
        print!("|   {:8.3e}   ",mode.phase_rel_dtemp.to_radians());
        print!("|  {:8.3e}  ",mode.rel_dg);
        print!("|    {:8.3e}    ", mode.phase_rel_dg.to_radians()); 
    }

    println!("\n+---+-------+--------------+");
    println!(  "| # | (l,m) | phase offset |");
    println!(  "+---+-------+--------------+");
    for (index,mode) in parameters.mode_data.iter().enumerate(){
        print!("| {} ",index+1);
        print!("| ({},{}) ",mode.l,mode.m);
        print!("|   {:8.3e}  \n",mode.phase_offset);
    }

    print!("- Ve: {:8.5} km/s ",parameters.star_data.v_omega);
    print!(" Vsini: {:8.5} km/s ",parameters.star_data.v_omega * parameters.star_data.inclination_angle.to_radians().sin());
    println!(" Inclination angle: {} degrees", parameters.star_data.inclination_angle);

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


    println!("STAR PARAMETERS");
    println!(" - Mass/Mass_sun: {:8.4}",parameters.star_data.mass);
    println!(" - Radius/Radius_sun: {:8.4}", parameters.star_data.radius);
    println!(" - T_eff (Kelvin): {:8.4}", parameters.star_data.effective_temperature);
    println!(" - Log(g_0): {:8.4}\n",log_g0);

    println!("RESOLUTION");
    println!(" - Delta theta: {}",theta_step);
    println!(" - Delta phi: {}\n",phi_step);

}