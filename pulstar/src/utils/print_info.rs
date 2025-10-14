use std::{f64::consts::PI, time::Instant};
use crate::PulstarConfig;
use temp_name_lib::{
    math_module::spherical_harmonics::norm_factor::ylmnorm, type_def::CYCLI2RAD,
    };
use polars::prelude::*;
pub fn print_report(now:&Instant,
    parameters: &PulstarConfig,
    time_points:usize){
    
    
    let new_path = std::path::PathBuf::from(format!("rasterized_star_{}tp.parquet",time_points));
    let mut file = std::fs::File::open(new_path).unwrap();
    let ddf = ParquetReader::new(&mut file).finish().unwrap();

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
        print!("|   {:8.5}   ",mode.frequency*CYCLI2RAD);
        print!("|  {:8.5}  ",2.0*PI/3.6e3/(mode.frequency*CYCLI2RAD));
        print!("|    {:8.5}    ",mode.rel_dr);
    }
    println!("\n+---+-------+------------+--------------+------------+---------------+");
    println!(  "| # | (l,m) |  Vp (km/s) |   K (user)   | K (theory) | Y_l^m norm    |");
    println!(  "+---+-------+------------+--------------+------------+---------------+");
    for (index,mode) in parameters.mode_data.iter().enumerate(){
        let k_theory = 74.437 * parameters.star_data.mass 
                /parameters.star_data.radius.powi(3)//r^3
                /((mode.frequency).powi(2));//freq^2
        print!("| {} ",index+1);
        print!("| ({},{}) ",mode.l,mode.m);
        print!("|  {:8.5}  ",vampl[index]);
        print!("|   {:8.5}   ",mode.k);
        print!("|  {:8.5}  ",k_theory);
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
    println!("+-----------------+-----------------+------------+-----------+--------------+-------------+");
    println!("| Min. Proj. Vtot | Max. Proj. Vtot |  Min. T    |  Max. T   | Min. log(g)  | Max. log(g) |");
    println!("+-----------------+-----------------+------------+-----------+--------------+-------------+");
    let vels = extract_column_as_vectorf64("velocity", &ddf);
    let tefs = extract_column_as_vectorf64("temperature", &ddf);
    let loggs = extract_column_as_vectorf64("log gravity", &ddf);

    let max_vel = get_max(&vels);
    let min_vel = get_min(&vels);
    let max_t = get_max(&tefs);
    let min_t = get_min(&tefs);
    let max_logg = get_max(&loggs);
    let min_logg = get_min(&loggs);
    print!("|     {:8.3e}    ",min_vel);
    print!("|     {:8.3e}    ",max_vel);
    print!("|  {:8.3e} ",min_t);
    print!("|  {:8.3e} ",max_t);
    print!("|  {:8.3e}   ",min_logg);
    print!("|  {:8.3e}   |\n",max_logg);

    println!("EXTREMA:");
    println!("Total veloc.:  {:8.3} -> {:8.3} (diff. = {:8.3} )",min_vel,max_vel,max_vel-min_vel);
    println!("T_eff :  {:8.3} -> {:8.3} (diff. = {:8.3} )",min_t,max_t,max_t-min_t);
    println!("- log(g) :  {:8.3} -> {:8.3} (diff. = {:8.3} )",min_logg,max_logg,max_logg-min_logg);


    println!("STAR PARAMETERS");
    println!(" - Mass/Mass_sun: {:8.4}",parameters.star_data.mass);
    println!(" - Radius/Radius_sun: {:8.4}", parameters.star_data.radius);
    println!(" - T_eff (Kelvin): {:8.4}", parameters.star_data.effective_temperature);
    println!(" - Log(g_0): {:8.4}\n",(max_logg+min_logg)*0.5);//This is a very lazy way of computing surface gravity...

    println!("RESOLUTION");
    println!(" - Delta theta: {}",theta_step);
    println!(" - Delta phi: {}\n",phi_step);

}



/// This function takes a polars data frame and returns all of the values from a given column that holds f64 values. 
/// ### Arguments: 
/// * `column_name` - a string slice that holds the name of a column. The column should hold f64 values.
/// * `df`- a polars DataFrame
/// ### Returns:
/// * `Vec<f64>` - a vector that contains all of the values on the column.
fn extract_column_as_vectorf64(column_name: &str,df:&DataFrame)->Vec<f64>{
    let column = df.column(column_name).unwrap();
    column.f64().unwrap().into_iter().flatten().collect()
}

fn get_max(vec:&Vec<f64>)->f64{
    vec.iter().fold(0.0, |accumulator,item| 
    if *item>accumulator {*item} else{accumulator})
}

fn get_min(vec:&Vec<f64>)->f64{
    vec.iter().fold(0.0, |accumulator,item| 
    if *item<accumulator ||*item>0.0 {*item} else{accumulator})
    
}