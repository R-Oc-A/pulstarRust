use std::env;
use polars::prelude::*;
use profile::intensity::IntensityGrids;
use temp_name_lib::type_def::N_FLUX_POINTS;
use profile::*;

fn main() {

   let env_args:Vec<String> = env::args().collect(); 

   if env_args.len() < 2usize {
    println!("Not enough arguments.");
    panic!("USAGE: profile -- <profile_input.txt>");
   }

   let profile_config=Config{
    lambda_0:412.0,
    lambda_f:413.0,
    delta_lbd:0.01,
    v_max:1.0e2,
   };


   let capacity = ((profile_config.lambda_f - profile_config.lambda_0)
        /profile_config.delta_lbd).floor() as usize + 1usize ;

    if capacity >= N_FLUX_POINTS as usize {panic!("Error, too many flux points requested.")}

    let mut wavelength:Vec<f64> = Vec::with_capacity(capacity);
    
    wavelength.push(profile_config.lambda_0);

    for i in 0..=capacity {//<- inclussive loop so wavelength[capacity]==λ_f.
            wavelength.push(profile_config.lambda_0 
                + profile_config.delta_lbd * i as f64 );
    }

    //let lambda_lower = wavelength[0] * (1.0 - profile_config.v_max/CLIGHT*1.0e3);
    //let lambda_upper= wavelength.last().unwrap() * (1.0 + profile_config.v_max/CLIGHT*10.0e3);

    let path= String::from("pulstar_output.parquet");
    let lf = LazyFrame::scan_parquet(path, Default::default()).unwrap();
    
    //Create a data frame with all of the intensity grids information
    let grid_id_df = IntensityGrids{
        temperatures:vec![20000.0,20000.0,21000.0,21000.0],
        log_g: vec![3.5,3.5,4.0,4.0],
        filenames: vec![String::from("name1"),String::from("name2"),String::from("name3"),String::from("name4")]
    }.create_dataframe_sorted();

    // Create a lazy frame that will be cloned inside the loop so that querys can be chained from here.
    let grid_id_lf_original = grid_id_df.lazy();
    
    //get vector of time_points
    let tf = lf.clone().select([col("time").unique(),]).collect().unwrap();
    let extract_time_series = tf.column("time").unwrap();
    let time_points:Vec<f64> = extract_time_series.f64().unwrap().into_iter().flatten().collect();
    //get vector of theta_points
    let theta_df = lf.clone().select([col("theta").unique()]).collect().unwrap();
    let theta_steps_serie = theta_df.column("theta").unwrap();
    let theta_steps:Vec<f64> = theta_steps_serie.f64().unwrap().into_iter().flatten().collect();

    let mut all_fluxes:Vec<f64> = Vec::new();
    let mut all_cont: Vec<f64> = Vec::new();
    let mut all_times: Vec<f64> = Vec::new();
    let mut all_wavelengths: Vec<f64> = Vec::new();
    //time loop    
    for phases in time_points.iter() {
        
        let expr = col("time").eq(lit(*phases));
        let theta_frame = lf.clone().filter(expr);

        let mut flux = vec![0.0;capacity];
        let mut cont = vec![0.0;capacity];
        let mut time_vec = vec![*phases;capacity];
        //theta loop
        for theta_step in theta_steps.iter(){
            
            let expr: Expr = col("theta").eq(lit(*theta_step));
            let phi_frame=theta_frame.clone().filter(expr);
            let expr = col("coschi").gt(lit(0.0));
            let pf_visible =phi_frame.filter(expr);

            let pf = append_doppler_shift(pf_visible).collect().unwrap();

            let doppler_shift_series = pf.column("doppler shift").unwrap();
            let area_series = pf.column("area").unwrap();
            let coschi_series = pf.column("coschi").unwrap();
            let temperature_series = pf.column("temperature").unwrap();
            let log_gravity_series = pf.column("log gravity").unwrap();

            let doppler_shift:Vec<f64> = doppler_shift_series.f64().unwrap().into_iter().flatten().collect();
            let area:Vec<f64> = area_series.f64().unwrap().into_iter().flatten().collect();
            let coschi:Vec<f64> = coschi_series.f64().unwrap().into_iter().flatten().collect();
            let temperature:Vec<f64> = temperature_series.f64().unwrap().into_iter().flatten().collect();
            let log_gravity:Vec<f64> = log_gravity_series.f64().unwrap().into_iter().flatten().collect();
            
            // phi loop
            for i in 0..coschi.len(){
                let grid_id_lf = grid_id_lf_original.clone();
                let fluxcont =return_flux_for_cell_thetaphi(
                    coschi[i],
                    temperature[i],
                    log_gravity[i],
                    doppler_shift[i],
                    area[i],
                    &wavelength,
                    grid_id_lf);
                let flux_thetaphi = fluxcont.0;
                let cont_thetaphi = fluxcont.1;
                //collect fluxex
                for (index,flux_item) in flux_thetaphi.into_iter().enumerate(){
                    flux[index] += flux_item;
                    cont[index] += cont_thetaphi[index]; 
                }//end collect fluxes

            }//end phi loop
        }//end theta loop
        all_fluxes.append(& mut flux);
        all_cont.append(& mut cont);
        all_times.append(& mut time_vec);
        all_wavelengths.extend(wavelength.iter().copied());
    }//end time loop
    
    //Write parquet file
    let output_file_name= String::from("wavelengths.parquet");
    utils::write_into_parquet(&output_file_name,
         &all_wavelengths, 
         &all_fluxes, 
         &all_cont, 
         &all_times);

    
}
