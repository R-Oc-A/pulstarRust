use std::{env,fs::File};
use polars::prelude::{LazyFrame,LazyGroupBy};
use polars::prelude::*;
use profile::intensity::IntensityGrids;
use temp_name_lib::type_def::{CLIGHT, N_FLUX_POINTS};
use profile::*;
use std::sync::Arc;
//use pulstar::utils::PulstarConfig;

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
    n_phases:10,
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

    let lambda_lower = wavelength[0] * (1.0 - profile_config.v_max/CLIGHT*1.0e3);
    let lambda_upper= wavelength.last().unwrap() * (1.0 + profile_config.v_max/CLIGHT*10.0e3);


    let mut flux = vec![0.0;capacity];
    let mut cont = vec![0.0;capacity];
    //open parquet file with pulsations
    //open txt file as csv !heel moeilijk

    //group by timepoints

    //Loop over the phases
//    for i in 0..=profile_config.n_phases{
//        let mut flux=vec![0.0;capacity];
//       let mut continuum = vec![0.0;capacity];
        


        //----------------------------------------
        // Start the loop over all points of the visible surface
        // Compute the local line profile in each point. 
        // Add the appropriate contribution to each point of the global line profile.
        
        // lambda_rest = wavelength[j]/(1.0_Vtotal/CLIGHT * 1.0e3)


        //write the profile lines!<- gold!
//    }
    let path= String::from("pulstar_output.parquet");
    let lf = LazyFrame::scan_parquet(path, Default::default()).unwrap();
    
    let grid_id_df = IntensityGrids{
        temperatures:vec![20000.0,20000.0,21000.0,21000.0],
        log_g: vec![3.5,3.5,4.0,4.0],
        filenames: vec![String::from("name1"),String::from("name2"),String::from("name3"),String::from("name4")]
    }.create_dataframe_sorted();
    
    let grid_id_lf= grid_id_df.lazy();
    //get vector of time_points
    let tf = lf.clone().select([col("time").unique(),]).collect().unwrap();
    let extract_time_series = tf.column("time").unwrap();
    let time_points:Vec<f64> = extract_time_series.f64().unwrap().into_iter().flatten().collect();
    let theta_df = lf.clone().select([col("theta").unique()]).collect().unwrap();
    let theta_steps_serie = theta_df.column("theta").unwrap();
    let theta_steps:Vec<f64> = theta_steps_serie.f64().unwrap().into_iter().flatten().collect();
    
    for phases in time_points.iter() {
        //get vector of theta points
        let expr = col("time").eq(lit(*phases));
        let theta_frame = lf.clone().filter(expr);

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

            for i in 0..coschi.len(){
                return_flux_for_cell_thetaphi(
                    coschi[i],
                    temperature[i],
                    log_gravity[i],
                    doppler_shift[i],
                    area[i],
                    &wavelength,
                    grid_id_lf.clone()
                );
            }
        }

        //write profile line
    }

    /*let mut output_schema= lf.collect_schema().unwrap();
    
    //output_schema.remove("theta".into());
    //output_schema.remove("phi".into());
    //output_schema.remove("time".into());

    output_schema.insert("flux".into(), DataType::List(Box::new(DataType::Float64)));
    output_schema.insert("continuum".into(), DataType::List(Box::new(DataType::Float64)));
    
    let analyzed = lf
    .filter(col("coschi").gt(lit(0.0)))
    .group_by(["time","theta","phi"])
    .apply(|group_df| {
        let doppler_shift = group_df.column("doppler shift")?.get(0)?.try_extract::<f64>()?;
        let temperature = group_df.column("temperature")?.get(0)?.try_extract::<f64>()?;
        let log_gravity = group_df.column("log gravity")?.get(0)?.try_extract::<f64>()?;
        let coschi = group_df.column("coschi")?.get(0)?.try_extract::<f64>()?;
        let area = group_df.column("area")?.get(0)?.try_extract::<f64>()?;
        
        let flux_vec = vec![0.0;capacity];
        let cont_vec = vec![0.0;capacity];
        
        let analysis_series = Series::new("flux".into(), flux_vec);
        let analysis_series2 = Series::new("continuum".into(),cont_vec);
        
        let res = group_df.with_columns(analysis_series);




    },output_schema)  
    .collect();
    */

    
}
