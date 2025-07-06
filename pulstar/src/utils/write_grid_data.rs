use polars::prelude::*;
use std::fs::File;
use std::io;

use crate::utils::PulstarConfig;

pub struct RasterStarOutput{//[Ricardo]: This is redundant, but parquet handles value repetition. 
    all_thetas: Vec<f64>,
    all_phis: Vec<f64>,
    all_times: Vec<f64>,
    all_vel: Vec<f64>,
    all_temp: Vec<f64>,
    all_logg: Vec<f64>,
    all_coschi: Vec<f64>,
    all_area: Vec<f64>,
}

impl RasterStarOutput{
    pub fn new_sphere(
        num_theta_steps:usize,
        num_phi_steps:usize,
        num_timepoints:usize,
    )->RasterStarOutput{
        let total_rows = num_theta_steps* num_phi_steps * num_timepoints;
        RasterStarOutput { 
            all_thetas: Vec::with_capacity(total_rows),
            all_phis: Vec::with_capacity(total_rows), 
            all_times: Vec::with_capacity(total_rows),
            all_vel: Vec::with_capacity(total_rows),
            all_temp: Vec::with_capacity(total_rows),
            all_logg: Vec::with_capacity(total_rows),
            all_coschi: Vec::with_capacity(total_rows),
            all_area: Vec::with_capacity(total_rows) }
    }
}

pub fn collect_output( RasterStar: & mut RasterStarOutput,
        theta_rad:f64,
        phi_rad:f64,
        time:f64,
        local_veloc: f64,
        local_temp: f64,
        local_logg:f64,
        coschi:f64,
        local_area:f64,){

            RasterStar.all_thetas.push(theta_rad);
            RasterStar.all_phis.push(phi_rad);
            RasterStar.all_times.push(time);
            RasterStar.all_vel.push(local_veloc);
            RasterStar.all_temp.push(local_temp);
            RasterStar.all_logg.push(local_logg);
            RasterStar.all_coschi.push(coschi);
            RasterStar.all_area.push(local_area);
        
}


pub fn write_output_to_parquet(parameters: &PulstarConfig, output: RasterStarOutput) -> PolarsResult<()>{
    let mut df = DataFrame::new(vec![
        Series::new("theta".into(), output.all_thetas).into(),
        Series::new("phi".into(), output.all_phis).into(),
        Series::new("time".into(), output.all_times).into(),
        Series::new("velocity".into(), output.all_vel).into(),
        Series::new("temperature".into(), output.all_temp).into(),
        Series::new("log gravity".into(), output.all_logg).into(),
        Series::new("coschi".into(), output.all_coschi).into(),
        Series::new("area".into(), output.all_area).into(), 
    ])?;

    //println!("DataFrame for Raster Star Output  created ");
    //println!("First 5 rows:");
    //println!("{}", df.head(Some(5usize)));
    let id = String::from(parameters.mode_config.n_modes.to_string());
    let filename= format!("raster_star_tempid{}.parquet",id);
    
    let file_parquet = File::create(&filename)?;

    ParquetWriter::new(file_parquet)
        .finish(& mut df)?;
    
    println!("Parquet file for Raster Star output succesfully created");
    
    //Reading test
    let mut file = std::fs::File::open(&filename).unwrap();

    let ddf = ParquetReader::new(&mut file).finish().unwrap();
    println!("---------------------------------------");
    println!("DataFrame for Raster Star Output  opened ");
    println!("First 5 rows:");
    println!("{}", ddf.head(Some(5usize)));

    Ok(())
}


