use polars::prelude::*;
use polars::prelude::LazyFrame;

use crate::RasterizedStar;

/// This structure is an abstraction of the star. It effectively creates a columnar data base with all of the surface cells and the associated quantities.
/// 
pub struct RasterizedStarOutput{
    /// This is the column that indicates all of the colatitude angles θ_i
    all_coords1: Vec<f64>,
    /// This is the column that indicates all of the azimuthal angles φ_i
    all_coords2: Vec<f64>,
    /// This is the time point collumn. This value is redundant as it's the same for all of the star however this redundancy is necesary for the parquet file. 
    all_times: Vec<f64>,
    /// Collumn that collects all of the observed variations in total velocity with respect to the observer.
    all_vel: Vec<f64>,
    /// Collumn with all of the observed variations in temperature
    all_temp: Vec<f64>,
    /// Collumn with all of the observed variations of log_g
    all_logg: Vec<f64>,
    /// Collumn with all of the cosines of the angle between unit vector normal to a cell and a unit vector pointing towards the observer. This is necessary for calculating intensities with the limb darkening law.
    all_coschi: Vec<f64>,
    /// Collumn with all of the (observed)variations of the cell's area caused by pulsations and with respect to the observer. This quantity is normalized such that the sum of all areas is equal to 1. 
    all_area: Vec<f64>,
}

impl RasterizedStarOutput{
    /// This method creates a new instance of a rasterized star 
    pub fn new_sphere(
        num_theta_steps:usize,
        num_phi_steps:usize,
        num_timepoints:usize,
    )->RasterizedStarOutput{
        let total_rows = num_theta_steps* num_phi_steps * num_timepoints;
        RasterizedStarOutput { 
            all_coords1: Vec::with_capacity(total_rows),
            all_coords2: Vec::with_capacity(total_rows), 
            all_times: Vec::with_capacity(total_rows),
            all_vel: Vec::with_capacity(total_rows),
            all_temp: Vec::with_capacity(total_rows),
            all_logg: Vec::with_capacity(total_rows),
            all_coschi: Vec::with_capacity(total_rows),
            all_area: Vec::with_capacity(total_rows) }
    }

    pub fn format_the_star(
        star: &RasterizedStar,        
    )->RasterizedStarOutput{
        let capacity = star.cells.len();
        let mut all_coords1:Vec<f64> = Vec::with_capacity(capacity);
        let mut all_coords2:Vec<f64> = Vec::with_capacity(capacity);
        let mut all_times:Vec<f64> = Vec::with_capacity(capacity);
        let mut all_vel:Vec<f64> = Vec::with_capacity(capacity);
        let mut all_temp:Vec<f64> = Vec::with_capacity(capacity);
        let mut all_logg:Vec<f64> = Vec::with_capacity(capacity);
        let mut all_coschi:Vec<f64> = Vec::with_capacity(capacity);
        let mut all_area:Vec<f64> = Vec::with_capacity(capacity);


        for cell in star.cells.iter(){
            if cell.coschi!=0.0{
                all_coords1.push(cell.coord_1);
                all_coords2.push(cell.coord_2);
                all_times.push(star.time_stamp);
                all_vel.push(cell.v_tot);
                all_temp.push(cell.t_eff);
                all_logg.push(cell.log_g);
                all_coschi.push(cell.coschi);
                all_area.push(cell.area);
            }
        }
        RasterizedStarOutput{all_area:all_area,
            all_coords1:all_coords1,
            all_coords2:all_coords2,
            all_coschi:all_coschi,
            all_logg:all_logg,
            all_temp:all_temp,
            all_times:all_times,
            all_vel:all_vel}
    }
}

/// This function creates a [DataFrame] out of a [RasterizedStarOutput]. WARNING: This function takes ownership of the RasterizedStarOutput.
/// 
/// ### Arguments:
/// * `star` - An instance of a [RasterizedStarOutput] that contains all of the local values of all of the cell surfaces.
/// ### Returns: 
///  This function returns a [PolarsResult] with the following variants:
/// * `Ok(DataFrame)` - In case the [DataFrame] was adequately created.
/// * `Err(PolarsError)` - Returning a [PolarsError] to the calling function. 
fn create_rasterized_star_dataframe(star: RasterizedStarOutput)->PolarsResult<DataFrame>{
    // The df! macro creates a new dataframe with the columns ("column header"=>values) ordered from left to right
    df!(
        "coord1" => star.all_coords1,
        "coord2" => star.all_coords2,
        "time" => star.all_times,
        "velocity" => star.all_vel,
        "temperature" => star.all_temp,
        "log gravity" => star.all_logg,
        "coschi" => star.all_coschi,
        "area" => star.all_area,
    )
}

/// This function opens the parquet file and creates a lazyframe out of the handle.
/// ### Arguments: 
/// * `path to parquet` - a [std::path::PathBuf] that indicates the path and name to the parquet file.
/// ### Returns:
///  This function returns a [PolarsResult] with the following variants:
/// * `Ok(LazyFrame)` - In case the [LazyFrame] was adequately created.
/// * `Err(PolarsError)` - Returning a [PolarsError] to the calling function. 
fn open_collecting_parquet_file_as_lazyframe(path_to_parquet: &std::path::PathBuf)->PolarsResult<LazyFrame>{
    LazyFrame::scan_parquet(path_to_parquet, ScanArgsParquet::default())
}


/// This function stacks the the LazyFrame of the rasterized star into the collection stored in the parquet file and creates a lazyframe out of the handle.
/// ### Arguments: 
/// * `star_lazyframe` - a [LazyFrame] created from the data frame of the rasterized star.
/// * `parquet_file_lazyframe` - a [LazyFrame] created from the parquet file
/// ### Returns:
///  This function returns a [PolarsResult] with the following variants:
/// * `Ok(LazyFrame)` - In case the [LazyFrame] was adequately created.
/// * `Err(PolarsError)` - Returning a [PolarsError] to the calling function. 
fn append_current_lf_into_collection_lf(star_lazyframe:LazyFrame,parquet_file_lazyframe:LazyFrame)->PolarsResult<LazyFrame>{
    concat(
        [parquet_file_lazyframe,star_lazyframe],
        UnionArgs::default()
    )
}

/// This function removes the parquet file that holds the old collection of rasterized stars
/// 
/// ### Arguments: 
/// * `path_to_parquet` - A [std::path::PathBuf] path to the old parquet file. 
/// ### Returns: 
/// This function returns a [Result] with the following variants:
/// * `Ok(_)` - if everything went ok.
/// * `Err(std::io::Error)` - where the error is passed to the calling function to indicate that it could not remove the file. 
fn remove_temp_parquet_file(time_points:u16)->Result<(), std::io::Error>{
    //let old_path = std::path::PathBuf::from(format!("rasterized_star_{}tp.parquet",time_points-1));
    let old_file =format!("rasterized_star_{}tp.parquet",time_points-1); 
    println!("deletting {}",old_file);
    std::fs::remove_file(old_file)?;
    Ok(())
}


/// This function writes the output of the pulstar binary into a parquet file. 
/// This will collect all of the computations for each time point
/// 
/// ### Arguments: 
/// * `star_output` - a [RasterizedStarOutput] instance that contains all of the surface cells local values. 
/// * `time_points` - a [u16] integer that indicates the time_point to be added. 
/// * `old_parquet_file` - a [std::path::PathBuf] that contains the path to the old parquet file.
/// The parquet file will have a name with the convention "rasterized_star_<time_points>tp.parquet",
/// meaning that it contains also all of the previous time_points computations.
/// ### Returns:
/// This function returns a [PolarsResult] with the following variants:
/// * `Ok(_)` - if the file was succesfully created and written. 
/// * `Err(PolarsError)` - returns a [PolarsError] to the calling function in case that there was an error creating a data frame for the rasterized star,
///  the output file couldn't be created, or  the `old_output_file`` file wasn't succesfully erased. 
pub fn write_output_to_parquet(
    star: &RasterizedStar,
    time_points:u16,
    ) -> PolarsResult<()>{
    
    let star_output = RasterizedStarOutput::format_the_star(star);


    let new_path = std::path::PathBuf::from(format!("rasterized_star_{}tp.parquet",time_points));
    let star_df = create_rasterized_star_dataframe(star_output)?;
    let star_lf = star_df.lazy();
    
    let lf_to_write = lazyframe_to_be_written(time_points, star_lf.clone())?;
    if let Ok(lf) = lf_to_write.sink_parquet(
        SinkTarget::Path(Arc::new(new_path.clone())),
        ParquetWriteOptions::default(), 
        
        None, 
        SinkOptions::default()){
            lf.collect()?;
        }else {eprint!("unable to sink to a parket in {} time_point",time_points)};
    //Reading test
    println!("{:#?} file created",new_path);
    let mut file = std::fs::File::open(new_path).unwrap();
    let ddf = ParquetReader::new(&mut file).finish().unwrap();
    println!("---------------------------------------");
    println!("DataFrame for Rasterized Star Output  opened ");
    println!("First 5 rows:");
    println!("{}", ddf.head(Some(5usize)));
    if time_points > 1u16
        {remove_temp_parquet_file(time_points)?}
    
    Ok(())
    /*ParquetWriter::new(file_parquet)
        .finish(& mut df)?;
    
    println!("Parquet file for Rasterized Star output succesfully created");
    

    Ok(())*/
}

///This function creates the [LazyFrame] that will be used to create the parquet file 
/// 
/// ### Arguments:
/// * `time_points` - a [u16] integer that indicates the time_point to be added. 
/// * `star_lf` - the [LazyFrame] of the rasterized star [DataFrame]
/// ### Returns:
/// * [LazyFrame] - This lazyframe will be sinked ([polars::prelude::LazyFrame::sink_parquet]) into a parquet file
fn lazyframe_to_be_written (time_points:u16,star_lf:LazyFrame)->PolarsResult<LazyFrame>{
    if time_points == 1{
        Ok(star_lf)
    }else{
        let old_path = std::path::PathBuf::from(format!("rasterized_star_{}tp.parquet",time_points-1));
        let old_lf = open_collecting_parquet_file_as_lazyframe(&old_path)?;
        Ok(append_current_lf_into_collection_lf(star_lf, old_lf)?)
    }
}