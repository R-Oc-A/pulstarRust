use std::path::PathBuf;
use polars:: prelude::*;
use crate::FluxOfSpectra;
pub mod csv_to_ndarray;



pub struct IntensityFlux{
    pub data_frame:DataFrame,
}

impl IntensityFlux{
    pub fn new()->Self{
        Self{data_frame: df!(
            "time"=>Vec::<f64>::new(),
            "wave length" => Vec::<f64>::new(),
            "flux" => Vec::<f64>::new(),
            "continuum" => Vec::<f64>::new(),
            "normalized flux" => Vec::<f64>::new()
        ).unwrap(),}
    }

    pub fn append_fluxes(self,fluxes:FluxOfSpectra)->Self{
        let flux_df=df!(
            "time" => fluxes.time,
            "wave length" => fluxes.wavelengths,
            "flux" => fluxes.flux,
            "continuum" => fluxes.continuum
        ).unwrap();

        // construct the mean flux expresion for the lazy data frame flux/cont
        let expr = (col("flux") / col("continuum")).alias("normalized flux");
        let flux_lf=flux_df.lazy().with_column(expr);
        let result_lf=append_current_lf_into_collection_lf(flux_lf, self.data_frame.lazy()).unwrap();

        IntensityFlux { data_frame:result_lf.collect().unwrap()}
    }

    pub fn write_output(self, time_points:u16)->PolarsResult<()>{
        // write lazy frame into parquet
        let new_path = PathBuf::from(
            format!("wavelengths_tp{}.parquet",time_points)
        );
        let lf_to_write = self.data_frame.lazy();

        if let Ok(lf) = lf_to_write.sink_parquet(
            SinkTarget::Path(Arc::new(new_path.clone())),
            ParquetWriteOptions::default(),
            None,
            SinkOptions::default()){
                lf.collect()?;
            }else {eprint!("unable to sink to a parket in {} time_point",time_points)};
    
        // print 5 rows of the parquet output
        let llf = LazyFrame::scan_parquet(new_path,
        ScanArgsParquet::default()).unwrap();
    
        println!("---------------------------------------");
        println!("DataFrame for  Profile output  opened ");
        println!("First 5 rows:");
        println!("{}", llf.filter(col("time").eq(lit(time_points))).collect()
        .unwrap().head(Some(5usize)));
    
        Ok(())
        
    }

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


/// This function stacks the the LazyFrame of the spectrum into the collection stored in the parquet file and creates a lazyframe out of the handle.
/// ### Arguments: 
/// * `spectrum_lazyframe` - a [LazyFrame] created from the data frame of the rasterized star.
/// * `parquet_file_lazyframe` - a [LazyFrame] created from the parquet file
/// ### Returns:
///  This function returns a [PolarsResult] with the following variants:
/// * `Ok(LazyFrame)` - In case the [LazyFrame] was adequately created.
/// * `Err(PolarsError)` - Returning a [PolarsError] to the calling function. 
fn append_current_lf_into_collection_lf(spectra_lazyframe:LazyFrame,parquet_file_lazyframe:LazyFrame)->PolarsResult<LazyFrame>{
    concat(
        [parquet_file_lazyframe,spectra_lazyframe],
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
    let old_file =format!("wavelengths_tp{}.parquet",time_points-1); 
    println!("deletting {}",old_file);
    std::fs::remove_file(old_file)?;
    Ok(())
}


/// This function creates a [DataFrame] out of a [RasterizedStarOutput]. WARNING: This function takes ownership of the RasterizedStarOutput.
/// 
/// ### Arguments:
/// * `fluxes` - An instance of [FluxOfSpectra] that contains all of the intensity values for a specific time point.
/// ### Returns: 
///  This function returns a [PolarsResult] with the following variants:
/// * `Ok(DataFrame)` - In case the [DataFrame] was adequately created.
/// * `Err(PolarsError)` - Returning a [PolarsError] to the calling function. 
fn create_spectra_dataframe(fluxes: FluxOfSpectra)->PolarsResult<DataFrame>{
    // The df! macro creates a new dataframe with the columns ("column header"=>values) ordered from left to right
    df!(
        "time" => fluxes.time,
        "wave length" => fluxes.wavelengths,
        "flux" => fluxes.flux,
        "continuum" => fluxes.continuum
    )
}


///This function creates the [LazyFrame] that will be used to create the parquet file 
/// 
/// ### Arguments:
/// * `time_points` - a [u16] integer that indicates the time_point to be added. 
/// * `flux_lf` - the [LazyFrame] of the [FluxOfSpectra] [DataFrame]
/// ### Returns:
/// * [LazyFrame] - This lazyframe will be sinked ([polars::prelude::LazyFrame::sink_parquet]) into a parquet file
fn lazyframe_to_be_written (time_points:u16,flux_lf:LazyFrame)->PolarsResult<LazyFrame>{
    if time_points == 1{
        Ok(flux_lf)
    }else{
        let old_path = std::path::PathBuf::from(format!("wavelengths_tp{}.parquet",time_points-1));
        let old_lf = open_collecting_parquet_file_as_lazyframe(&old_path)?;
        Ok(append_current_lf_into_collection_lf(flux_lf, old_lf)?)
    }
}


/// This function writes all of the data calculated into a single parquet. It may be modified into a version that handles memory more efficiently as
/// this one currently stores all of the calculations into big arrays. 
/// 
/// ### Arguments: 
/// * `output_file_name` - a borrowed string slice that indicates the name of the output file
/// * `wavelength_vec` - a borrowed vector that contains all of the observed wavelengths, these are repeated for every time point
/// * `all_flux` - a borrowed vector that contains all of the intensity fluxes computed for the observed wavelengths. These fluxes are normalized and dimensionless.
/// * `all_cont` - a borrowed vector that contains all of the continuum fluxes computed for the observed wavelengths. These fluxes are normalized and dimensionless.
/// * `all_time` - a borrowed vector that contains the time points repeated on a redundant way. This redundancy is handle by the parquet file.
/// ### Returns: 
/// This function doesn't return a value however, it creates an output_file.parquet that has the information stored as columns with the headers 
/// * `|time|wave length|flux|continuum|mean flux|`
/// where `mean flux` is `flux / continuum`.
pub fn write_into_parquet(
    time_points:u16,
    fluxes: FluxOfSpectra,
)->PolarsResult<()>{
    //Initialize output_df
    let time_0 = fluxes.time[0];
    let output_df = create_spectra_dataframe(fluxes).expect("something went very wrong while creating output");
    let lf = output_df.lazy();
    
    // construct the mean flux expresion for the lazy data frame flux/cont
    let expr = (col("flux") / col("continuum")).alias("normalized flux");
    let flux_lf = lf.with_column(expr);//<--Here's an error, it should be something like select... let's see How I fix it. However it's weird That it functions when I just present the fluxes...

    // write lazy frame into parquet
    let new_path = PathBuf::from(
        format!("wavelengths_tp{}.parquet",time_points)
    );
    let lf_to_write = lazyframe_to_be_written(time_points,
         flux_lf.clone())?;

    if let Ok(lf) = lf_to_write.sink_parquet(
        SinkTarget::Path(Arc::new(new_path.clone())),
        ParquetWriteOptions::default(),
        None,
        SinkOptions::default()){
            lf.collect()?;
        }else {eprint!("unable to sink to a parket in {} time_point",time_points)};

    // print 5 rows of the parquet output
    let llf = LazyFrame::scan_parquet(new_path,
    ScanArgsParquet::default()).unwrap();

    println!("---------------------------------------");
    println!("DataFrame for  Profile output  opened ");
    println!("First 5 rows:");
    println!("{}", llf.filter(col("time").eq(lit(time_0))).collect()
    .unwrap().head(Some(5usize)));
    
    if time_points > 1u16
    {remove_temp_parquet_file(time_points)?};

    Ok(())
}

