use std::path::PathBuf;
use polars:: prelude::*;


/// This function writes all of the data calculated into a single parquet. It may be modified into a version that handles memory more efficiently as
/// this one currently stores all of the calculations into big arrays. 
/// 
/// ### Arguments: 
/// 
/// `output_file_name` - a borrowed string slice that indicates the name of the output file
/// 
/// `wavelength_vec` - a borrowed vector that contains all of the observed wavelengths, these are repeated for every time point
/// 
/// `all_flux` - a borrowed vector that contains all of the intensity fluxes computed for the observed wavelengths. These fluxes are normalized and dimensionless.
/// 
/// `all_cont` - a borrowed vector that contains all of the continuum fluxes computed for the observed wavelengths. These fluxes are normalized and dimensionless.
/// 
/// `all_time` - a borrowed vector that contains the time points repeated on a redundant way. This redundancy is handle by the parquet file.
/// 
/// ### Returns: 
/// 
/// This function doesn't return a value however, it creates an output_file.parquet that has the information stored as columns with the headers 
/// 
/// `|time|wave length|flux|continuum|mean flux|`
/// 
/// where `mean flux` is `flux / continuum`.
pub fn write_into_parquet(
    output_file_name: &str,
    wavelength_vec:&[f64],
    all_flux: &[f64],
    all_cont: &[f64],
    all_time: &[f64]
){
    let time_0 = all_time[0];
    // create a data frame with all of the data, then create a lazy frame
    let output_df = df!(
        "time" => all_time,
        "wave length" => wavelength_vec,
        "flux" => all_flux,
        "continuum" => all_cont,
    ).expect("something went very wrong");
    
    let lf = output_df.lazy();
    
    // construct the mean flux expresion for the lazy data frame flux/cont
    let expr = col("flux") / col("continuum").alias("mean flux");
    let output_lf = lf.with_column(expr);

    // write lazy frame into parquet
    let mut path = PathBuf::new();
    path.push(output_file_name);
    let stream_lf = output_lf
    .sink_parquet(SinkTarget::Path(Arc::new(path)),
    ParquetWriteOptions::default(), 
     None,
     SinkOptions::default()).unwrap();
    stream_lf.collect().unwrap();

    // print 5 rows of the parquet output
    let llf = LazyFrame::scan_parquet(output_file_name,
    ScanArgsParquet::default()).unwrap();

    println!("---------------------------------------");
    println!("DataFrame for  Profile output  opened ");
    println!("First 5 rows:");
    println!("{}", llf.filter(col("time").eq(lit(time_0))).collect()
    .unwrap().head(Some(5usize)));

    // wait for some one to ask "what's up doc?"
}

