use std::path::PathBuf;
use polars:: prelude::*;

pub fn write_into_parquet(
    output_file_name: &str,
    wavelength_vec:&[f64],
    all_flux: &[f64],
    all_cont: &[f64],
    all_time: &[f64]
){

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
    println!("{}", llf.collect()
    .unwrap().head(Some(5usize)));

    // wait for some one to ask "what's up doc?"
}

