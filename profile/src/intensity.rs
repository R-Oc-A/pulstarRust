pub struct IntensityGrids{
    temperatures:Vec<f64>,
    log_g:Vec<f64>,
    filenames:Vec<String>,
}

impl IntensityGrids{
    pub fn create_dataframe_sorted(self)->DataFrame{
        let filenames_str:Vec<&str> = self.filenames.iter().map(|s| s.as_str()).collect();

        let temperature_series = Series::new("temperature".into(),self.temperatures);
        let log_gravity_serie=Series::new("log_gravity".into(),self.log_g);
        let file_name_serie = Series::new("file name".into(),filenames_str);
        
        let df = DataFrame::new(vec![temperature_series.into(),
            log_gravity_serie.into(),
            file_name_serie.into()]).unwrap()
            .sort(["temperature","log_gravity"], 
            SortMultipleOptions::default()
            .with_order_descending_multi([false,false]));

       df.unwrap()
    }

    pub fn sort_intensitygrid(self)->IntensityGrids{
        let df = self.create_dataframe_sorted();
        let temperature_series = df.column("temperature").unwrap();
        let log_gravity_series = df.column("log_gravity").unwrap();
        let file_name_series = df.column("file name").unwrap();

        let temperatures:Vec<f64> = temperature_series.f64().unwrap().into_iter().flatten().collect();
        let log_g:Vec<f64> = log_gravity_series.f64().unwrap().into_iter().flatten().collect();
        let filenames_str:Vec<&str> = file_name_series.str().unwrap().into_iter().flatten().collect();

        let filenames:Vec<String> = filenames_str.iter().map(|s| s.to_string()).collect();

        IntensityGrids { temperatures, log_g, filenames}
    }
}

pub fn read_intensity_grid_file(path:&str) -> PolarsResult<LazyFrame> {
    //let file = File::open(&path)?;

    let schema =  Schema::from_iter(vec![
		Field::new("wavelength".into(), DataType::Float64),
		Field::new("a".into(), DataType::Float64),
		Field::new("b".into(), DataType::Float64),
		Field::new("c".into(), DataType::Float64),
		Field::new("d".into(), DataType::Float64),
		Field::new("ac".into(), DataType::Float64),
		Field::new("bc".into(), DataType::Float64),
		Field::new("cc".into(), DataType::Float64),
		Field::new("dc".into(), DataType::Float64)
    ]); 
   
    let lf= LazyCsvReader::new(path)
    .with_separator(b' ')
    .with_has_header(false)
    .with_schema(Some(Arc::new(schema))).finish()?;
    
    Ok(lf)

}