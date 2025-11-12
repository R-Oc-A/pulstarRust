use super::*;

impl GridsData{
    /// This function extracts all of the values outside from the data frames into vectors. With this setup the number of allocations is greatly reduced. 
    /// ### Arguments: 
    /// -dfs: [IntensityDataFrames] The in-between strucutre that  contains all of the data from the grids to be used for interpolation. 
    /// ### Returns: 
    /// - An instance of [GridsData] that holds all of the relevant information of the dataframes. 
    pub fn construct_from_dataframes(dfs: IntensityDataFrames)->Self{
        //Initialize the limb darkening coefficients for flux and continuum. 
        let mut limb_coef_flux:Vec<Vec<Vec<f64>>>=Vec::new();
        let mut limb_coef_cont:Vec<Vec<Vec<f64>>>=Vec::new();
        // Initialize the vectors that will hold the flux and continuum specific intensities. 
        let mut flux:Vec<Vec<f64>> = Vec::new();
        let mut continuum: Vec<Vec<f64>> = Vec::new();
        // Fill the wavelength array. The grids that we'll be using are homogeneous with a spacing of  0.001 nm
        let grid_wavelengths = extract_column_as_vectorf64("wavelength", &dfs.intensity_dfs[0]);
        // Initialize the array that will hold the relevant indices for the interpolation procedure. This is used so that the query is done only once.  
        let indices = vec![0usize;grid_wavelengths.len()*2];
        // This array holds the 4 grids that  will be used to perform the interpolation.
        //[Ricardo] Variations of Log G and effective temperature shouldn't be so great as to load more than 4 grids, however I leave this possibility coded
        let grids_indices = vec![0usize;4];

        // Read each dataframe and extract the columns
        for dataframe in dfs.intensity_dfs.iter(){
            flux.push(vec![0.0;grid_wavelengths.len()]);
            continuum.push(vec![0.0;grid_wavelengths.len()]);
            let a = extract_column_as_vectorf64("a",dataframe);
            let b = extract_column_as_vectorf64("b",dataframe);
            let c = extract_column_as_vectorf64("c",dataframe);
            let d = extract_column_as_vectorf64("d",dataframe);
            let ac = extract_column_as_vectorf64("ac",dataframe);
            let bc = extract_column_as_vectorf64("bc",dataframe);
            let cc = extract_column_as_vectorf64("cc",dataframe);
            let dc = extract_column_as_vectorf64("dc",dataframe);

            let limb_flux = vec![a,b,c,d];
            let limb_cont = vec![ac,bc,cc,dc];
            
            // Fill the nested vector structures
            limb_coef_flux.push(limb_flux);
            limb_coef_cont.push(limb_cont);
        }

        //Construct the instance of GridsData
        GridsData{
            temperature_vector: dfs.temperature_vector,
            log_g_vector:dfs.log_g_vector,
            grid_wavelengths:grid_wavelengths,
            limb_coef_cont:limb_coef_cont,
            limb_coef_flux:limb_coef_flux,
            flux:flux,
            continuum:continuum,
            row_indices: indices,
            grids_indices:grids_indices,
        }
    }
}

impl GridsData{
    /// This method returns the four indices that construct the parameter space of where to find the intensity [DataFrame]s
    /// 
    /// ### Arguments: 
    /// * `temperature` - The temperature of a surface cell
    /// * `log_gravity` - The log gravity value of  a surface cell 
    /// ### Returns:
    /// * This function returns the four indices of the intensity data frames that will be used for interpolation
    /// 
    /// (upper_temp,lower log_g)------------(upper_temp,upper_logg)
    /// 
    /// ------------------(temperature,log_gravity)----------------
    /// 
    /// (lower_temp,lower log g)------------(lower_temp,upper_logg)
    pub fn select_grids(&mut self,
    temperature:f64,
    log_gravity:f64)->PolarsResult<()>{
        
        //temp_vec = self.temperature_vector.clone();
        // lower_temperature_val
        let l_temperature_val = self.temperature_vector.iter()
            .fold(0.0,
            |accumulator, temp| { if *temp <= temperature { *temp } else {accumulator}});
        // upper_temperature_val
        let u_temperature_val = self.temperature_vector.iter()
            .fold(0.0,
            |accumulator, temp| { 
                if *temp <= temperature{ *temp } 
                else if accumulator<=temperature && temperature <= *temp{*temp}
                else {accumulator}
                });
        

        self.grids_indices.fill(0); 
        let mut counter = 0;
        let mut not_found_yet_lb =true;
        let mut not_found_yet_ub =true;

        for (index,temperature_val) in self.temperature_vector.iter().enumerate(){
			// look for l_boundary
			if *temperature_val==l_temperature_val && not_found_yet_lb{
			    if let Some(value) = self.log_g_vector.get(index+1){
			        if self.log_g_vector[index]<=log_gravity && log_gravity<=*value && not_found_yet_lb{
			            self.grids_indices[counter]=index;
			            self.grids_indices[counter+1] = index+1;
                        counter+=1;
			            not_found_yet_lb = false;
			        }
			    }
			}
			// look for u_boundary
			if *temperature_val == u_temperature_val && not_found_yet_ub{
			    if let Some(value) = self.log_g_vector.get(index+1){
			        if self.log_g_vector[index]<=log_gravity && log_gravity<=*value && not_found_yet_ub{
                        self.grids_indices[counter] = index;
                        self.grids_indices[counter+1] = index+1;
                        counter +=2;
			            not_found_yet_ub = false;
			        }
			    }
			}
        };
        if counter==3{ Ok(())}
        else{
            Err(PolarsError::
                InvalidOperation(
                    ErrString::new_static("Unable to load grid files for a surface cell; temperature and log g out of range. Load more grids.")))
        }
    }


    
}


/// This function is used to create a polars LazyFrame out of the intensity grid file in order to perform 
/// column wise operations faster. 
/// The Kurukz intensity files contain the coefficients $a_k$ useful to calculate the intensity via the limb darkening law
/// given by 
/// $ I_{\lambda}(\mu) =a_0 \sum_{k=1}^{3} a_k \left( 1 - \mu^k  \right) $
/// where $\mu$ is the $cos(\xi)$
///  
/// ### Arguments:
/// * `path` - a reference to a string slice that contains the relative path to the intensity grid file
/// ### Returns:
/// * `PolarsResult<LazyFrame>` - where the lazy frame has as headers
/// * `|wavelength|a|b|c|d|ac|bc|cc|dc|`
// This function returns a lazy frame with the apropiate column names of the intensity grid file
// the column names are 
// "wavelength","a","b","c","d","ac","bc","cc","dc"
// Input: path of the intensity grid file as a str
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

/// This function constructs a polars expression [Expr] that filters out wavelengths that are greater than or less than
/// observed (requested) ones. This reduces memory consumption and computation time.
/// 
/// ### Arguments: 
/// * `shifted_wavelengths` - The observed wavelengths
/// * `threshold` - a f64 value that specifies how close wavelengths need to be.
/// ### Returns:
/// *`Option<Expr>` -  where `Expr` is polars expression  that filters out unrelevant wavelengths from the lazyframe of an intensity grid file.
/// 
pub fn filter1_if_contains_wavelenghts(
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64)->Option<Expr>{

    
    let epsilon = 0.01;

    let min_wavelength = minval_rel_dopplershift 
        * wavelengths.iter()
        .fold(wavelengths[0],
            |accumulator,wavelength_val| 
            if*wavelength_val < accumulator {*wavelength_val} else {accumulator}
        )-epsilon;
    let max_wavelength = maxval_rel_dopplershift 
        * wavelengths.iter()
        .fold(wavelengths[0],
            |accumulator,wavelength_val| 
            if*wavelength_val > accumulator {*wavelength_val} else {accumulator}
        )+epsilon;    

    let filter_lower_expr = col("wavelength").gt(lit(min_wavelength));
    let filter_greater_expr = col("wavelength").lt(lit(max_wavelength));

    let combined_filter_exp = filter_lower_expr.or(filter_greater_expr);
    Some(combined_filter_exp)
}

// There might be necessary to apply another filter of the kind 
// Include wavelength only if observed_wavelength*min_rel_dopplershift - epsilon < wavelength <observed_wavelength * max_rel_dopplershift + epsilon
pub fn filter2_sift_wavelengths(
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64)->Option<Expr> {

    let epsilon = 1.0e-3;    
    let mut combined_expresion: Option<Expr> = None;
    
    for wavelength in wavelengths.iter(){
        let lb_wavelength = col("wavelength").gt(lit(wavelength*minval_rel_dopplershift - epsilon));
        let ub_wavelength = col("wavelength").lt(lit(wavelength*maxval_rel_dopplershift + epsilon));
        
        let current_mask = lb_wavelength.and(ub_wavelength);
        combined_expresion = match combined_expresion{
            Some(expression) => {Some(expression.or(current_mask))}
            None => {Some(current_mask)}
        }
    }    
    combined_expresion
}



///  This function materializes all of the filtered intensity data frames that will be used throughout the full program. It 
///  also fills an instance of the [IntensityDataFrames] 
/// 
/// ### Arguments:
/// * `grids_db` - a [DataFrame] that contains the name of the relevant gridfiles and their associated temperature and log_g values.
/// * `wavelengths` - a &[f64] vector (slice) that contains the observed wavelengths
/// * `max_rel_dopplershift` - a [f64] value that contains the maximum relative dopplershift
/// * `min_rel_dopplershift` - a [f64] value that contians the minimum relative dopplershift
/// 
/// ### Returns:
/// * an Instance of [IntensityDataFrames] where the temperature, log_g are ordered and it contains the [DataFrame]s of the intensity grid files
/// 
pub fn parse_relevant_intensity_grids(
    grids_db: DataFrame,
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64,
)->IntensityDataFrames{
    
    let temperature_vector = extract_column_as_vectorf64("temperature", &grids_db);
    let log_g_vector = extract_column_as_vectorf64("log_gravity", &grids_db);
    let filenames = extract_column_as_vector_string("file name", &grids_db);

    let is_low_resolution = (wavelengths[1]-wavelengths[0])>5.0e-3;//if the requested wavelengths are separated by more than 0.005 nm
    let mut intensity_dfs:Vec<DataFrame> = Vec::new();
    for name in filenames{
        let lf = read_intensity_grid_file(&name).expect("Unable to read intensity grid file");
        println!("loading {} grid file",name);

        let temp_df =lf.filter(
            filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift).unwrap()
            ).collect().expect("unable to produce dataframe using the intensity grid files");
        if is_low_resolution {
            intensity_dfs.push(
                sifted_wavelengths_dataframe(wavelengths,
                    temp_df.lazy(),
                    maxval_rel_dopplershift,
                    minval_rel_dopplershift).expect("unable to produce dataframe using the intensity grid files")
                );
        }
        else {
            intensity_dfs.push(temp_df);
        }
    }
    IntensityDataFrames{
        temperature_vector:temperature_vector,
        log_g_vector: log_g_vector,
        intensity_dfs:intensity_dfs,
    }
}


/// This function returns an intensity dataframe with only relevant wavelengths
/// The main purpose of this function is when you have a  low_resolution wavelength spectra
/// with bigger number of waves than 500. 
/// 
/// [polars] has a difficult time dealing with the combined expresion of more than 500 conditionals so the strategy implemented here
/// is to sepparate an intensity dataframe into chunks and sift the relevant wavelengths by chunks. 
/// 
/// ### Arguments:
/// * `wavelengths` - a [Vec<f64>] collection that contains the requested wavelengths.
/// * `intensity_df` - An intensity [DataFrame] about to be filtered
/// ### Returns: 
/// This function returns a [PolarsResult] with the following variants
/// * `Ok(DataFrame)` -where the binded data frame contains the sifted wavelengths
fn sifted_wavelengths_dataframe(
    wavelengths:&[f64],
    temp_lf: LazyFrame,
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64
)->PolarsResult<DataFrame>{
    //separate wavelengths into chunks of 500
    let chunk_size = 500usize;

    //initiallize collecting dataframe
    let mut collecting_df = df!("empty"=>[0.0]).unwrap();

    //loop over chunks
    for (iteration,chunk_of_wavelengths) in wavelengths.chunks(chunk_size).enumerate(){

        let partially_sifted_lf = temp_lf.clone().filter(filter2_sift_wavelengths(
            chunk_of_wavelengths,
            maxval_rel_dopplershift,
            minval_rel_dopplershift).unwrap());

        let new_df = 
        if iteration == 0{
            partially_sifted_lf.collect()?
        }
        else{
            let old_lf = collecting_df.lazy();
            concat([old_lf,partially_sifted_lf],
            UnionArgs::default())?.collect()?
        };

        collecting_df = new_df;        
    }
    Ok(collecting_df)
    //initialize the collecting dataframe
    //There's an iterator function called .chunks(usize) I shall look into it
    //loop over chunks
        //create data frame with the filetered wavelenghts from chunk
        //append with collecting data frame
            //create a lazyframe with thefiltered wavelengths dataframe 
            // create a lazyframe witht the collecting dataframe
    //
    //
}
