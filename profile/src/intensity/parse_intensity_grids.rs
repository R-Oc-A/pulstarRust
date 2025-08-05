use super::*;
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
pub fn filter1_if_contains_wavelenghts(shifted_wavelengths:&[f64])->Option<Expr>{

    let max_wavelength = shifted_wavelengths.last().unwrap().clone();
    let min_wavelength = shifted_wavelengths.get(0).unwrap().clone();

    let filter_lower_expr = col("wavelength").gt(lit(min_wavelength));
    let filter_greater_expr = col("wavelength").lt(lit(max_wavelength));

    let combined_filter_exp = filter_lower_expr.or(filter_greater_expr);
    Some(combined_filter_exp)
}



/// This is a function that extracts only the relevant wavelengths as [Vec<f64>]. 
/// 
/// ### Arguments: 
/// * `shifted_wavelengts` - A reference to a [Vec<f64>] that contains the doppler shifted observed wavelengths.
/// * `grid_dataframe` - A [DataFrame] that contains the wavelengths from the grid file. This DataFrame has exluded all the wavelengths outside of range. 
/// ### Returns: 
/// * `relevant_df` - A [DataFrame] that only contains the wavelengths that'll be used for the [interpolate] step.
pub fn extract_relevant_wavelengths(
    shifted_wavelengths:&[f64],
    grid_dataframe:DataFrame)->DataFrame{
    
    let df = &grid_dataframe;
    // Extract df's collumns into vectors; Wavelengths are ordered from lower to greater so we'll profit from that.
        let df_wavelengths=extract_column_as_vectorf64("wavelength", &grid_dataframe);
        let df_a = extract_column_as_vectorf64("a", &grid_dataframe);
        let df_b = extract_column_as_vectorf64("b", &grid_dataframe);
        let df_c = extract_column_as_vectorf64("c", &grid_dataframe);
        let df_d = extract_column_as_vectorf64("d", &grid_dataframe);
        let df_ac = extract_column_as_vectorf64("ac", &grid_dataframe);
        let df_bc = extract_column_as_vectorf64("bc", &grid_dataframe);
        let df_cc = extract_column_as_vectorf64("cc", &grid_dataframe);
        let df_dc = extract_column_as_vectorf64("dc", &grid_dataframe);
    // Get relevant indices of the df_wavelength vector into an index vector
        let indices = extract_important_indices(shifted_wavelengths, &df_wavelengths);
    // Construct new vectors using only the relevant indices
        let wavelengths = construct_new_vec_using_indices(df_wavelengths, &indices);
        let a = construct_new_vec_using_indices(df_a, &indices);
        let b = construct_new_vec_using_indices(df_b, &indices);
        let c = construct_new_vec_using_indices(df_c, &indices);
        let d = construct_new_vec_using_indices(df_d, &indices);
        let ac = construct_new_vec_using_indices(df_ac, &indices);
        let bc = construct_new_vec_using_indices(df_bc, &indices);
        let cc = construct_new_vec_using_indices(df_cc, &indices);
        let dc = construct_new_vec_using_indices(df_dc, &indices);

    // Construct a data frame using the relevant indices
    let relevant_df= df![
        "wavelength"=>wavelengths,
        "a" => a,
        "b" => b,
        "c" => c,
        "d" => d,
        "ac" => ac,
        "bc" => bc,
        "cc" => cc,
        "dc" => dc
    ].unwrap();
    relevant_df
}

/// This function is useful to get the relevant indices out of the wavelength vector. 
///
/// ### Arguments:
/// * `shifted_wavelengts` - A reference to a [Vec<f64>] that contains the doppler shifted observed wavelengths.
/// * `df_wavelengts` - A reference to a [Vec<f64>] that contains the grid's wavelengths extracted from the grid [DataFrame].
/// ### Returns:
/// * `indices` - A [Vec<usize>] that contains the indices of the df_wavelengths that'll be used for the [interpolate] step.
fn extract_important_indices(
    shifted_wavelengths:&[f64],
    df_wavelengths: &[f64]
)->Vec<usize>{
    let mut indices:Vec<usize> = Vec::new();
    for wavelength in shifted_wavelengths.iter(){
        let index = search_geq(df_wavelengths, *wavelength);
        //extract last index
        /*if let Some(last_index)= indices.last(){
            if df_wavelengths[*last_index]!=df_wavelengths[index-1]{
                indices.push(index-1)
            }
        }else{
            indices.push(index-1);
        }*/
        indices.push(index-1);
        indices.push(index);
    }
    indices
}

/// This function filters a vector using a collection of indices
/// 
/// ### Arguments: 
/// * `original_vector` - A [Vec<f64>] that contains the unfiltered data from the grid files
/// * `index_collection` - A [ &[usize] ] that contains the relevant indices.
/// ### Returns:
/// * `filtered_vector` - A [Vec<f64>] that contains the filtered data from the grid files. 
fn construct_new_vec_using_indices(
    original_vector:Vec<f64>,
    index_collection:&[usize]
)->Vec<f64>{
    let mut filtered_vector:Vec<f64>= Vec::new();
    for index in index_collection.iter(){
        filtered_vector.push(original_vector[*index]);
    }
    filtered_vector
}

