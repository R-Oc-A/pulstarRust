use ndarray::Data;

use super::*;

pub mod joris_grids;


fn add_padding_for_wavelength(wavelength:&[f64],
    maxval_rel_dopplershift:f64
)->Vec<f64>{

    let mut padded_wavelength_vec:Vec<f64> = Vec::new();
    let min_wavelength = if let Some( start) = wavelength.get(0){
        start * (1.0-maxval_rel_dopplershift)
    }else{
        panic!("What are you trying to do? Your wavelength array has length 0")
    };
    let max_wavelength = if let Some( end) = wavelength.last(){
        end * (1.0+maxval_rel_dopplershift)
    }else{
        panic!("What are you trying to do? Your wavelength array has length 0")
    };
    let last_wavelength = if let Some(value) = wavelength.last(){value}else{panic!()};
    
    let d_lambda = (wavelength[1]-wavelength[0]);
    
    let left_padding = {
        let mut index = 0usize;
        while (min_wavelength < wavelength[0]-f64::from(index)*d_lambda ){
            index += 1;
        }
        index
    };
    
    let right_padding = {
        let mut index = 0usize;
        while (max_wavelength > last_wavelength + f64::from(index)*d_lambda ){
        index += 1;
        }
        index
    };

    for index in (1..=left_padding).rev(){
        padded_wavelength_vec.push(wavelength[0] - f64::from(index)*d_lambda);
    }

    for lambda in wavelength.inter(){
        padded_wavelength_vec.push(lambda);
    }
    
    for index in 1..=right_padding{
        padded_wavelength_vec.push(last_wavelength + f64::from(index)*d_lambda);
    }

    padded_wavelength_vec

}


fn make_df_w0(wavelength:&[f64])->DataFrame{
    df!(
        "wavelength" => wavelength,
    ).unwrap()
}

fn join_into_forward_backward(left:LazyFrame, right:LazyFrame, is_forward:bool)->LazyFrame{
    let (suffix,fillstrategy) = if is_forward{
        (format!("_forward"),FillNullStrategy::Forward(None))
    }else{
        (format!("_backward"),FillNullStrategy::Backward(None))
    };

    let mut names:Vec<String> = Vec::new();
    for i in 1..=7 {names.push(format!("mu{}_s",i))}
    for i in 1..=7 {names.push(format!("mu{}_c",i))}
    let new_names:Vec<String> = names.iter().map(|x| format!("{}{}",*x.clone(),suffix)).collect();

    let mut renamed:Vec<Expr> = vec![col("wavelength")];    
    renamed.push(col("wavelength"));
    for i in 1..names.len(){
        renamed.push(col(names[i].clone()).alias(new_names[i]))
    }
    renamed.push(col("wavelength").alias(format!("original_wavelength{}",suffix)));
    let extra_lf = left.clone().select(renamed);


    let include_lf = extra_lf.clone().join(
        right.clone(),
        [col("wavelength")],
        [col("wavelength")],
        JoinArgs::new(JoinType::Full).with_coalesce(JoinCoalesce::CoalesceColumns)
    ).sort(["wavelenght"],Default::default());


    let mut last_exprs:Vec<Expr> = Vec::new();
    last_exprs.push(col("wavelength"));
    
    for name in new_names.into_iter(){
        last_expr.push(col(name).fill_null_with_strategy(strategy))
    }
    last_expr.push(col(format!("original_wavelength{}",suffix)));
    let filled_null_lf = include_lf.clone().select(last_exprs);

    filled_null_lf
}

fn select_only_df_w0(left:LazyFrame,right:LazyFrame)->LazyFrame{
    left.join(
        right.clone(),
        [col("wavelength")],
        [col("wavelength")],
        JoinArgs::new(JoinType::Inner)
    )
}

fn append_fractional_distance(left:LazyFrame,right:LazyFrame)->LazyFrame{

    left.join(
        right.clone(),
        [col("wavelength")],
        [col("wavelength")],
        JoinArgs::new(JoinType::Inner)
    ).with_columns([
        ((col("wavelength") - col("original_wavelength_backward"))
        /(col("original_wavelength_forward")- col("original_wavelength_backward")))
        .alias("fractional_distance")
    ])
}

fn linear_interpolation_full(lf:LazyFrame)->LazyFrame{
    let mut names:Vec<String> = Vec::new();
    for i in 1..=7 {names.push(format!("mu{}_s",i))}
    for i in 1..=7 {names.push(format!("mu{}_c",i))}
    let mut exprs:Vec<Expr> = vec![col("wavelength")];
    for index in 1..=7{
        let expression:Expr = 
        (col(format!("mu{}_s_forward"))*col("fractional_distance")
        + col(format!("mu{}_s_backward")) * ( lit(1.0) - col("fractional_distance")))
        .alias(format!("mu{}_s",index));
        exprs.push(expression);
    }
    for index in 1..=7{
        let expression:Expr = 
        (col(format!("mu{}_c_forward"))*col("fractional_distance")
        + col(format!("mu{}_c_backward")) * ( lit(1.0) - col("fractional_distance")))
        .alias(format!("mu{}_c",index));
        exprs.push(expression);
    }

    lf.select(exprs)
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

    let filter_lower_expr = col("wavelengths").gt(lit(min_wavelength));
    let filter_greater_expr = col("wavelengths").lt(lit(max_wavelength));

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
        let lb_wavelength = col("wavelengths").gt(lit(wavelength*minval_rel_dopplershift - epsilon));
        let ub_wavelength = col("wavelengths").lt(lit(wavelength*maxval_rel_dopplershift + epsilon));
        
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
pub fn filter_wavelength_range(
    grids_lf: LazyFrame,
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64,
)->LazyFrame{
    
    //Nadya's grids are in Angstroms while Joris's are in nm. To check if the requested wavelengths are low resolution i.e. dλ,1e-3nm, I need to specify that or use the same wavelenght units.
    let is_low_resolution=false;
    let combined_expresion=match is_low_resolution{
        true => {filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift).or(
            filter2_sift_wavelengths(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift)
        )}
        false => {filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift)}
    };

    match combined_expresion{
        Some(expresion)=>{grids_lf.filter(expresion)}
        None=>{panic!("unable to produce dataframe using the intensity grid files")}
        
    }
}

pub fn sift_data_frame(
    grids_lf:LazyFrame,
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64,
)->LazyFrame{

    let df_w0 = df!(
        "wavelength"=>wavelengths.clone()
    ).unwrap();

    // a new lazyframe that's a copy of the original one but it has a forward and backward original wavelengths
    let grids_lf_forward = grids_lf.clone().select(
        [
            col("wavelength"),
            col("therest"),
            col("wavelength").alias("original wavelength forward"),
            col("wavelength").alias("original wavelength backward"),
        ]
    );

    let mut join_args=JoinArgs::default();
    join_args.how = JoinType::Full;

    let step5_lf = grids_lf_forward.join(df_w0.lazy(),
    col("wavelength"),
    col(""),
    join_args
    );



    let grids_lf_backward = grids_lf.clone();




}



//add interpolating profile test for each fractional coordinate.
