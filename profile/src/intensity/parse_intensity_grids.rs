use super::*;

pub mod joris_grids;

fn add_padding_for_wavelength(wavelength:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64,
)->Vec<f64>{

    let mut padded_wavelength_vec:Vec<f64> = Vec::new();

    let min_wavelength = if let Some( start) = wavelength.get(0){
        start * (minval_rel_dopplershift)
    }else{
        panic!("What are you trying to do? Your wavelength array has length 0")
    };
    let max_wavelength = if let Some( end) = wavelength.last(){
        end * (maxval_rel_dopplershift)
    }else{
        panic!("What are you trying to do? Your wavelength array has length 0")
    };
    let last_wavelength = if let Some(value) = wavelength.last(){value}else{panic!()};
    
    let d_lambda = wavelength[1]-wavelength[0];

    let left_padding = {
        let mut index = 0usize;
        while min_wavelength < (wavelength[0]-(index as f64)*d_lambda){
            //println!("min_wavelength is {}",min_wavelength);
            //println!("current padding index {}",index);
            index += 1;
        }
        index
    };

    
    let right_padding = {
        let mut index = 0usize;
        while max_wavelength > (last_wavelength + (index as f64)*d_lambda){
        index += 1;
        }
        index
    };

    for index in (1..=left_padding).rev(){
        padded_wavelength_vec.push(wavelength[0] - (index as f64)*d_lambda);
    }

    for lambda in wavelength.iter(){
        padded_wavelength_vec.push(*lambda);
    }
    
    for index in 1..=right_padding{
        padded_wavelength_vec.push(last_wavelength + (index as f64)*d_lambda);
    }
    padded_wavelength_vec
}


fn make_df_w0(wavelength:&[f64])->DataFrame{
    df!(
        "wavelength" => wavelength,
    ).unwrap()
}

fn join_into_forward_backward(left:LazyFrame,
    right:LazyFrame,
    is_forward:bool,
    is_global:bool)->LazyFrame{
    let (suffix,fillstrategy) = if is_forward{
        (format!("_forward"),FillNullStrategy::Forward(None))
    }else{
        (format!("_backward"),FillNullStrategy::Backward(None))
    };

    
    let mut names:Vec<String> = Vec::new();
    if is_global{
        for i in 1..=7 {names.push(format!("mu{}_s",i))}
        for i in 1..=7 {names.push(format!("mu{}_c",i))}
    }else{
        names.push(format!("mu_avg_s"));
        names.push(format!("mu_avg_c"));
    }

    let new_names:Vec<String> = names.iter().map(|x| format!("{}{}",x.clone(),suffix)).collect();

    let mut renamed:Vec<Expr> = vec![col("wavelength")];
    //if !is_global{
    //    renamed.push(col("pixel_id"))
    //};
    for i in 0..names.len(){
        renamed.push(col(names[i].clone()).alias(new_names[i].clone()))
    }
    renamed.push(col("wavelength").alias(format!("original_wavelength{}",suffix)));
    let extra_lf = left.clone().select(renamed);
    //Here it is 
    let include_lf = extra_lf.clone().join(
        right.clone(),
        [col("wavelength")],
        [col("wavelength")],
        JoinArgs::new(JoinType::Full).with_coalesce(JoinCoalesce::CoalesceColumns)
    );
    let include_lf = include_lf.sort(["wavelength"], 
Default::default());

    
    let mut last_exprs:Vec<Expr> = Vec::new();
    last_exprs.push(col("wavelength"));

    if !is_global{last_exprs.push(col("pixel_id"))};//fill_null_with_strategy(fillstrategy))};

    for name in new_names.into_iter(){
        last_exprs.push(col(name).fill_null_with_strategy(fillstrategy))
    }
    last_exprs.push(col(format!("original_wavelength{}",suffix))
    .fill_null_with_strategy(fillstrategy));

    let filled_null_lf = include_lf.clone().select(last_exprs);
    filled_null_lf

}

fn select_only_df_w0(left:LazyFrame,right:LazyFrame)->LazyFrame{
    let selected_lf = left.clone().join(
        right.clone(),
        [col("wavelength")],
        [col("wavelength")],
        JoinArgs::new(JoinType::Inner)
    );
    selected_lf
}

fn append_fractional_distance(left:LazyFrame,right:LazyFrame,is_global:bool)->LazyFrame{

    let columns_to_be_joined:Vec<Expr> = if is_global{vec![col("wavelength")]}else{vec![col("wavelength"),col("pixel_id")]};

    let joined_lf = left.clone().join(
        right.clone(),
        &columns_to_be_joined,
        &columns_to_be_joined,
        JoinArgs::new(JoinType::Inner)
    );

    
    let interpolation_mask_lf = joined_lf.clone().with_columns([
        col("wavelength").neq( col(format!("original_wavelength_forward"))).alias("requires_interpolation")]);

    let fractional_distance_lf = interpolation_mask_lf.clone().with_columns([
        when( col("requires_interpolation"))
        .then(
        (col("wavelength") - col("original_wavelength_forward"))
        /(col("original_wavelength_backward")- col("original_wavelength_forward"))
        )
        .otherwise(lit(0.0))
        .alias("fractional_distance")
    ]);

    fractional_distance_lf
}

fn linear_interpolation_full(lf:LazyFrame,is_global:bool)->LazyFrame{
    let mut exprs:Vec<Expr> = Vec::new();

    let expr = |forward:&str,backward:&str,name:&str|->Expr {
        (col(forward) * col("fractional_distance")
        + col(backward) * (lit(1.0) - col("fractional_distance") ) )
        .alias(name)
    };
    
    exprs.push(col("wavelength"));
    if !is_global{exprs.push(col("pixel_id"))};
    if is_global{
        for char in ["s","c"]{
            for index in 1..=7{
                let name = format!("mu{}_{}",index,char);
                let forward = format!("mu{}_{}_forward",index,char);
                let backward = format!("mu{}_{}_backward",index,char);
                exprs.push(
                    expr(&forward,&backward,&name)
                )
            }
        }
    }else{
        for char in ["s","c"]{
            let name  = format!("mu_avg_{}",char);
            let forward = format!("mu_avg_{}_forward",char);
            let backward = format!("mu_avg_{}_backward",char);
            exprs.push(
                expr(&forward,&backward,&name)
            )
        }
    }

    let result_lf= lf.clone().select( exprs);

    result_lf
}


/// This is the second filter, it's used so that the loaded grids can be linearly interpolated by bulk, 
/// in the sense that no extra queries should be implemented to look for the appropriate grid values that encompas an observed wavelength.
/// 
/// ### Arguments: 
/// * `wavelengths` - a [&[f64]] collection that holds the observed wavelengths.
/// * `maxval_rel_dopplershift` - a [f64] value that stores the maximum value of the Doppler shift.
/// * `minval_rel_dopplershift` - a [f64] value that stores the minimum value of the Doppler shift.
/// ### Returns:
/// * A [DataFrame] that has the structure
/// 
/// `|wavelength|mu1_s|..|mu7_s|mu1_c|..|mu7_c|`
/// and that has the same resolving power of the observable wavelength. 
pub fn sift_dataframe(
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64,
    grids_lf:LazyFrame,
)->DataFrame{

    let padded_wavelength = add_padding_for_wavelength(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift);

    let df_w0 = make_df_w0(&padded_wavelength);
    let lf_w0 = df_w0.lazy();
  
    let forward = join_into_forward_backward(grids_lf.clone(), lf_w0.clone(), true,true);
    let backward = join_into_forward_backward(grids_lf.clone(), lf_w0.clone(), false,true);
    let forward_df = forward.collect().unwrap();
    let backward_df = backward.collect().unwrap();

    let forward = select_only_df_w0(forward_df.lazy().clone(), lf_w0.clone());
    let backward = select_only_df_w0(backward_df.lazy().clone(), lf_w0.clone());
    let lf_fractional_distance = append_fractional_distance(forward.clone(), backward.clone(),true);

    let linear_lf = linear_interpolation_full(lf_fractional_distance.clone(),true);

    linear_lf.collect().unwrap()

}

// I need here to fixe column names
// But I should do it carefully

pub fn wavelength_interpolation(
    shifted_wavelengths:LazyFrame,
    grids_lf: LazyFrame,
    )->LazyFrame{

    //let df_w0 = make_df_w0(shifted_wavelengths);
    let lf_w0 = shifted_wavelengths.clone().select([col("pixel_id"),col("wavelength")]);
    let forward = join_into_forward_backward(grids_lf.clone(), lf_w0.clone(), true,false);
    let backward = join_into_forward_backward(grids_lf.clone(), lf_w0.clone(), false,false);
    let forward = select_only_df_w0(forward.clone(), lf_w0.clone());
    let backward = select_only_df_w0(backward.clone(), lf_w0.clone());
    let lf_fractional_distance = append_fractional_distance(forward.clone(), backward.clone(),false);
    let linear_lf = linear_interpolation_full(lf_fractional_distance.clone(),false); 
    
    //println!("linear_lf {:?}",linear_lf.clone().collect().unwrap());

    linear_lf.select([col("pixel_id"),
        col("wavelength"),
        col("mu_avg_s"),
        col("mu_avg_c")])

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
fn filter1_if_contains_wavelenghts(
    wavelengths:&[f64],
    maxval_rel_dopplershift:f64,
    minval_rel_dopplershift:f64)->Option<Expr>{

    
    let epsilon = 10.0 * (wavelengths.get(0).expect("empty wavelength array")
    -wavelengths.get(1).expect("wavelength array of only one element not supported"));

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
    
    let filter1 = filter1_if_contains_wavelenghts(wavelengths, maxval_rel_dopplershift, minval_rel_dopplershift);

    match filter1{
        Some(expresion)=>{grids_lf.filter(expresion)}
        None=>{panic!("unable to produce dataframe using the intensity grid files")}
        
    }
}

