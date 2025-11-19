use polars::prelude::*;

/// This module is used to transform a Grid from joris into a regular grid. 
/// 
/// 

/// This structure holds the information contained on Nadya's specific intensity grids. 


//Nadya's mu values 0.9636 0.8864 0.7071 0.5976 0.4629 0.2673

//As with Joris'like grids, I'm going to divide this module into three parts

// The first part will consist of the methods to produce a [DataFrame]
// that stores the contents of nadyalike intensity grids. 
// 
// The second part will consist of the methods to fill in the spectral grid data structure
//
// The third part will consist on the methods to produce the linear interpolation


//--------------------------------------------------
//---------FIRST PART-------------------------------
//--------------------------------------------------
//I have a lot of functions to open the intensity grids. So let's leave it like that. 
// I'll focus on parsing them. 


        /*match self{
            Self::Joris { temperature:_, log_gravity:_, filename:_ }=>{
                collection.append(& mut extract_column_as_vectorf64("a", &grid_df));
                collection.append(& mut extract_column_as_vectorf64("b", &grid_df));
                collection.append(& mut extract_column_as_vectorf64("c", &grid_df));
                collection.append(& mut extract_column_as_vectorf64("d", &grid_df));
                collection.append(& mut extract_column_as_vectorf64("ac", &grid_df));
                collection.append(& mut extract_column_as_vectorf64("bc", &grid_df));
                collection.append(& mut extract_column_as_vectorf64("cc", &grid_df));
                collection.append(& mut extract_column_as_vectorf64("dc", &grid_df));
            }*/
//

pub fn convert_joris_grid_to_regular_grid(joris_lf:LazyFrame)->LazyFrame{
    let mu_values=[0.9636,0.8864,0.8018,0.7071,0.5976,0.4629,0.2623];

    joris_lf.clone().select([
        col("wavelengths"),
        append_specific_intensity(mu_values[0]).alias("mu1_s"),
        append_specific_intensity(mu_values[1]).alias("mu2_s"),
        append_specific_intensity(mu_values[2]).alias("mu3_s"),
        append_specific_intensity(mu_values[3]).alias("mu4_s"),
        append_specific_intensity(mu_values[4]).alias("mu5_s"),
        append_specific_intensity(mu_values[5]).alias("mu6_s"),
        append_specific_intensity(mu_values[6]).alias("mu7_s"),
        append_continuum_intensity(mu_values[0]).alias("mu1_c"),
        append_continuum_intensity(mu_values[1]).alias("mu2_c"),
        append_continuum_intensity(mu_values[2]).alias("mu3_c"),
        append_continuum_intensity(mu_values[3]).alias("mu4_c"),
        append_continuum_intensity(mu_values[4]).alias("mu5_c"),
        append_continuum_intensity(mu_values[5]).alias("mu6_c"),
        append_continuum_intensity(mu_values[6]).alias("mu7_c"),
    ])
    
}

fn append_specific_intensity(mu:f64)->Expr{
    col("a") + col("b") * lit(1.0-mu) + col("c") * lit(1.0-mu.powi(2)) + col("d") * lit(1.0-mu.powi(3))
}
fn append_continuum_intensity(mu:f64)->Expr{
    col("ac") + col("bc") * lit(1.0-mu) + col("cc") * lit(1.0-mu.powi(2)) + col("dc") * lit(1.0-mu.powi(3))
}