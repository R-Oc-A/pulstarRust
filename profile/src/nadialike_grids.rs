//Things to do 
// I need to change the struct GridsData to have the same structure as the one here below
// This list approach I have for the intensity grids might not be the best as 
// I think it's best to stick to 4 intensity grids that cover the whole parameter space
// I also have to enumerate the intensity grids. 
// Later it might be necessary to 

use ndarray;
use polars::prelude::*;
/// This structure holds the information contained on Nadya's specific intensity grids. 
#[derive(Clone)]
pub struct SpectralGrid {
    mu_values:[f64;7],
    t_eff:[f64;2],
    log_g:[f64;2],
    grid_values:ndarray::Array3<f64>,
    wavelengths: Vec<f64>
}


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


//

