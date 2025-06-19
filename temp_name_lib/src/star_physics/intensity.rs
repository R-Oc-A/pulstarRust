use nalgebra as na;

use crate::utils::MathErrors;

enum GridPointStatus{
    NotAvailable,
    NotLoaded,
    Loaded,
}

enum IntensityType{
    Total,
    ContinuumOnly
}

pub struct Intensity {
    temp_eff:Vec<f64>,
    log_g:Vec<f64>,
    lambda:Vec<f64>,
    file_name:String,//[Ricardo]: I suppose you have different intensity files.
    
    temp_eff_index:usize,
    log_g_index:usize,
    lambda_index:usize,

    lambda_lower:f64,
    lambda_upper:f64,
    limb_darkening_coef:na::Vector4<f64>,
    grid_status:GridPointStatus,
    kind:IntensityType,
}


impl Intensity{

    pub fn read_file(){}
    pub fn read_info(){}
    
    pub fn set_wavelength_range(&mut self,lower:f64,upper:f64){
        match lower<upper{
            false => {panic!("lower_bound bigger than upper_bound")}
            true => {self.lambda_lower=lower;
                     self.lambda_upper=upper;}
        }
    }
    pub fn interpolate(){}
}


pub fn hunt(vec:&[f64],x:f64)->Result<usize,MathErrors>{
    //let mut higher:usize = 0;
    //let mut middle:usize = 0;
    //let mut index:usize=0;
    //higher = vec.si
    //is bisection useful here???? I think both are O(n)

    match x<vec[0] || x>*vec.last().expect("empty vector"){
        true => {Err(MathErrors::OutOfBounds)}
        false =>{
            let mut index:usize=0;
            for (n,item) in vec.iter().enumerate(){
                if x>*item {
                    index = n;
                    break;
                }
            }
            Ok(index)
        }
    }
}

