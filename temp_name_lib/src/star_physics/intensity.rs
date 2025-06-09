use nalgebra as na;

enum GridPointStatus{
    NotAvailable,
    NotLoaded,
    Loaded,
}

enum IntensityType{
    Total,
    ContinuumOnly
}

pub struct Intensity<'a> {
    temp_eff:Vec<f64>,
    log_g:Vec<f64>,
    lambda:Vec<f64>,
    file_name:&'a mut String,//[Ricardo]: I suppose you have different intensity files.
    
    temp_eff_index:usize,
    log_g_index:usize,
    lambda_index:usize,

    lambda_lower:f64,
    lambda_upper:f64,
    limb_darkening_coef:na::Vector4<f64>,
    grid_status:GridPointStatus,
    kind:IntensityType,
}