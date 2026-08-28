use polars::frame::DataFrame;
use serde::{Deserialize,Serialize};
use super::GaussianProfile;
use std::f64::consts::PI;
use temp_name_lib::type_def::{CLIGHT, GRAVCONSTANT, MASSSUN, RADIUSSUN};

#[derive(Deserialize,Serialize)]
pub struct GaussianProfileInit{
    ///intrinsic width
    sigma:f64,
    ///Equivalent width
    eq_w:f64,
    ///Alpha W
    alpha_w:f64,
    ///Zero point shift of velocity
    zero_point_shift:f64,
    ///Central Wavelength
    central_wavelength:f64,
    //leftmost wavelength
    left_wavelength:f64,
    //Rightmost wavelength
    right_wavelength:f64,
    //Lambda_resolution
    step:f64,
    // Star temperature
    t_eff:f64,
    //Star mass
    mass:f64,
    //Star Radius
    radius:f64,
}
//Things I need to initialize:
// * fl_in_ul
// * y_gauss
// * wavelength
// * sigma_sqrtpi_sqrt2
// * sigma_sqrt2_pow2
// * Star temperature
// * Star logg
// * Limb darkening coefficients
// * Output df

pub fn init_profile(toml_string:&str)->GaussianProfile{

    let config = GaussianProfileInit::read_from_toml(toml_string);
    let wavelength = config.init_wavelength_arr();
    let y_gauss = vec![0.0;wavelength.len()];
    let fl_in_ul = vec![0.0;wavelength.len()];

    let eq_w_lmbd = config.central_wavelength - from_kms_to_lambda(config.eq_w, config.central_wavelength);
    let sigma_lmbd = config.central_wavelength - from_kms_to_lambda(config.sigma, config.central_wavelength);


    let sigma_sqrtpi_sqrt2 = 1.0/(sigma_lmbd * PI.sqrt() * 2.0f64.sqrt());
    let sigma_sqrt2_pow2 = 0.5/sigma_lmbd.powi(2);
    //surface gravity
    let logg = 4.438 + config.mass.log10() - 2.0*config.radius.log10();
    // Init limb_coeffs


    GaussianProfile{
        fl_in_ul:fl_in_ul,
        continuum:0.0,
        y_gauss:y_gauss,
        wavelength:wavelength,
        sigmag_sqrt2_pow2:sigma_sqrt2_pow2,
        sigmag_sqrtpi_sqrt2:sigma_sqrtpi_sqrt2,
        eq_w:eq_w_lmbd,
        alpha_w:config.alpha_w,
        zero_point_shift:config.zero_point_shift,
        central_wavelength:config.central_wavelength,
        t_eff:config.t_eff,
        log_g:logg,
        output:DataFrame::empty(),
        time_point:0.0,
        limb: super::LimbDarkeningCoefficients([0.0;4])
    }

}
fn from_kms_to_lambda(vel:f64,lambda_0:f64)->f64{
    lambda_0 * (1.0/(1.0-vel/CLIGHT *1.0e3))
}
impl GaussianProfileInit{
    
    fn read_from_toml(toml_string:&str)->Self{
        toml::from_str(toml_string).expect("error parsing toml for profile config") /*{
        
            Ok(profile)=>{profile}
            Err(error)=>{panic!("error parsing toml for profile config{}",error)}
        }*/
    }

    fn init_wavelength_arr(&self)->Vec<f64>{
        let npts = ((self.right_wavelength-self.left_wavelength)/self.step).floor() as usize;
        let mut wavelength_arr:Vec<f64> = Vec::with_capacity(npts+1);
        for i in 0..=npts{
            wavelength_arr.push(self.left_wavelength + self.step * (i as f64));
        }
        wavelength_arr
    }


    
}