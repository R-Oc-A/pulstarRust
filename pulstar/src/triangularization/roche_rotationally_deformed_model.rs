use crate::triangularization::roche_rotationally_deformed_model;
use crate::{MeshConfig::DSphere, PulstarConfig, na::Vector3};
use temp_name_lib::{type_def::{GRAVCONSTANT, MASSSUN, RADIUSSUN,PI}, utils::MACHINE_PRECISION};


//TODO: Check units
pub fn unperturbed_local_g_en_teff(parameters:&PulstarConfig,point_coords:&Vector3<f64>,g0:f64,t_eff:f64)->(f64,f64){
    let b = compute_b_dimensionless_quantity(parameters);
    let stheta= point_coords.z/point_coords.norm();//assuming cartesian coordinates
    let mut theta = stheta.asin();
    if theta >=0.0{
                        theta = (theta-0.5*PI).abs();
                    }else{
                        theta = theta.abs() + 0.5*PI;
                    };
    let s2theta  = (2.0*theta).sin();
    //let stheta = point_coords.y.sin(); //assuming spherical coordinates.
    //let s2theta = (2.0*point_coords.y).sin();//assuming spherical coordinates
    let accuracy = MACHINE_PRECISION;
    let x = compute_r_rp_ratio(b, stheta, accuracy);

    let beta = 0.25;
    
    let g = g0 * (
        (1.0/x.powi(2) - 2.0*b*x*stheta.powi(2)).powi(2)
        + (b * x * s2theta).powi(2)
    ).sqrt();

    let temperature = t_eff*(g/g0).powf(beta);

    (g.log10(), temperature)

}

pub fn compute_b_dimensionless_quantity(parameters:&PulstarConfig)->f64{
    //let rotation_frequency=//parameters.get_rotation_frequency()*CYCLI2RAD;//in Hz
    let mut rotation_frequency = roche_rotationally_deformed_model::rotation_frequency(parameters);
    let G= GRAVCONSTANT;//in something
    let R = parameters.star_data.radius * RADIUSSUN; //in m I guess, and this is the polar radius
    let M = parameters.star_data.mass*MASSSUN;
    rotation_frequency *= (G*M/R.powi(3)).sqrt();
    let u_pow2 = 27.0/8.0 * R.powi(3)*rotation_frequency.powi(2)/(G*M);
    4.0*u_pow2/27.0
}

//Adapted from Wenjin Huang code
// ratio radius/polar_radius
fn compute_r_rp_ratio(beta:f64,stheta:f64,accuracy:f64)->f64{
    if stheta.abs()<=MACHINE_PRECISION {return 1.0*stheta.signum()}
    let a = beta * stheta.powi(2);
    let mut x = (1.0-2.0*a)/(1.0-3.0*a);
    let mut eps = 1.0;
    
    while (eps>accuracy){
        let x_new = (x + (x-1.5)/(3.0*a*x.powi(2) - 1.0))/1.5;
        eps = (x_new -x).abs();
        x = x_new;
    }
    x
}

/// Returns the rotation frequency in Hz (s^{-1})
fn rotation_frequency(parameters:&PulstarConfig)->f64{
    let mut rotation_frequency = match parameters.mesh{
        DSphere{ triangle_length:_,rotation_frequency:omega}=>{omega},
        _=>{panic!("this function should only be called for the Deformed star case")}
    };
    let G= GRAVCONSTANT;//in
    let R = parameters.star_data.radius * RADIUSSUN; //in m I guess, and this is the polar radius
    let M = parameters.star_data.mass*MASSSUN;
    rotation_frequency *= (G*M/ R.powi(3)).sqrt();
    rotation_frequency
}
///Returns the rotation velocity in km/s. This value is dependant on the colatitudinal angle θ
///### Arguments:
/// * `parameters` - An instance of [PulstarConfig] that contains the global information of the star. 
/// * `θ` - Colatitudinal angle
/// ### Returns:
/// * `rotational_velocity` - A [f64] value that contains the rotational velocity of a given point on a surface of the star assuming solid like rotation (equal frequency of rotation) depending on the radius of the star. 
pub fn get_rotation_velocity(parameters:&PulstarConfig,theta:f64)->f64{
    let rotation_frequency = roche_rotationally_deformed_model::rotation_frequency(
        parameters
    );
    let b = compute_b_dimensionless_quantity(parameters);
    let stheta = theta.sin();
    let accuracy = 1.0e-8;
    let x = compute_r_rp_ratio(b, stheta, accuracy);
    let mut r_sintheta = x*parameters.star_data.radius * RADIUSSUN * stheta;//radius in m
    r_sintheta *= 1.0e-3;//radius in km
    r_sintheta * rotation_frequency
}