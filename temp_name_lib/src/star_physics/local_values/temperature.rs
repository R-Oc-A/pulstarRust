use crate::{star_physics::local_values::local_value, type_def::{Config,Eigenfunctions}, utils::MathErrors};

///This function returns the local time dependent temperature (in units of temp_eff0).
///where temp_eff0 is the equilibrium effective temperature.
pub fn local_temp(temperature: &[Eigenfunctions],
                parameters:&Config,
                theta:f64,
                phi:f64,
                temp_eff0:f64)->Result<f64,MathErrors>{
    local_value(temperature, parameters, theta, phi, temp_eff0)
}