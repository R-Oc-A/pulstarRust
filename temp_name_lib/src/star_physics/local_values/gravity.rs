use super::local_value;
use crate::{type_def::{Config,Eigenfunctions}};
use crate::utils::MathErrors;
///This function returns the lograithm (base10) of the local time dependent gravity
///(in units of g0). log_g0 is the logarithm(base 10) of the equilibrium gravity.
///  
/// INPUT: . pulspar: pulsation parameters of the (different) mode(s).
///        . g[].ampl: contains in fact not the amplitude for delta g/g but 
///                    contains a factor which the amplitude of xi_r/r should be
///                    multiplied with to obtain the true amplitude of delta g/g.
///          g[].phasedif: contains the phaselag (i.e. the phasedifference with 
///                        xi_r/r) expressed in radians, i.e. 
///                        phase(dg/g) - phase(dr/r) )
///        
///        . thetarad : theta coordinate (radians)
///        . phirad   : phi coordinate (radians)
///
///
pub fn log_gravity(g: &[Eigenfunctions],
                    parameters: &Config,
                    theta:f64,
                    phi:f64,
                    log_g0:f64) -> Result<f64,MathErrors>{

    local_value(g, parameters, theta, phi, log_g0)
}














