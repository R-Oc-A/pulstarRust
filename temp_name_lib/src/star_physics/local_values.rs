use crate::type_def::{Eigenfunctions,Config};
use crate::joris_math::geometry::displacement::displacement;
use crate::joris_math::spherical_harmonics::norm_factor::ylmnorm;
use crate::utils::MathErrors;

pub mod gravity;
pub mod temperature;

//Tells if the measurement has non trivial values. Internal use only.
fn non_zero_ampl(g: &[Eigenfunctions])-> bool{
    let mut res=false;
    for item in g{
        match item.ampl != -1.0 {
            true =>{
                res = true;
                break;}
            false =>{res = false;}
        }
    }
    return res;
}


//This function calculates the local time dependent values of the
//Eigenfunctions (Local gravity or Temperature) in units of some
//base value. Internal use only.
pub fn local_value(value: &[Eigenfunctions],
                   parameters: &Config,
                   theta:f64,
                   phi:f64,
                   base_unit:f64)-> Result<f64,MathErrors> {

    match non_zero_ampl(value){
        false => {Ok(base_unit)}

        true =>{
            let mut res = 0.0;
            let sintheta=theta.sin();
            let costheta=theta.cos();
            for (n,item) in value.iter().enumerate(){
                if item.ampl != 0.0{
                    let ampl_radial= parameters.rel_deltar[n] 
                                        * ylmnorm(parameters.l[n], parameters.m[n]);
                    let ampl_tangential= parameters.rel_deltar[n] 
                                            * parameters.k[n] * 
                                            ylmnorm(parameters.l[n], parameters.m[n]);
                    
                    let dif_phase=Config{
                        n_modes:0,
                        l:vec![parameters.l[n]],
                        m:vec![parameters.m[n]],
                        rel_deltar:vec![parameters.rel_deltar[n]],
                        k:vec![parameters.k[n]],
                        phase:vec![parameters.phase[n]+value[n].phasedif],
                    };

                    let ds = displacement(&dif_phase, 
                        0,
                        sintheta,
                        costheta,
                        phi,
                        ampl_radial,
                        ampl_tangential) ? ;

                    res += item.ampl * ds.coords[0];
                }
            }

            Ok(res)
        }

    }
}