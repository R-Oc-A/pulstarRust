use crate::type_def::Config;
use crate::joris_math::spherical_harmonics::{plmcos::plmcos,
                                        norm_factor::ylmnorm};
use crate::utils::MathErrors;

///Computes the derivatives of Δr/r0 with respect to ϕ in the point with spherical
///coordinates θ,ϕ
pub fn d_dr_rdphi(parameters: &Config,
                        index: usize,
                        sintheta: f64,
                        costheta:f64,
                        phi: f64) -> f64{

    let r_dr=parameters.rel_deltar[index];
    let l = parameters.l[index];
    let m = parameters.m[index];
    let phase = parameters.phase[index];
    
    r_dr * ylmnorm(l,m) * (-m as f64)
    * plmcos(l, m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).sin()
}

///Computes the derivatives of Δϕ with respect to ϕ in the point with spherical
///coordinates θ,ϕ
pub fn d_dphi_dphi(parameters: &Config,
                        index: usize,
                        sintheta: f64,
                        costheta: f64,
                        phi: f64) -> Result<f64,MathErrors>{
    
    let machine_precision_value= sintheta.abs()<1.0e-8;

    match machine_precision_value{  
        false => {
        let r_dr=parameters.rel_deltar[index];
        let l = parameters.l[index];
        let m = parameters.m[index];
        let phase = parameters.phase[index];
        let k = parameters.k[index];

        Ok(r_dr * k * ylmnorm(l, m) * (-(m as f64).powi(2))
        * plmcos(l, m.abs() as u16, sintheta, costheta)
        * (phase + (m as f64) * phi).cos()
        /(sintheta.abs().powi(2)) )
        }

        true =>{
        Err(MathErrors::DivisionByZero) //will pass the error in order for the calling function to do something
        }
    }
}