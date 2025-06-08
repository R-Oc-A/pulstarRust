use crate::type_def::Config;
use crate::joris_math::spherical_harmonics::d_plmcos_dtheta::deriv1_plmcos_dtheta as d_plmcos_dtheta;
use crate::joris_math::spherical_harmonics::d_plmcos_dtheta::deriv2_plmcos_dtheta as d2_plmcos_dtheta;
use crate::joris_math::spherical_harmonics::norm_factor::ylmnorm;

///Computes the derivatives of Δr/r0 with respect to θ in the point with spherical
///coordinates θ,ϕ
pub fn d_dr_rdtheta(parameters: &Config,
                        index: usize,
                        sintheta: f64,
                        costheta: f64,
                        phi: f64) -> f64{
    let r_dr=parameters.rel_deltar[index];
    let l = parameters.l[index];
    let m = parameters.m[index];
    let phase = parameters.phase[index];
    
    r_dr*ylmnorm(l,m)
    * d_plmcos_dtheta(l,m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).cos()
}

///Computes the derivatives of Δθ with respect to θ in the point with spherical
///coordinates θ,ϕ
pub fn d_dtheta_dtheta(parameters: &Config,
                        index: usize,
                        sintheta: f64,
                        costheta: f64,
                        phi: f64) -> f64{

    let r_dr=parameters.rel_deltar[index];
    let k = parameters.k[index];
    let l = parameters.l[index];
    let m = parameters.m[index];
    let phase = parameters.phase[index];
    
    r_dr*ylmnorm(l,m)*k
    * d2_plmcos_dtheta(l,m.abs() as u16,sintheta,costheta)
    * (phase + (m as f64) * phi).cos()
}