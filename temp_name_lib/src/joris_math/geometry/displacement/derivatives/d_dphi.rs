use crate::type_def::Config;
use crate::joris_math::spherical_harmonics::{plmcos::plmcos,
                                        norm_factor::ylmnorm};


///Computes the derivatives of Δr/r0 with respect to ϕ in the point with spherical
///coordinates θ,ϕ
pub fn d_dr_rdphi(parameters: &Config,
                        index: usize,
                        theta: f64,
                        phi: f64) -> f64{
    let r_dr=parameters.rel_deltar[index];
    let l = parameters.l[index];
    let m = parameters.m[index];
    let phase = parameters.phase[index];
    
    r_dr*ylmnorm(l,m) * (-m as f64)
    * plmcos(l, m.abs() as u16, theta.sin(), theta.cos())
    * (phase + (m as f64) * phi).sin()
}
voi aqyu
///Computes the derivatives of Δϕ with respect to ϕ in the point with spherical
///coordinates θ,ϕ
pub fn d_dphi_dphi(parameters: &Config,
                        index: usize,
                        theta: f64,
                        phi: f64) -> f64{
    let r_dr=parameters.rel_deltar[index];
    let l = parameters.l[index];
    let m = parameters.m[index];
    let phase = parameters.phase[index];
    
    r_dr*ylmnorm(l,m) * (-m as f64)
    * plmcos(l, m.abs() as u16, theta.sin(), theta.cos())
    * (phase + (m as f64) * phi).sin()
}