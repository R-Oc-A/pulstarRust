use crate::joris_math::spherical_harmonics::d_plmcos_dtheta::deriv1_plmcos_dtheta;
use crate::joris_math::spherical_harmonics::norm_factor::ylmnorm;
use crate::type_def::{Config, VectorComponents,VectorBase};
use crate::joris_math::spherical_harmonics::plmcos::plmcos;
use crate::utils::{self, MathErrors};
use nalgebra as na;


///Computes v_r, v_θ, v_φ of the pulsational velocity.
pub fn v_puls(parameters:&Config, //pulsation parameters
                index:usize,//radial order
                sintheta:f64,
                costheta:f64,
                phi:f64,
                velocity_amplitude:&[f64],
                ) -> Result<VectorComponents,MathErrors>{

    match sintheta.abs() <= utils::MACHINE_PRECISION{
        true =>{Err(MathErrors::DivisionByZero)}

        false => {
            let l = parameters.l[index];
            let m = parameters.m[index];
            let phase = parameters.phase[index];
            let k = parameters.k[index];

            let v_r = velocity_amplitude[index] * ylmnorm(l, m)
                  * plmcos(l, m.abs() as u16, sintheta, costheta)
                  * (phase + (m as f64) * phi).sin();
            let v_theta = velocity_amplitude[index] * k
                   * ylmnorm(l, m)
                   * deriv1_plmcos_dtheta(l, m.abs() as u16, sintheta, costheta)
                   * (phase + (m as f64) * phi).sin();
            let v_phi = velocity_amplitude[index] * k
                   * ylmnorm(l, m)
                   * (-(m as f64))
                   * plmcos(l, m.abs() as u16, sintheta, costheta)
                   * (phase * (m as f64) * phi).cos()
                   / sintheta;

            Ok(VectorComponents{
                base: VectorBase::Spherical,
                coords: na::Vector3::new(v_r,v_theta,v_phi)
            })
        }
    }

}

