use super::plmcos::plmcos;
/// Computes the derivative of P_l^m(cos(theta)) with respect to theta, 
/// using a recurrence relation.
///
/// # Arguments
/// 
/// * `l` - The degree l >= 0
/// * `m` - The azimuthal number m, 0 <= m <= l
/// * `sintheta`: sin(theta), can not be 0
/// * `costheta`: cos(theta)
///
pub fn deriv1_plmcos_dtheta(l: u16, m: u16, sintheta: f64, costheta: f64) -> f64 {

    (- f64::from(l+1) * costheta * plmcos(l,m,sintheta,costheta) + f64::from(l-m+1) * plmcos(l+1,m,sintheta,costheta)) / sintheta
}



/// Computes the 2nd derivative of P_l^m(cos(theta)) with respect to theta,
/// using a recurrence relation to do so.
///
/// # Arguments
/// 
/// * `l` - The degree l >= 0
/// * `m` - The azimuthal number m, 0 <= m <= l
/// * `sintheta`: sin(theta)
/// * `costheta`: cos(theta)
///

pub fn deriv2_plmcos_dtheta(l: u16, m: u16, sintheta: f64, costheta: f64) -> f64 {

    let inv_sqr_sintheta = 1.0 / (sintheta * sintheta);

    f64::from(l+1) * (1.0 + f64::from(l+2) * costheta*costheta*inv_sqr_sintheta) * plmcos(l,m,sintheta,costheta)
    - 2.0 * f64::from(l-m+1) * f64::from(l+2) * costheta * inv_sqr_sintheta * plmcos(l+1,m,sintheta,costheta)
    + f64::from(l-m+1) * f64::from(l-m+2) * inv_sqr_sintheta * plmcos(l+2,m,sintheta,costheta)
}



