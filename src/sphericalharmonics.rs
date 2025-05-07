use is_odd::IsOdd;




/// Computes the associated Legendre function P_l^m(x) defined by
///     1/2^l/(l!)*(1-x^2)^(m/2) \frac{d^(l+m)}{dx^(l+m)}(x^2-1)^l
/// with x = cos(theta). It is an adapted version of the routine 
/// plgndr() in Numerical Recipes in C, 1992, Press et al., where
/// the factor (-1)^m was removed.
///
/// # Arguments
/// 
/// * `l` - The degree l >= 0
/// * `m` - The azimuthal number m, 0 <= m <= l
/// * `theta`: an angle [rad]
/// 
pub fn plmcos(l: u16, m: u16, theta: f64) -> f64 {

    // The following array [1..13] ([0] is dummy) contains (2 n - 1)!! where j!!
    // denotes the product of all odd integers less than or equal to j.
    // E.g. oddfac[5] = (2 * 5 - 1)!! = 9!! = 945.0

    const MAX_ODDFAC_ARG: usize = 13;
    const MAX_ODDFAC_ARG_U16: u16 = MAX_ODDFAC_ARG as u16;
    const ODDFAC: [f64; MAX_ODDFAC_ARG+1] =  [ 0.0, 1.0, 3.0, 15.0, 105.0, 945.0, 10395.0, 135135.0, 2027025.0, 34459425.0,
                                               654729075.0, 13749.310575e6, 316234.143225e6, 7905853.580625e6 ];

    // Only allow valid values of m

    if m > l {
        panic!("plmcos: m > l");
    }

    // First compute P_m^m(costheta)

    let sintheta = theta.sin();
    let costheta = theta.cos();

    let mut pmmcostheta: f64 = match m {
        0 => 1.0,
        1 => sintheta,
        2 => 3.0 * sintheta * sintheta,
        3 => 15.0 * sintheta * sintheta * sintheta,
        4 => 105.0 * sintheta * sintheta * sintheta * sintheta,
        5..=MAX_ODDFAC_ARG_U16 => ODDFAC[m as usize] * sintheta.powf(m as f64),
        _ => {
                let mut oddfactors = *ODDFAC.last().unwrap();
                for i in (2*(MAX_ODDFAC_ARG+1)-1..=(2 * m as usize - 1)).step_by(2) {
                    oddfactors *= i as f64;
                }
                oddfactors
             },
    };

    // If l == m we're already done

    if l == m {
        return pmmcostheta;
    }

    // If not, use a recurrence relation to compute Plm(costheta)

    let mut pm1mcostheta: f64 = costheta * pmmcostheta * (2 * m + 1) as f64;            // P_{m+1}^m(costheta)

    if l == m+1 {
        return pm1mcostheta;
    }
    else {
        let mut plmcostheta: f64 = 0.0;
        for i in m+2..=l {
            plmcostheta = (costheta * f64::from(2*i-1) * pm1mcostheta - f64::from(i + m - 1) * pmmcostheta) / f64::from(i - m);
            pmmcostheta = pm1mcostheta;
            pm1mcostheta = plmcostheta;
        }

        return plmcostheta;
    }

}









/// Computes the derivative of P_l^m(cos(theta)) with respect to theta, 
/// using a recurrence relation.
///
/// # Arguments
/// 
/// * `l` - The degree l >= 0
/// * `m` - The azimuthal number m, 0 <= m <= l
/// * `theta`: an angle [rad], can not be 0 or a multiple of pi.
///
pub fn deriv1_plmcos(l: u16, m: u16, theta: f64) -> f64 {

    (- f64::from(l+1) * theta.cos() * plmcos(l,m,theta) + f64::from(l-m+1) * plmcos(l+1,m,theta)) / theta.sin()
}









/// Computes the 2nd derivative of P_l^m(cos(theta)) with respect to theta,
/// using a recurrence relation to do so.
///
/// # Arguments
/// 
/// * `l` - The degree l >= 0
/// * `m` - The azimuthal number m, 0 <= m <= l
/// * `theta`: an angle [rad], can not be 0 or a multiple of pi.
///


pub fn deriv2_plmcos(l: u16, m: u16, theta: f64) -> f64 {

    let costheta = theta.cos();
    let sintheta = theta.sin();
    let inv_sqr_sintheta = 1.0 / (sintheta * sintheta);

    f64::from(l+1) * (1.0 + f64::from(l+2) * costheta*costheta*inv_sqr_sintheta) * plmcos(l,m,theta)
    - 2.0 * f64::from(l-m+1) * f64::from(l+2) * costheta * inv_sqr_sintheta * plmcos(l+1,m,theta)
    + f64::from(l-m+1) * f64::from(l-m+2) * inv_sqr_sintheta * plmcos(l+2,m,theta)
}








/// Computes the normalisation factor N_l^m of the spherical harmonic Y_l^m so that:
///          Y_l^m = N_l^m * P_l^{|m|}(cos(theta)) * e^(i m phi)
///
/// # Arguments
///
/// * `l` - The degree l >= 0
/// * `m` - The azimuthal number m, 0 <= m <= l
///
pub fn ylmnorm(l: u16, m: i16) -> f64 {


    // The following array PRECOMPUTED[0..4][0..4] contains: 
    //   if  (0 <= m <= l <= 4): 
    //     PRECOMPUTED[l][m] = sqrt{(2*l + 1)/(4 pi) * (l - m)! / (l + m)!}
    //   else
    //     PRECOMPUTED[l][m] = 0.0

    const PRECOMPUTED: [[f64; 5]; 5] = [[0.282094792, 0.0,         0.0,          0.0,          0.0],
                                        [0.488602512, 0.345494149, 0.0,          0.0,          0.0],
                                        [0.630783131, 0.257516135, 0.128758067,  0.0,          0.0],
                                        [0.746352665, 0.215453456, 0.0681323651, 0.0278149216, 0.0],
                                        [0.846284375, 0.189234939, 0.0446031029, 0.0119206807, 0.00421459707]];

    const INV4PI: f64 = 0.07957747155;      // 1/(4 pi)


    // Panic if m is not between -l and l

    if m.abs() as u16 > l {
        panic!("ylmnorm: |m| > l");
    }

    // For small values of l (<= 4) use the build-in table, but take into account the Condon-Shortley 
    // phase factor, which is 1 for negative m, 1 for positive even m, and -1 for positive odd m.

    if l <= 4 {
        if m > 0 && m.is_odd() { 
            return -1.0 * PRECOMPUTED[l as usize][m.abs() as usize];
        } else {
            return PRECOMPUTED[l as usize][m.abs() as usize];
        }
    } 

    // For other values of l, we first compute the division (l - |m|)!/(l + |m|)!           
    // Note that the number in the faculty in the numerator is always <= than the number in the faculty in the denominator. 
    // So all the factors in the numerator cancel away with some factors in the denominator.                                              

    let mut fac_division: f64 = 1.0;                  // (l - |m|)! / (l + |m|)!
    if m != 0 {
        for i in (l - m.abs() as u16 + 1)..=(l + m.abs() as u16) {
            fac_division *= f64::from(i); 
        }
        fac_division = 1.0 / fac_division; 
    }

    // Then the rest of the Y_l^m norm, including the Condon-Shortley phase factor

    if m > 0 && m.is_odd() {
        return - (INV4PI * f64::from(2*l+1) * fac_division).sqrt();
    } else {
        return   (INV4PI * f64::from(2*l+1) * fac_division).sqrt();

    }
}













// pub fn lnfac(n: u32) -> f64 {
//
// }










#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_plmcos() {
        let theta = 2.6;                                             // [rad]
        let computed = plmcos(4,4,theta);
        let expected = 105.0 * theta.sin().powf(4.0);
        assert_approx_eq!(computed, expected, 1.0e-10);

        let theta = -0.34;
        let computed = plmcos(3, 2, theta);
        let expected = 15.0 * theta.cos() * theta.sin() * theta.sin(); 
        assert_approx_eq!(computed, expected, 1.0e-10);
    }

    #[test]
    #[should_panic]
    fn test_plmcos_panic() {
        let theta = 2.6;
        plmcos(4,5,theta);
    }
}

