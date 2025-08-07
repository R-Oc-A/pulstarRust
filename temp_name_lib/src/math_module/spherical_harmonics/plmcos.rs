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
/// * `sintheta`: sin(theta)
/// * `costheta`: cos(theta)
/// 
pub fn plmcos(l: u16, m: u16, sintheta: f64, costheta: f64) -> f64 {

    // Only allow valid values of m

    assert!(m <= l, "plmcos: m > l");

    // The following array [1..13] ([0] is dummy) contains (2 n - 1)!! where j!!
    // denotes the product of all odd integers less than or equal to j.
    // E.g. oddfac[5] = (2 * 5 - 1)!! = 9!! = 945.0

    const MAX_ODDFAC_ARG: usize = 13;
    const MAX_ODDFAC_ARG_U16: u16 = MAX_ODDFAC_ARG as u16;
    const ODDFAC: [f64; MAX_ODDFAC_ARG+1] =  [ 0.0, 1.0, 3.0, 15.0, 105.0, 945.0, 10395.0, 135135.0, 2027025.0, 34459425.0,
                                               654729075.0, 13749.310575e6, 316234.143225e6, 7905853.580625e6 ];

    // First compute P_m^m(costheta)

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







#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    #[ignore]
    fn test_plmcos() {
        let theta: f64 = 2.6;                                             // [rad]
        let (sintheta, costheta) = (theta.sin(), theta.cos());
        let computed = plmcos(4, 4, sintheta, costheta);
        let expected = 105.0 * theta.sin().powf(4.0);
        assert_approx_eq!(computed, expected, 1.0e-10);

        let theta: f64 = -0.34;
        let (sintheta, costheta) = (theta.sin(), theta.cos());
        let computed = plmcos(3, 2, sintheta, costheta);
        let expected = 15.0 * theta.cos() * theta.sin() * theta.sin(); 
        assert_approx_eq!(computed, expected, 1.0e-10);
    }

    #[test]
    #[ignore]
    #[should_panic]
    fn test_plmcos_panic() {
        let theta: f64 = 2.6;
        let (sintheta, costheta) = (theta.sin(), theta.cos());
        plmcos(4, 5, sintheta, costheta);
    }
}