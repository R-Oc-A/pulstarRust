use is_odd::IsOdd;
/// Computes the normalisation factor N_l^m of the spherical harmonic Y_l^m so that:
///          Y_l^m = N_l^m * P_l^{|m|}(cos(theta)) * e^(i m phi)
///
/// # Arguments
///
/// * `l` - The degree l >= 0
/// * `m` - The azimuthal number m, 0 <= m <= l
///
pub fn ylmnorm(l: u16, m: i16) -> f64 {

    // Panic if m is not between -l and l

    assert!(m.abs() as u16 <= l, "ylmnorm: |m| > l");

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