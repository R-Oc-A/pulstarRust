use crate::math_module::lnfac::lnfac;
use crate::math_module::binomial::binomial;
use is_odd::IsOdd;
/// Computes the quantity d^{(l)}_{km}(angle) as defined by Condon & Odabasi 
/// (1980, "Atomic Structure" ISBN 0 521 21859 4). 
///
/// These function can be used to transfrom a spherical harmonic Y_l^m between a reference frame
/// where the z'-axis points towards the observer and a reference frame were the z-axis is the rotation axis:
/// Y_l^m(\theta,\phi) = \sum_{k=-l}{+l} d^{(l)}_{km}(i) Y_l^k(\theta', \phi')
///
/// # Arguments:
/// * l - degree of pulsation: l >= 0
/// * k - -l <= k <= l
/// * m - azimuthal number: -l <= m <= l
/// * angle - angle in radians
///
pub fn dlkm(l: u32, k: i32, m: i32, angle: f64) -> f64 {

    // Although l should be a positive integer of type u32, we will often need to subtract or
    // compare it with i32 values, for which Rust complains. So cast it once to an i32 value.  

    let l: i32 = l.try_into().unwrap();

    // Verify the range of the arguments k and m

    assert!(k.abs() <= l, "dlkm(): |k| > l");
    assert!(m.abs() <= l, "dlkm(): |m| > l");

    // Start computing the sum. Begin with determining the lower and upper boundaries.


    let lower = if -m-k > 0 { -m-k } else { 0 };
    let upper = if l-m > l-k { l-k } else { l-m };

    let cos_half_angle = f64::cos(angle/2.0);
    let sin_half_angle = f64::sin(angle/2.0);
    let mut sum: f64 = 0.0;
    for r in lower..=upper { 
        let term: f64 = binomial(l+m, l-k-r) * binomial(l-m, r) 
                * cos_half_angle.powf(f64::from(k+m+2*r)) * sin_half_angle.powf(f64::from(2*l-m-k-2*r));
        if (l-m-r).is_odd() {
            sum -= term;
        } else {
            sum += term;
        }
    }

    // Now multiply with the big square root factor. I use exp() because I have only ln(n!) available.

    sum *= (0.5 * (lnfac(l+k) + lnfac(l-k) - lnfac(l+m) - lnfac(l-m))).exp();

    return sum;
}