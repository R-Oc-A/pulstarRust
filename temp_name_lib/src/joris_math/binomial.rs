use super::lnfac::lnfac;

///  Compute the binomial number (n k) = n!/(k! (n-k)!) for 0 <= k <= n.
///
pub fn binomial(n: i32, k: i32) -> f64 {

    // Do a range check for the arguments

    assert!(n >= 0, "binomial(): n < 0");
    assert!(k >= 0, "binomial(): k < 0");
    assert!(k <= n, "binomial(): k > n");

    // The array BINOMIALTABLE[0..10][0..10] contains the binomial number 
    // (n k) = n!/(k! (n-k)!), for 0 <= k <= n <= 10. 
    // E.g. binomial[5][3] = 5!/(3! (5-3)!) = 10.

    const BINOMIALTABLE: [[f64; 11]; 11] = 
      [ [1.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0],
        [1.0,  1.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0],
        [1.0,  2.0,  1.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0],
        [1.0,  3.0,  3.0,   1.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0],
        [1.0,  4.0,  6.0,   4.0,   1.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0],
        [1.0,  5.0, 10.0,  10.0,   5.0,   1.0,   0.0,   0.0,  0.0,  0.0, 0.0],
        [1.0,  6.0, 15.0,  20.0,  15.0,   6.0,   1.0,   0.0,  0.0,  0.0, 0.0],
        [1.0,  7.0, 21.0,  35.0,  35.0,  21.0,   7.0,   1.0,  0.0,  0.0, 0.0],
        [1.0,  8.0, 28.0,  56.0,  70.0,  56.0,  28.0,   8.0,  1.0,  0.0, 0.0],
        [1.0,  9.0, 36.0,  84.0, 126.0, 126.0,  84.0,  36.0,  9.0,  1.0, 0.0],
        [1.0, 10.0, 45.0, 120.0, 210.0, 252.0, 210.0, 120.0, 45.0, 10.0, 1.0] ];

    // For small numbers n, use the pre-tabulated value, otherwise compute
    // manually

    if n <= 10 {
        return BINOMIALTABLE[n as usize][k as usize];
    } else {
        // The individual factorials in a binomial overflow long before the coefficient
        // itself will. So work with logarithms. The floor function cleans roundoff error 
        // for smaller values of n and k.

        return (0.5f64 + (lnfac(n) - lnfac(k) - lnfac(n-k)).exp()).floor();
    }
} 
