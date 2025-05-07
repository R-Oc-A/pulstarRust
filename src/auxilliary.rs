
/// Compute ln(n!)
///
/// This code is a combination and a variant of the functions 
/// gammln() and factrl() in Numerical Recipes in C., Press et al., 1992
///
pub fn lnfac(n: i32) -> f64{

    assert!(n >= 0, "lnfac(): n < 0");

    // The array lnfactable[0..20] contains ln(n!)
    // E.g. lnfactable[5] = ln(5!) = 4.78749174278.

    const LNFACTABLE: [f64; 21] =
      [ 0.0, 0.0, 0.69314718056, 1.79175946923, 3.17805383035,
        4.78749174278, 6.57925121201, 8.52516136107, 10.6046029027,
        12.8018274801, 15.1044125731, 17.5023078459, 19.9872144957,
        22.5521638531, 25.1912211827, 27.8992713838, 30.6718601061,
        33.5050734501, 36.395445208,  39.3398841872, 42.3356164608
      ];

    // The array cof[0..5] contains info to compute the gamma function

    const COF: [f64; 6] = [ 76.18009172947146, -86.50532032941677, 24.01409824083091,
                            -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 ];

    // For small values of n use the pre-computed table. For large values, compute the result manually.

    if n <= 20 {
        return LNFACTABLE[n as usize];
    } else { 
        let x: f64 = f64::from(n) + 1.0;
        let mut y = x;
        let mut temp: f64 = x + 5.5;
        temp -= (x + 0.5) * temp.ln();
        let mut ser: f64 = 1.000000000190015;
        for j in 0..=5 { 
            y += 1.0;
            ser += COF[j]/y;
        }

        return -temp + (2.5066282746310005 * ser / x).ln();
    }
} 










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

        return (0.5 + (lnfac(n) - lnfac(k) - lnfac(n-k)).exp()).floor();
    }
} 

