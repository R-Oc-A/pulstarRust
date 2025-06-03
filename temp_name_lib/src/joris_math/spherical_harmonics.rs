pub mod plmcos;
pub mod d_plmcos_dtheta;
pub mod norm_factor;
pub mod dlkm_function;



































/*
#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
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
    #[should_panic]
    fn test_plmcos_panic() {
        let theta: f64 = 2.6;
        let (sintheta, costheta) = (theta.sin(), theta.cos());
        plmcos(4, 5, sintheta, costheta);
    }
}

*/