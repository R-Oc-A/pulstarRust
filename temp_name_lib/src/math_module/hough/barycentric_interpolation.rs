use scirs2_core::Array1;
use scirs2_interpolate::advanced::barycentric::BarycentricInterpolator;
use crate::type_def::PI;


/// Construct the gird points in the colatitude space
pub fn construct_mu_array_on_grid(npts:usize)->Array1<f64>{
    let step = PI/((npts -1) as f64);
    let mut mu_vec:Vec<f64> = Vec::with_capacity(npts);
    
    for i in 0..npts{
        let phase = step * i as f64;
        mu_vec.push(phase.cos())
    }
    scirs2_core::Array1::<f64>::from_vec(mu_vec)
}

fn construct_interpolator(mu_collocation_points:&Array1<f64>,hough_f:&Array1<f64>)->BarycentricInterpolator<f64>{
    BarycentricInterpolator::new(
        &mu_collocation_points.view(),&hough_f.view(),1usize
    ).expect("Error while creating barycentric interpolator")
}

/// interpolate hough functions into grid points.
pub fn interpolate_hough_function(mu:&Array1<f64>,mu_collocation_points:&Array1<f64>,hough_f:&Array1<f64>)->Array1<f64>{
    let interp =construct_interpolator(mu_collocation_points, hough_f);
    let hough_f_on_grid = interp.evaluate_array(&mu.view()).expect("error while interpolating hough function on grid points");
    hough_f_on_grid
}



