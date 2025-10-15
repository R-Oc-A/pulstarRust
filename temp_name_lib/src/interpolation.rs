/// This structure contains all the relevant information to produce 
/// multilinear interpolation. 
/// 
pub struct ParameterSpaceHypercube{
    /// This vector contains all the pairs of values that define the hypercube in the parameter space.
    pub fractional_coordinates:Vec<f64>,

    /// This vector holds the pairs fractional distances that could be computed when requesting for an interpolation. 
    /// 
    /// We'll be using by convention the following:
    /// Let (x1,x2,...,xN) be the coordinates of the point on parameter space where one wants to interpolate. Then
    /// the fractional_distance is  (xi-Xi_L)/(Xi_R - Xi_L) where X_I,X_I+1 are the position of the vertex of the hypercube along the i-direction.
    pub fractional_distances:Vec<f64>,

    /// This [Vec] must contain the 2^N values of the vertices of the hypercube on the parameter space. 
    pub corner_values:Vec<f64>,

}

impl ParameterSpaceHypercube{
    /// Creates a new instance of a hypercube in parameter space
    fn new(dimension: usize)->Self{
        let mut fractional_coordinates:Vec<f64> = Vec::with_capacity(dimension);
        let mut fractional_distances:Vec<f64> = Vec::with_capacity(dimension);
        let mut corner_values:Vec<f64> =  Vec::with_capacity(2usize.pow(dimension));

        for i in 0usize..dimension {
            fractional_coordinates
        }

    }

    /// Fills the coordinates and values of the hypercube.
    fn refill

}

/// Multilinear interpolation
/// 
pub fn multilinear_interpolation(
    point_coordinates:&[f64],
    hypercube_corners:&[f64],
    vertex_on_param_space:&[f64],
)->f64{
 0.0;
}


fn compute_weights(
    point_coordinates:&[f64],
    hypercube_corners:&[f64]
)->Vec<f64>{
    let mut weight:Vec<f64> = Vec::with_capacity(capacity)