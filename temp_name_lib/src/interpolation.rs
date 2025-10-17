use std::num::ParseIntError;

use crate::utils::MathErrors;

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

    /// Intermediate Interpolation values. 
    pub partial_interpolations:Vec<f64>,

}

impl ParameterSpaceHypercube{
    /// Creates a new instance of a hypercube in parameter space
    fn new(dimension: usize)->Self{
        let mut fractional_coordinates:Vec<f64> = Vec::new();
        let mut fractional_distances:Vec<f64> = Vec::new();
        let mut corner_values:Vec<f64> =  Vec::new();

        for i in 0usize..2usize*dimension {
            fractional_coordinates.push(0.0);
        }
        for i in 0usize..dimension{
            fractional_distances.push(0.0);
        }
        for i in 0usize..2usize.pow(dimension){
            corner_values.push(0.0);
        }
        for i in 0usize..(2.pow(dimension)-1){
            partial_interpolations.push(0.0);
        }

    }

    /// Fills in the coordinates of the hypercube. 
    fn fill_coordinates(&mut self,pairs:&[f64])->Result<(),MathErrors>{
        if pairs.len()!= self.fractional_coordinates.len(){Err(MathErrors::NotAdequateNumberOfElements)}
        else{
            for (index,item) in pairs.iter().enumerate(){
                self.fractional_coordinates[index]=*item;
            }
            Ok(())
        }
    }

    fn get_fractional_distances( &mut self, coords_in_param_space:&[f64])->Result<(),MathErrors>{
       if coords_in_param_space.len()!= self.fractional_distances{Err(MathErrors::NotAdequateNumberOfElements)} 
       else{
            for (index,item) in self.fractional_coordinates.chunks(2usize).enumerate(){
                let x_l = item[0];
                let x_r = item[1];
                self.fractional_distances[index] = (coords_in_param_space[index]-x_l)/(x_r-x_l)
            }
       }
    }

    fn fill_vertices_data (&mut self, vertices_data:&[f64])->Result<(),MathErros>{
        if vertices_data.len() != self.corner_values.len(){Err(MathErrors::NotAdequateNumberOfElements)}
        else {
            for (index,item) in vertices_data.iter().enumerate(){
                self.corner_values[index] = *item;
            }
            Ok(())
        }
    }


    // so this is what I want to do. I will take a slice of the partial interpolation and fill it, then I will take the next slice
    // Since this is done in a loop, it should go like this
    // start_index=2^i-1?
    // end_index = 2^i?
    // slice_i = 2[start_index..=end_index]
    // and then store the intermediate interpolations. 
    // With this I would be able to iterate. 
    pub fn multilinear_interpolation(&mut self, coords_in_param_space:&[f64])->Result<_,f64>{
        self.get_fractional_distances(coords_in_param_space)?;        
        let mut partial_interpolation = self.corner_values;

        for i in 1usize..coords_in_param_space.len() {
            let dimension_counter=coords_in_param_space.len()+1-i;
            partial_interpolation = multilinear_interpolation(dimension_counter, partial_interpolation, self.fractional_distances[i])
        }
        
    }


}


pub fn multilinear_interpolation(dimension:usize,partial_interpolation:&[f64],fractional_distances:&f64)->Vec<f64>{
    let mut new_partial_interpolation:Vec<f64> = Vec::with_capacity(2usize.pow(dimension-1));
    for (index,pair) in partial_interpolation.chunks(2usize).enumerate(){
        new_partial_interpolation.push(linear_interpolation(pair, fractional_distances))
    }
}

pub fn linear_interpolation(values:&[f64],fractional_distance:f64)->f64{
    values[0]*fractional_distance +  values[1]*(1.0-fractional_distance)
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