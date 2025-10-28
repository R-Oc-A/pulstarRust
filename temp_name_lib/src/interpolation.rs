use crate::utils::MathErrors;

/// This structure contains all the relevant information to produce 
/// multilinear interpolation. 
/// 
pub struct ParameterSpaceHypercube{
    /// This vector contains all the pairs of values that define the hypercube in the parameter space.
    pub fractional_coordinates:Vec<[f64;2]>,

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
    pub fn new(dimension: usize)->Self{
        let mut fractional_coordinates:Vec<[f64;2]> = Vec::new();
        let mut fractional_distances:Vec<f64> = Vec::new();
        let mut corner_values:Vec<f64> =  Vec::new();
        let mut partial_interpolations:Vec<f64> = Vec::new();

        for _i in 0usize..dimension {
            fractional_coordinates.push([0.0;2]);
            fractional_distances.push(0.0);
        }
        for _i in 0usize..2usize.pow(dimension as u32){
            corner_values.push(0.0);
        }
        for _i in 0usize..(2usize.pow(dimension as u32 +1)){
            partial_interpolations.push(0.0);
        }

        Self { fractional_coordinates: fractional_coordinates,
               fractional_distances: fractional_distances,
               corner_values: corner_values,
               partial_interpolations: partial_interpolations }

    }

    /// Fills in the coordinates of the hypercube. 
    pub fn fill_coordinates(&mut self,pairs:&[[f64;2]])->Result<(),MathErrors>{
        if pairs.len()!= self.fractional_coordinates.len(){Err(MathErrors::NotAdequateNumberOfElements)}
        else{
            for (index,item) in pairs.iter().enumerate(){
                self.fractional_coordinates[index]=*item;
            }
            Ok(())
        }
    }

    fn get_fractional_distances( &mut self, coords_in_param_space:&[f64])->Result<(),MathErrors>{
       if coords_in_param_space.len()!= self.fractional_distances.len(){Err(MathErrors::NotAdequateNumberOfElements)} 
       else{
            for (index,item) in self.fractional_coordinates.iter().enumerate(){
                let x_l = item[0];
                let x_r = item[1];
                self.fractional_distances[index] = (coords_in_param_space[index]-x_l)/(x_r-x_l)
            }
            Ok(())
       }
    }

    pub fn fill_vertices_data (&mut self, vertices_data:&[f64])->Result<(),MathErrors>{
        if vertices_data.len() != self.corner_values.len(){Err(MathErrors::NotAdequateNumberOfElements)}
        else {
            for (index,item) in vertices_data.iter().enumerate(){
                self.corner_values[index] = *item;
            }
            Ok(())
        }
    }


    // So this is what I want to do. I will take a slice of the partial interpolation and fill it, then I will take the next slice
    // Since this is done in a loop, it should go like this
    // start_index=2^i-1?
    // end_index = 2^i?
    // slice_i = 2[start_index..=end_index]
    // and then store the intermediate interpolations. 
    // With this I would be able to iterate. 
    pub fn multilinear_interpolation(&mut self, coords_in_param_space:&[f64])->Result<f64,MathErrors>{
        self.get_fractional_distances(coords_in_param_space)?;        

        let dimension = coords_in_param_space.len();
        for index in 0..2usize.pow(dimension as u32){
            self.partial_interpolations[index] = self.corner_values[index];
        }

        //Something like this but I still need to think on some(most) details 
        let mut start_index = 0usize;
        for dimension_step in 1..=dimension{
            let dimension_counter = dimension - dimension_step + 1 ;//should go from dimension to 1 in steps by 1
            let end_index = start_index+2usize.pow(dimension_counter as u32);
            
            let (previous_interpolation,current_interpolation) = self.partial_interpolations.split_at_mut(end_index);
            for (index,pair) in previous_interpolation.chunks(2usize).enumerate(){
                current_interpolation[index]=linear_interpolation(pair, self.fractional_distances[dimension_step-1]);
            }
            start_index = end_index;//???not sure
        }

        Ok(self.partial_interpolations[start_index])
    }
}

fn linear_interpolation(values:&[f64],fractional_distance:f64)->f64{
    values[0]*fractional_distance +  values[1]*(1.0-fractional_distance)
}



#[cfg(test)]
mod tests {
    use ndarray::{Array3, array};

    use super::*;
    //quizas valdria la pena  pensar en utilizar un diccionario  para estas abstracciones
    struct SampleGrids{
        mu_values:[f64;7],
        log_g:[f64;2],
        t_eff:[f64;2],
        grid_values:ndarray::Array3<f64>,
        wavelengths: Vec<f64>
    }

    impl SampleGrids{
        ///This interpolation methods are developed to work using a grid of specific intensities calculated by Nadya. 
        ///Thus here we construct a sample grid of 4 wavelengths, 2 mu values, and 2 values of log_gravity and effective temperature. 
        fn nadya_sample()->Self {

            let mu_values = [0.2673, 0.4629, 0.5976, 0.7071, 0.8018, 0.8864, 0.9636];//[0.9636, 0.8864, 0.8018, 0.7071, 0.5976, 0.4629, 0.2673];
            let t_eff = [21000.0,24000.0];
            let log_g= [3.5,4.5];

            let wavelengths = vec![400.0,400.1];

            let arr:Array3<f64> = array![
                [
                    
                    //2DMatrix
                ],//Fourth point on paramspace Teff=24000,logg=4.5,
            ];

            SampleGrids { mu_values:mu_values, log_g:log_g, t_eff:t_eff, grid_values: arr, wavelengths:wavelengths }
        }
        
        fn fill_hypercube(t_eff:f64, log_g:f64, mu_val:f64, wavelength:f64)->ParameterSpaceHypercube{
            let sample_grid = Self::nadya_sample();
            let mut coordinates_in_parameter_space:Vec<[f64;2]> = Vec::new();

            let find_mu_index = |x:f64,mu_vals:&[f64]|->usize{
                let index:usize =0;
                for (n,mu_val) in mu_vals.iter().enumerate(){
                    if x>=mu_val || index < 6 { index = n}
                }
                index
            };
            coordinates_in_parameter_space.push(sample_grid.log_g.clone());
            coordinates_in_parameter_space.push(sample_grid.t_eff.clone());
            let index = find_mu_index(mu_val,sample_grid.mu_values);
            
            coordinates_in_parameter_space.push(sample_grid.mu_values[index..=index+1].clone());
        }
    }
    //This is a test to see if the linear interpolation works well for  a point on a line. 
    #[test]
    fn interpolate_traveltime(){
        let point_a=1.0;
        let point_c = 7.0;
        let travel_time_a = 5.0;
        let travel_time_c = 18.0;

        let point_b = 0.5*(point_a + point_c)//point b is the midpoint between a and c.
        ;

        let travel_time_b = travel_time_a + (point_b-point_a) * (travel_time_c-travel_time_a)/(point_c-point_a);

        let values= vec![travel_time_a,travel_time_c];
        let fractional_distance = (point_b-point_a )/(point_c-point_a);

        assert_eq!(travel_time_b,linear_interpolation(&values, fractional_distance))
    }

    // Now working with a grid of data

    

    
    
}





