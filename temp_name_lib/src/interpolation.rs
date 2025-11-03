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

    /// Intermediate Interpolation values. The multilinear-interpolation will be done on an iterative process on the dimensions.  
    pub partial_interpolations:Vec<f64>,

}

impl ParameterSpaceHypercube{
    /// Creates a new instance of a hypercube in parameter space. This method is intendet to be used to create a mutable instance. 
    /// ### Arguments:
    /// `dimension` - A [usize] value indicating the number of dimensions on the parameter space. 
    /// ### Returns: 
    /// -A new instance of the [ParameterSpaceHypercube] filled with zeroes.
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
    /// ### Arguments: 
    /// * `pairs` - A &[[f64];2] reference containing the coordinates in parameter space ordered from lower to upper. 
    /// ### Returns: 
    /// * This method returns a [Result] with [Ok()] variant or an [[Err] ([MathErrors])] variant indicating that the number of pairs does not match the number of dimensions in parameter space. 
    pub fn fill_coordinates(&mut self,pairs:&[[f64;2]])->Result<(),MathErrors>{
        if pairs.len()!= self.fractional_coordinates.len(){Err(MathErrors::NotAdequateNumberOfElements)}
        else{
            for (index,item) in pairs.iter().enumerate(){
                self.fractional_coordinates[index]=*item;
            }
            Ok(())
        }
    }

    /// This method is used to compute the fractional distances to a point in parameter space.
    /// ### Arguments: 
    /// * `coords_in_param_space` - a &[[f64]] slice that contains the coordinates of a poin inside the [ParameterSpaceHypercube] where we want to know the result of the interpolation.
    /// ### Returns: 
    /// * This method returns a [Result] with a [Ok()] variant in case the fractional distances where calculated correctly and an 
    /// [Err] variant  in case the coordinates of the point are not well indicated or the point is outside the domain of the [ParameterSpaceHypercube].
    fn get_fractional_distances( &mut self, coords_in_param_space:&[f64])->Result<(),MathErrors>{
       if coords_in_param_space.len()!= self.fractional_distances.len(){
        println!("here're the values {},{}",coords_in_param_space.len(),self.fractional_distances.len());
        Err(MathErrors::NotAdequateNumberOfElements)} 
       else{
            for (index,item) in self.fractional_coordinates.iter().enumerate(){
                let x_l = item[0];
                let x_r = item[1];
                if coords_in_param_space[index]<x_l || coords_in_param_space[index]>x_r{
                    println!("One of the coordinates of the point in parameter space is out of bounds");
                    println!("coordinate value = {}, left bound = {}, right bound = {}",coords_in_param_space[index],x_l,x_r);
                    return Err(MathErrors::OutOfBounds)
                };
                
                self.fractional_distances[index] = (coords_in_param_space[index]-x_l)/(x_r-x_l);
            }
            Ok(())
       }
    }

    /// This function fills the data contained in the vertices of the [ParameterSpaceHypercube].
    /// ### Arguments:
    /// * `vertices_data` - a &[[f64]] reference that contains the values in the vertices. 
    /// ### Returns:
    /// * This method returns a [Result] with an [Ok()] variant in case the vertices in the [ParameterSpaceHypercube] were filled correctly and an
    /// [Err] ([MathErrors]) variant in case the number of elements provided to fill the data was not adequate. 
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

    ///This method performs multilinear interpolation for a point in the parameter space. The implementation was based on the algorithm described on the [Numerical Recipes](https://numerical.recipes/) book. 
    /// ### Arguments:
    /// * `Coords_in_param_space`- A &[[f64]] reference that contains the coordinates of a point in the parameter space. 
    /// ### Returns: 
    /// This function returns a [Result] with the following variants: 
    /// * [Ok] ([f64]) - where the binded value is the result of the interpolation. 
    /// * [Err] ([MathErrors]) - If there was a problem with the slice containing the coordinates on the parameter space. 
    pub fn multilinear_interpolation(&mut self, coords_in_param_space:&[f64])->Result<f64,MathErrors>{
        //Compute the fractional distances for all of the dimensions. 
        self.get_fractional_distances(coords_in_param_space)?;        
        
        // get the dimensions of the parameterspace
        let dimension = coords_in_param_space.len();

        // The first iteration of linear interpolations uses all of the data contained on the vertices of the hypercube. So first I produce a copy of the values contained there to perform the loop.
        let end_index = 2usize.pow(dimension as u32);
        let (slice0,_slice1)=self.partial_interpolations.split_at_mut(end_index);
        slice0.copy_from_slice(& self.corner_values[0..end_index]);
        //for index in 0..2usize.pow(dimension as u32){
        //    self.partial_interpolations[index] = self.corner_values[index];
        //}
        
        //Something like this but I still need to think on some(most) details 
        let mut start_index = 0usize;
        for dimension_step in 1..=dimension{

            let dimension_counter = dimension - dimension_step + 1 ;//should go from dimension to 1 in steps by 1
            let end_index = start_index+2usize.pow(dimension_counter as u32);
            
            let (old_interpolations,current_interpolation) = self.partial_interpolations.split_at_mut(end_index);
            let (_oldest_interpolations,previous_interpolation) = old_interpolations.split_at_mut(start_index);

            println!("start_index {}, end index {}",start_index,end_index);
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
                    [0.00100261, 0.00115899, 0.00124444, 0.0013054 , 0.00135344, 0.00139335, 0.00142767, 0.00099939, 0.00115426, 0.00123879, 0.00129903, 0.00134648, 0.00138588, 0.00141975],
                    [0.1, 0.00115898, 0.00124444, 0.00130539, 0.00135343, 0.00139335, 0.00142767, 0.00099949, 0.00115441, 0.00123897, 0.00129923, 0.00134669, 0.00138611, 0.00141998] 
                    //2DMatrix
                ],//First point on paramspace Teff=21000,logg=3.5,
                [
                    [0.0009855 , 0.00112527, 0.00120154, 0.00125596, 0.00129886, 0.00133452, 0.00136519, 0.00097379, 0.00110768, 0.00118042, 0.00123221, 0.00127298, 0.00130682, 0.00133589],
                    [0.0009855 , 0.00112527, 0.00120153, 0.00125595, 0.00129886, 0.00133452, 0.00136518, 0.00097372, 0.00110766, 0.00118043, 0.00123225, 0.00127303, 0.00130689, 0.00133597]
                ],//Second point on paramspace Teff=21000, logg=4.5,
                [
                    [0.00133685, 0.0015275 , 0.00163161, 0.00170613, 0.00176507, 0.00181421, 0.00185659, 0.00133438, 0.00152362, 0.00162682, 0.00170064, 0.001759  , 0.00180763, 0.00184957],
                    [0.00133685, 0.0015275 , 0.0016316 , 0.00170612, 0.00176507, 0.00181421, 0.00185659, 0.0013345 , 0.0015238 , 0.00162702, 0.00170086, 0.00175924, 0.00180789, 0.00184983]
                ],//Third point on paramspace Teff=24000, logg=3.5,
                [
                    [0.00131438, 0.00148353, 0.00157459, 0.00163924, 0.00169009, 0.0017323 , 0.00176859, 0.00130271, 0.00146613, 0.00155372, 0.00161575, 0.00166444, 0.0017048 , 0.00173945],
                    [0.00131438, 0.00148352, 0.00157458, 0.00163924, 0.00169009, 0.0017323 , 0.00176858, 0.00130296, 0.00146649, 0.00155414, 0.00161621, 0.00166493, 0.00170532, 0.00174   ]
                ]//Fourth point on paramspace Teff = 24000, logg=4.5
            ];
            SampleGrids { mu_values:mu_values, log_g:log_g, t_eff:t_eff, grid_values: arr, wavelengths:wavelengths }
        }
        
        fn fill_hypercube()->ParameterSpaceHypercube{
            let sample_grid = Self::nadya_sample();
            let mut coordinates_in_parameter_space:Vec<[f64;2]> = Vec::new();

            let find_mu_index = |x:f64,mu_vals:&[f64]|->usize{
                let mut index:usize =0;
                for (n,mu_val) in mu_vals.iter().enumerate(){
                    if x<=*mu_val { index = n;
                    break;}}
                index-1
            };
            coordinates_in_parameter_space.push(sample_grid.t_eff.clone());
            coordinates_in_parameter_space.push(sample_grid.log_g.clone());
            
            let mut slice:[f64;2] = [0.0;2];
            
            let wavelength:Vec<f64>=vec![4000.0,4000.1];
            let index_wavelengths = 0usize;
            slice.copy_from_slice(&wavelength[index_wavelengths..=index_wavelengths+1]);
            coordinates_in_parameter_space.push(slice);
            
            let mu_val = 0.66;
            let index = find_mu_index(mu_val,&sample_grid.mu_values);
            slice.copy_from_slice(&sample_grid.mu_values[index..=index+1]);
            coordinates_in_parameter_space.push(slice);
            let mut corner_values:Vec<f64> = Vec::with_capacity(16usize);
    
            for i in 0..2usize{//teff
                for j in 0..2usize{//logg
                    let grid_number = 2*i+j;
                    for k in 0..2usize{//lambda
                        for l in 0..2usize{//mu
                            corner_values.push(sample_grid.grid_values[[grid_number,index_wavelengths+k,index+l]]);
                        }
                    }
                }
            }

            ParameterSpaceHypercube { fractional_coordinates: coordinates_in_parameter_space, fractional_distances: vec![0.0;4], corner_values:corner_values, partial_interpolations: vec![0.0;2usize.pow(5)] }
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
    #[test]
    fn are_corner_values_filled_appropriately(){
        let hypercube = SampleGrids::fill_hypercube();
        let ex_sample = SampleGrids::nadya_sample();
        let arr = ex_sample.grid_values;
        let corner_values_manually_chosen = vec![
            arr[[0,0,2]],//Teff=21000, Logg = 3.5,wavelenght=4000.00,mu =0.5976
            arr[[0,0,3]],//Teff=21000, Logg = 3.5,wavelength=4000.00,mu = 0.7071
            arr[[0,1,2]],//Teff=21000, Logg = 3.5,wavelenght=4000.01,mu =0.5976
            arr[[0,1,3]],//Teff=21000, Logg = 3.5,wavelength=4000.01,mu = 0.7071
            arr[[1,0,2]],//Teff=21000, Logg = 4.5,wavelenght=4000.00,mu =0.5976
            arr[[1,0,3]],//Teff=21000, Logg = 4.5,wavelength=4000.00,mu = 0.7071
            arr[[1,1,2]],//Teff=21000, Logg = 4.5,wavelenght=4000.01,mu =0.5976
            arr[[1,1,3]],//Teff=24000, Logg = 4.5,wavelength=4000.01,mu = 0.7071
            arr[[2,0,2]],//Teff=24000, Logg = 3.5,wavelenght=4000.00,mu =0.5976
            arr[[2,0,3]],//Teff=24000, Logg = 3.5,wavelength=4000.00,mu = 0.7071
            arr[[2,1,2]],//Teff=24000, Logg = 3.5,wavelenght=4000.01,mu =0.5976
            arr[[2,1,3]],//Teff=24000, Logg = 3.5,wavelength=4000.01,mu = 0.7071
            arr[[3,0,2]],//Teff=24000, Logg = 3.5,wavelenght=4000.00,mu =0.5976
            arr[[3,0,3]],//Teff=24000, Logg = 3.5,wavelength=4000.00,mu = 0.7071
            arr[[3,1,2]],//Teff=24000, Logg = 3.5,wavelenght=4000.01,mu =0.5976
            arr[[3,1,3]],//Teff=24000, Logg = 3.5,wavelength=4000.01,mu = 0.7071
        ];
        let corner_values_if_done_right = vec![
            0.00124444,
            0.0013054,
            0.00124444,
            0.00130539,
            0.00120154,
            0.00125596,
            0.00120153,
            0.00125595,
            0.00163161,
            0.00170613,
            0.0016316,
            0.00170612,
            0.00157459,
            0.00163924,
            0.00157458,
            0.00163924
        ];
        assert_eq!(corner_values_manually_chosen,corner_values_if_done_right);
        assert_eq!(hypercube.corner_values,corner_values_if_done_right);
    }

    fn manual_grid_interpolation()->f64{

        let coordinates:Vec<f64> = vec![22000.0,4.32,4000.03,0.66];

        let mut hypercube = SampleGrids::fill_hypercube();

        hypercube.get_fractional_distances(&coordinates).unwrap();

        //first 16 partial lineal interpolations; this are done on temperature.
        let mut d = hypercube.fractional_distances[0];
        let mut c1: Vec<f64> = Vec::with_capacity(8usize);
        for chunk in hypercube.corner_values.chunks(2usize){
            c1.push(linear_interpolation(chunk, d));
        }
        if c1.len()!= 8 {panic!("not appropriate size, step 1")};
        
        let mut c2:Vec<f64>  = Vec::with_capacity(4);
        d = hypercube.fractional_distances[1];
        for chunk in c1.chunks(2usize){
            c2.push(linear_interpolation(chunk, d));
        }
        if c2.len() != 4 {panic!("not appropriate size, step 2")};

        d = hypercube.fractional_distances[2];
        let mut c3:Vec<f64> = Vec::with_capacity(2usize);
        for chunk in c2.chunks(2usize){
            c3.push(linear_interpolation(chunk,d))
        }
        if c3.len() != 2 {panic!("not appropriate size, step 3")};

        d = hypercube.fractional_distances[3];
        let c4 = linear_interpolation(&c3, d);

        c4
    }
    #[test]
    
    fn multilinear_interpolation_test(){
        let coordinates:Vec<f64> = vec![22000.0,4.32,4000.03,0.66];

        let mut hypercube = SampleGrids::fill_hypercube();

        assert_eq!(hypercube.multilinear_interpolation(&coordinates).unwrap(),manual_grid_interpolation())

    }    
}





