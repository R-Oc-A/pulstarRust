use crate::utils::MathErrors;



pub trait LinearlyInterpolatable{
    fn linear_interpolation(left:&Self,right:&Self,fractional_distance:f64)->Self;
}

/// Here is the implementation to perform multilinear interpolation on an [f64] variable. 
pub mod float64;
/// Here is the implementation to perform multilinear interpolation on a [polars::LazyFrame], where all of its columns contain [f64] elements
pub mod polars_lazyframe;


/// This structure contains all the relevant information to produce 
/// multilinear interpolation. 
/// 
pub struct ParameterSpaceHypercube<T>
where
    T:LinearlyInterpolatable+Clone,
{
    /// This vector contains all the pairs of values that define the hypercube in the parameter space.
    pub fractional_coordinates:Vec<[f64;2]>,

    /// This vector holds the pairs fractional distances that could be computed when requesting for an interpolation. 
    /// 
    /// We'll be using by convention the following:
    /// Let (x1,x2,...,xN) be the coordinates of the point on parameter space where one wants to interpolate. Then
    /// the fractional_distance is  (xi-Xi_L)/(Xi_R - Xi_L) where X_I,X_I+1 are the position of the vertex of the hypercube along the i-direction.
    pub fractional_distances:Vec<f64>,

    /// This [Vec] must contain the 2^N values of the vertices of the hypercube on the parameter space. 
    pub corner_values:Vec<T>,

    /// Intermediate Interpolation values. The multilinear-interpolation will be done on an iterative process on the dimensions.  
    pub partial_interpolations:Vec<T>,
}

impl <T> ParameterSpaceHypercube<T>
where
    T:LinearlyInterpolatable + Clone
{
    /// Creates a new instance of a hypercube in parameter space. This method is intendet to be used to create a mutable instance. 
    /// ### Arguments:
    /// `dimension` - A [usize] value indicating the number of dimensions on the parameter space. 
    /// ### Returns: 
    /// -A new instance of the [ParameterSpaceHypercube] filled with zeroes.
    pub fn new(dimension: usize)->Self{
        let mut fractional_coordinates:Vec<[f64;2]> = Vec::new();
        let mut fractional_distances:Vec<f64> = Vec::new();

        for _i in 0usize..dimension {
            fractional_coordinates.push([0.0;2]);
            fractional_distances.push(0.0);
        }

        let corner_values:Vec<T> =  Vec::with_capacity(2usize.pow(dimension as u32));
        let partial_interpolations:Vec<T> = Vec::with_capacity(2usize.pow(dimension as u32+1));

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
    /// * `vertices_data` - a &[T] reference that contains the values in the vertices. 
    /// ### Returns:
    /// * This method returns a [Result] with an [Ok()] variant in case the vertices in the [ParameterSpaceHypercube] were filled correctly and an
    /// [Err] ([MathErrors]) variant in case the number of elements provided to fill the data was not adequate. 
    pub fn fill_vertices_data (&mut self, vertices_data:&[T])->Result<(),MathErrors>{
        if ! self.corner_values.is_empty(){
            if vertices_data.len() != self.corner_values.len(){Err(MathErrors::NotAdequateNumberOfElements)}
            else {
                for (index,item) in vertices_data.iter().enumerate(){
                    self.corner_values[index] = item.clone();
                }
                Ok(())
            }
        }else{
            let dimension = self.fractional_coordinates.len();
            if vertices_data.len() != 2usize.pow(dimension as u32){Err(MathErrors::NotAdequateNumberOfElements)}
            else{
                for item in vertices_data.iter(){
                    self.corner_values.push(item.clone());
                }
                for _index in 0..2usize.pow(dimension as u32 + 1){
                    self.partial_interpolations.push(self.corner_values[0].clone());
                }
                Ok(())
            }
        }
    }


    // So this is what I want to do. I will take a slice of the partial interpolation and fill it, then I will take the next slice
    // Since this is done in a loop, it should go like this
    // start_index=2^i-1?
    // end_index = 2^i?
    // slice_i = 2[start_index..=end_index]
    // and then store the intermediate interpolations. 
    // With this I would be able to iterate. 

    /// This method performs multilinear interpolation for a point in the parameter space. The implementation was based on the algorithm described on the [Numerical Recipes](https://numerical.recipes/) book. 
    /// The implementation thus relies on doing linear interpolations in a loop on the dimensions.
    /// 
    /// For two dimensions, if the corner values are `C00,C01,C10,C11`, and the fractional distances are `t=(x-xl)/(xr-xl)` and `u =(y-yl)/(yr-yl)`, we have
    /// 
    /// `C1 = C00*t + C01 *(1-t); //<- Linear interpolation`
    /// 
    /// `C2 = C10*t + C11 *(1-t);`
    /// 
    /// `multi_interpolation = C1*u + C2*(1-u);`
    /// 
    /// ### Arguments:
    /// * `Coords_in_param_space`- A &[[f64]] reference that contains the coordinates of a point in the parameter space. 
    /// ### Returns: 
    /// This function returns a [Result] with the following variants: 
    /// * [Ok] ([T]) - where the binded value is the result of the interpolation. 
    /// * [Err] ([MathErrors]) - If there was a problem with the slice containing the coordinates on the parameter space. 
    pub fn multilinear_interpolation(&mut self, coords_in_param_space:&[f64])->Result<T,MathErrors>{
        //Compute the fractional distances for all of the dimensions. 
        self.get_fractional_distances(coords_in_param_space)?;        
        
        // get the dimensions of the parameterspace
        let dimension = coords_in_param_space.len();

        // The first iteration of linear interpolations uses all of the data contained on the vertices of the hypercube. So first I produce a copy of the values contained there to perform the loop.
        // We first copy the vertices values into the partial_interpolations vector. 
        let mut end_index = 2usize.pow(dimension as u32);
        let (slice0,_slice1)=self.partial_interpolations.split_at_mut(end_index);//Here we disregard slice1 because it will contain the coefficient of the following linear interpolations. 
        slice0.clone_from_slice(& self.corner_values[0..end_index]);
        
        //Begin the loop on dimensions:
        let mut start_index = 0usize;
        for dimension_step in 1..=dimension{

            let dimension_counter = dimension - dimension_step + 1 ;//should go from dimension to 1 in steps by 1
            end_index = start_index+2usize.pow(dimension_counter as u32);
            
            let (old_interpolations,current_interpolation) = self.partial_interpolations.split_at_mut(end_index);
            let (_oldest_interpolations,previous_interpolation) = old_interpolations.split_at_mut(start_index);

            //println!("start_index {}, end index {}",start_index,end_index);
            for (index,pair) in previous_interpolation.chunks(2usize).enumerate(){
                current_interpolation[index]=T::linear_interpolation(&pair[0],&pair[1], self.fractional_distances[dimension_step-1]);
            }
            start_index = end_index;
        }//end loop
        //println!("end index{}",end_index);
        Ok(self.partial_interpolations[end_index].clone())
    }
}
