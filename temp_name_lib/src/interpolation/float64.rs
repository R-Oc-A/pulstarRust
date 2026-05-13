use super::*;
impl LinearlyInterpolatable for f64{
    fn linear_interpolation(left:&f64,right:&f64,fractional_distance:f64)->Self{
        left*fractional_distance +right*(1.0-fractional_distance) 
    }
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
        
        fn fill_hypercube()->ParameterSpaceHypercube<f64>{
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

        assert_eq!(travel_time_b,f64::linear_interpolation(&values[0],&values[1], fractional_distance))
    }

    // Now working with a grid of data
    #[test]
    fn corner_values_filled_appropriately(){
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
            arr[[1,1,3]],//Teff=21000, Logg = 4.5,wavelength=4000.01,mu = 0.7071
            arr[[2,0,2]],//Teff=24000, Logg = 3.5,wavelenght=4000.00,mu =0.5976
            arr[[2,0,3]],//Teff=24000, Logg = 3.5,wavelength=4000.00,mu = 0.7071
            arr[[2,1,2]],//Teff=24000, Logg = 3.5,wavelenght=4000.01,mu =0.5976
            arr[[2,1,3]],//Teff=24000, Logg = 3.5,wavelength=4000.01,mu = 0.7071
            arr[[3,0,2]],//Teff=24000, Logg = 4.5,wavelenght=4000.00,mu =0.5976
            arr[[3,0,3]],//Teff=24000, Logg = 4.5,wavelength=4000.00,mu = 0.7071
            arr[[3,1,2]],//Teff=24000, Logg = 4.5,wavelenght=4000.01,mu =0.5976
            arr[[3,1,3]],//Teff=24000, Logg = 4.5,wavelength=4000.01,mu = 0.7071
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
            c1.push(f64::linear_interpolation(&chunk[0],&chunk[1], d));
        }
        if c1.len()!= 8 {panic!("not appropriate size, step 1")};
        
        let mut c2:Vec<f64>  = Vec::with_capacity(4);
        d = hypercube.fractional_distances[1];
        for chunk in c1.chunks(2usize){
            c2.push(f64::linear_interpolation(&chunk[0],&chunk[1], d));
        }
        if c2.len() != 4 {panic!("not appropriate size, step 2")};

        d = hypercube.fractional_distances[2];
        let mut c3:Vec<f64> = Vec::with_capacity(2usize);
        for chunk in c2.chunks(2usize){
            c3.push(f64::linear_interpolation(&chunk[0],&chunk[1],d))
        }
        if c3.len() != 2 {panic!("not appropriate size, step 3")};

        d = hypercube.fractional_distances[3];
        let c4 = f64::linear_interpolation(&c3[0],&c3[1], d);

        c4
    }
    #[test]    
    fn multilinear_interpolation_works_on_grid_sample(){
        let coordinates:Vec<f64> = vec![22000.0,4.32,4000.03,0.66];

        let mut hypercube = SampleGrids::fill_hypercube();

        assert_eq!(hypercube.multilinear_interpolation(&coordinates).unwrap(),manual_grid_interpolation())

    }

    #[test]
    fn multilinear_interpolation_works_on_1d(){
        let point_a=1.0;
        let point_c = 7.0;
        let point_b = 0.5 * (point_a+point_c);
        let travel_time_a = 5.0;
        let travel_time_c = 18.0;

        let fractional_distance = (point_b-point_a )/(point_c-point_a);

        let coordinates:Vec<[f64;2]>=vec![[point_a,point_c]];
        let corner_values:Vec<f64>=vec![travel_time_a,travel_time_c];

        let mut hypercube = ParameterSpaceHypercube::new(1);

        hypercube.fractional_coordinates = coordinates;
        hypercube.fill_vertices_data(&corner_values).unwrap();
        //hypercube.get_fractional_distances(&[point_b]);

        assert_eq!(f64::linear_interpolation(&corner_values[0],&corner_values[1], fractional_distance),hypercube.multilinear_interpolation(&[point_b]).unwrap())
    }

}


