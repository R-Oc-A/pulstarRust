use polars::prelude::*;
use super::*;

trait Add{
    fn add(self,other:Self)->Self;
}

trait Mul{
    fn mul(self,rhs:f64)->Self;
}


impl Add for LazyFrame{
    fn add(self, other:Self) -> Self{
        let schema_left = self.clone().collect_schema().unwrap();
        let schema_right = other.clone().collect_schema().unwrap();
        if schema_left == schema_right {
            let field_names:Vec<String> = schema_left.iter_names().map(|x| x.to_string()).collect();
            let self_field_columns:Vec<Expr> = field_names.iter().map(|x|col(x)).collect();
            let mut other_field_columns:Vec<Expr> = field_names.iter().map(|x| col(x).alias(format!("{}_other",x))).collect();
            other_field_columns[0]= self_field_columns[0].clone();
            let other_to_join = other.clone().select(other_field_columns.clone());


            let lf = self.clone().join(other_to_join,
            [self_field_columns[0].clone()],
            [other_field_columns[0].clone()],
            JoinArgs::default());

            let mut expression_of_addition:Vec<Expr> = Vec::new();
            for (index, expression) in self_field_columns.into_iter().enumerate(){
                if index == 0usize{
                    //for wavelength
                    expression_of_addition.push(expression)
                }else{
                    expression_of_addition.push((expression + other_field_columns[index].clone()).alias(field_names[index].clone()));
                }
            };
            lf.select( expression_of_addition)
        }else{
            panic!("LazyFrames to be added have different schemas")
        }     
    }
}

impl Mul for LazyFrame{
    fn mul(self, rhs:f64)->Self{
        let schema = self.clone().collect_schema().unwrap();
        let field_names:Vec<String> = schema.iter_names().map(|x|x.to_string()).collect();
        let self_field_columns:Vec<Expr> = field_names.iter().map(|x| col(x)).collect();
        let mut expression_of_multiplication:Vec<Expr> = Vec::new();
        for (index,expr) in self_field_columns.into_iter().enumerate(){
            if index == 0{
                expression_of_multiplication.push(expr);
            }else{
                expression_of_multiplication.push(
                    expr * lit(rhs)
                )
            }
        }
        self.clone().select(expression_of_multiplication)
    }
}

impl LinearlyInterpolatable for LazyFrame{
    fn linear_interpolation(left:&Self,right:&Self,fractional_distance:f64)->Self {
        let lf = (left.clone().mul(fractional_distance))
        .add( right.clone().mul(1.0-fractional_distance));
        lf
    }
}

mod tests{
    use ndarray::{Array3, array, s};
    use super::*;
    struct SampleGrids{
        mu_values:[f64;7],
        log_g:[f64;2],
        t_eff:[f64;2],
        grid_values:Vec<LazyFrame>,
        wavelengths: Vec<f64>
    }

    impl SampleGrids{
        ///This interpolation methods are developed to work using a grid of specific intensities calculated by Nadya. 
        ///Thus here we construct a sample grid of 4 wavelengths, 2 mu values, and 2 values of log_gravity and effective temperature. 
        fn nadya_vec_df()->Vec<DataFrame>{
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
            let vec_df:Vec<DataFrame> = vec![
                df!(
                    "wavelength" => vec![4000.0,4000.1],
                    "mu1_flux" => arr.slice(s![0,0..=1,0]).to_vec(),
                    "mu2_flux" => arr.slice(s![0,0..=1,1]).to_vec(),
                    "mu3_flux" => arr.slice(s![0,0..=1,2]).to_vec(),
                    "mu4_flux" => arr.slice(s![0,0..=1,3]).to_vec(),
                    "mu5_flux" => arr.slice(s![0,0..=1,4]).to_vec(),
                    "mu6_flux" => arr.slice(s![0,0..=1,5]).to_vec(),
                    "mu7_flux" => arr.slice(s![0,0..=1,6]).to_vec(),
                    "mu1_cont" => arr.slice(s![0,0..=1,7]).to_vec(),
                    "mu2_cont" => arr.slice(s![0,0..=1,8]).to_vec(),
                    "mu3_cont" => arr.slice(s![0,0..=1,9]).to_vec(),
                    "mu4_cont" => arr.slice(s![0,0..=1,10]).to_vec(),
                    "mu5_cont" => arr.slice(s![0,0..=1,11]).to_vec(),
                    "mu6_cont" => arr.slice(s![0,0..=1,12]).to_vec(),
                    "mu7_cont" => arr.slice(s![0,0..=1,13]).to_vec(),
                ).unwrap(),
                df!(
                    "wavelength" => vec![4000.0,4000.1],
                    "mu1_flux" => arr.slice(s![1,0..=1,0]).to_vec(),
                    "mu2_flux" => arr.slice(s![1,0..=1,1]).to_vec(),
                    "mu3_flux" => arr.slice(s![1,0..=1,2]).to_vec(),
                    "mu4_flux" => arr.slice(s![1,0..=1,3]).to_vec(),
                    "mu5_flux" => arr.slice(s![1,0..=1,4]).to_vec(),
                    "mu6_flux" => arr.slice(s![1,0..=1,5]).to_vec(),
                    "mu7_flux" => arr.slice(s![1,0..=1,6]).to_vec(),
                    "mu1_cont" => arr.slice(s![1,0..=1,7]).to_vec(),
                    "mu2_cont" => arr.slice(s![1,0..=1,8]).to_vec(),
                    "mu3_cont" => arr.slice(s![1,0..=1,9]).to_vec(),
                    "mu4_cont" => arr.slice(s![1,0..=1,10]).to_vec(),
                    "mu5_cont" => arr.slice(s![1,0..=1,11]).to_vec(),
                    "mu6_cont" => arr.slice(s![1,0..=1,12]).to_vec(),
                    "mu7_cont" => arr.slice(s![1,0..=1,13]).to_vec(),
                ).unwrap(),
                df!(
                    "wavelength" => vec![4000.0,4000.1],
                    "mu1_flux" => arr.slice(s![2,0..=1,0]).to_vec(),
                    "mu2_flux" => arr.slice(s![2,0..=1,1]).to_vec(),
                    "mu3_flux" => arr.slice(s![2,0..=1,2]).to_vec(),
                    "mu4_flux" => arr.slice(s![2,0..=1,3]).to_vec(),
                    "mu5_flux" => arr.slice(s![2,0..=1,4]).to_vec(),
                    "mu6_flux" => arr.slice(s![2,0..=1,5]).to_vec(),
                    "mu7_flux" => arr.slice(s![2,0..=1,6]).to_vec(),
                    "mu1_cont" => arr.slice(s![2,0..=1,7]).to_vec(),
                    "mu2_cont" => arr.slice(s![2,0..=1,8]).to_vec(),
                    "mu3_cont" => arr.slice(s![2,0..=1,9]).to_vec(),
                    "mu4_cont" => arr.slice(s![2,0..=1,10]).to_vec(),
                    "mu5_cont" => arr.slice(s![2,0..=1,11]).to_vec(),
                    "mu6_cont" => arr.slice(s![2,0..=1,12]).to_vec(),
                    "mu7_cont" => arr.slice(s![2,0..=1,13]).to_vec(),
                ).unwrap(),
                df!(
                    "wavelength" => vec![4000.0,4000.1],
                    "mu1_flux" => arr.slice(s![3,0..=1,0]).to_vec(),
                    "mu2_flux" => arr.slice(s![3,0..=1,1]).to_vec(),
                    "mu3_flux" => arr.slice(s![3,0..=1,2]).to_vec(),
                    "mu4_flux" => arr.slice(s![3,0..=1,3]).to_vec(),
                    "mu5_flux" => arr.slice(s![3,0..=1,4]).to_vec(),
                    "mu6_flux" => arr.slice(s![3,0..=1,5]).to_vec(),
                    "mu7_flux" => arr.slice(s![3,0..=1,6]).to_vec(),
                    "mu1_cont" => arr.slice(s![3,0..=1,7]).to_vec(),
                    "mu2_cont" => arr.slice(s![3,0..=1,8]).to_vec(),
                    "mu3_cont" => arr.slice(s![3,0..=1,9]).to_vec(),
                    "mu4_cont" => arr.slice(s![3,0..=1,10]).to_vec(),
                    "mu5_cont" => arr.slice(s![3,0..=1,11]).to_vec(),
                    "mu6_cont" => arr.slice(s![3,0..=1,12]).to_vec(),
                    "mu7_cont" => arr.slice(s![3,0..=1,13]).to_vec(),
                ).unwrap()
            ];
            vec_df
        }
        fn nadya_sample(vec_lf:&[LazyFrame])->Self{

            let mu_values = [0.2673, 0.4629, 0.5976, 0.7071, 0.8018, 0.8864, 0.9636];//[0.9636, 0.8864, 0.8018, 0.7071, 0.5976, 0.4629, 0.2673];
            let t_eff = [21000.0,24000.0];
            let log_g= [3.5,4.5];

            let wavelengths = vec![4000.0,4000.1];
            //let vec_lf = vec_df.iter().map(|x| {x.lazy().clone()}).collect();

            SampleGrids {mu_values:mu_values, log_g:log_g, t_eff:t_eff, grid_values: vec_lf.into(), wavelengths:wavelengths }
        }
        
        fn fill_hypercube(sample_grid:& SampleGrids)->ParameterSpaceHypercube<LazyFrame>{
            //let sample_grid = Self::nadya_sample();
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
            
            //let mut slice:[f64;2] = [0.0;2];
            
            //let wavelength:Vec<f64>=vec![4000.0,4000.1];
            //let index_wavelengths = 0usize;
            //slice.copy_from_slice(&wavelength[index_wavelengths..=index_wavelengths+1]);
            //oordinates_in_parameter_space.push(slice);
            
            let mu_val = 0.66;
            let index = find_mu_index(mu_val,&sample_grid.mu_values);
            let names=vec![
                format!("wavelength"),
                format!("mu{}_flux",index),
                format!("mu{}_flux",index+1),
                format!("mu{}_cont",index),
                format!("mu{}_cont",index+1)];
            let cols:Vec<Expr> = names.iter().map(|x| col(x)).collect();
            //slice.copy_from_slice(&sample_grid.mu_values[index..=index+1]);
            //coordinates_in_parameter_space.push(slice);
            let mut corner_values:Vec<LazyFrame> = Vec::with_capacity(4usize);
    
            for i in 0..2usize{//teff
                for j in 0..2usize{//logg
                    let grid_number = 2*i+j;
                    corner_values.push(sample_grid.grid_values[grid_number].clone().
                    select(cols.clone())
                    );
                }
            }

            ParameterSpaceHypercube { fractional_coordinates: coordinates_in_parameter_space, fractional_distances: vec![0.0;2], corner_values:corner_values, partial_interpolations: vec![DataFrame::empty().lazy();2usize.pow(3)] }
        }
    }
    
    fn manual_grid_interpolation()->DataFrame{
        let coordinates:Vec<f64> = vec![22000.0,4.32];//,4000.03,0.66];
        let nadya_vec_df = SampleGrids::nadya_vec_df();
        let mut vec_lf:Vec<LazyFrame> = Vec::new();
        for item in nadya_vec_df.iter(){
            let item_lf = item.clone().lazy();
            vec_lf.push(item_lf.clone());
        }
        
        let sample_grid = SampleGrids::nadya_sample(&vec_lf);
        let mut hypercube = SampleGrids::fill_hypercube(&sample_grid);

        hypercube.get_fractional_distances(&coordinates).unwrap();

        //first 4 partial lineal interpolations; this are done on temperature.
        let mut d = hypercube.fractional_distances[0];
        let mut c1: Vec<LazyFrame> = Vec::with_capacity(2usize);
        for chunk in hypercube.corner_values.chunks(2usize){
            c1.push(LazyFrame::linear_interpolation(&chunk[0],&chunk[1], d));
        }
        if c1.len()!= 2 {panic!("not appropriate size, step 1")};
        d = hypercube.fractional_distances[1];        
        let c2 =LazyFrame::linear_interpolation(&c1[0],&c1[1], d);
        
        c2.collect().unwrap()

    }

    #[test]
    fn automatic_and_manual_equal(){
        let coordinates:Vec<f64> = vec![22000.0,4.32];//,4000.03,0.66];
        let nadya_vec_df = SampleGrids::nadya_vec_df();
        println!("nadyaframes created");
        let mut vec_lf:Vec<LazyFrame> = Vec::new();
        for item in nadya_vec_df.iter(){
            let item_lf = item.clone().lazy();
            vec_lf.push(item_lf.clone());
        }
        println!("vec_lf created");
        let sample_grid = SampleGrids::nadya_sample(&vec_lf);
        let mut hypercube = SampleGrids::fill_hypercube(&sample_grid);
        println!("hypercube created");
        let df_manual = manual_grid_interpolation();
        println!("manual df = {:?}",df_manual.head(Some(2usize)));

        let df = hypercube.multilinear_interpolation(&coordinates).unwrap().collect().unwrap();

        println!("automatic df = {:?}",df.head(Some(2usize)));
        assert_eq!(df,manual_grid_interpolation())
    }

}

