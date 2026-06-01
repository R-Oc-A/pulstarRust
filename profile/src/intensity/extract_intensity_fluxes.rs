use super::*;
use temp_name_lib::{interpolation::ParameterSpaceHypercube, utils::MathErrors};

impl SpectralGrid{

    pub fn new_hypercube(& self,dimension:usize)->Result<ParameterSpaceHypercube<LazyFrame>,MathErrors>{
        let mut cube = ParameterSpaceHypercube::<LazyFrame>::new(dimension);
        let (temps,log_g) = (self.t_eff.clone(),self.log_g.clone());
        match dimension{
            2usize => {cube.fill_coordinates(&vec![temps,log_g])?;}
            _ => {panic!("there's no intensity grids for variable metalicity and chemical abundances.")}
        }
        //initializing memory for the corner values of the parameter space. 
        let dummy_df = DataFrame::empty();
        let dummy_lf = dummy_df.lazy();
        let dummy_coordinate_values:Vec<LazyFrame> = vec![dummy_lf.clone();2usize.pow(dimension as u32)];
        let dummy_partial_interpolation:Vec<LazyFrame> = vec![dummy_lf.clone();2usize.pow(dimension as u32 +1)-1];
        cube.corner_values = dummy_coordinate_values;
        cube.partial_interpolations = dummy_partial_interpolation;
        Ok(cube)
    }

    fn find_mu_index(&self, mu:f64)->(usize,f64){
        //self.mu_values must be ordered from lower to greater for this function to work.
        let (index,_) = self.mu_values.iter().enumerate().fold((0usize,self.mu_values[0]),
        |acc,(index,x)|{if mu>*x{(index,*x)}else{acc}});
        
        let fractional_distance_mu = if index == 6usize || index == 0usize{//because there's only 8 mu values
            0.0
        }else{
            (mu - self.mu_values[index])/(self.mu_values[index+1]-self.mu_values[index])
        };
        (index,fractional_distance_mu)
    }

    /*
    /// This function is used to store for the observed wavelength the indices of the wavelengths in [GridsData] that will be used for interpolation. 
    /// This function relies on the bisection algorithm to perform the query.
    pub fn extract_important_rows(&mut self,global_flux: &mut FluxOfSpectra){
        self.row_indices.fill(0);

        let mut counter = 0usize;
        /*for shifted_wavelength in global_flux.shifted_wavelength.iter() {
            self.row_indices[counter] = search_geq(&self.wavelengths, *shifted_wavelength)-1;
            self.row_indices[counter+1] = self.row_indices[counter]+1;
            counter += 2usize;
        }*/
    }*/

}
impl FluxOfSpectra{


/// This function that returns the intensity flux and the continuum flux interpolated from the intensity grids
/// for each surface cell of the rasterized sphere.
/// ### Arguments:
/// * `cell` -  A &[SurfaceCell] that contains the `coschi`-Projection of the normal vector of the surface cell with the unit vector towards the observer;
///  `temperature` - temperature value over the surface cell of the rasterized star (in kelvin);
///  `log_gravity` - log g value over the surface cell of the rasterized star;
///  `relative doppler shift` - relative doppler wavelength shift, this is of course related to the velocity.
///  'area' - area of the surface cell of the rasterized star;
/// * `spectral_grid` - a reference to an instance of [SpectralGrid] that contains all of the relevant information parsed from the intensity grid files, or the neural network regressor. 
/// * `hypercube' - a reference to an instance of [ParameterSpaceHypercube]. The coordinates and values of the vertices in the parameter space that are used to perform the multilinear interpolation. 
/// ### Returns:
/// - This function adds the contribution of the observed specific intensities by a surface cell. 
pub fn collect_flux_from_cell(& mut self,
    cell: & SurfaceCell,
    spectral_grid: &mut SpectralGrid,
    hypercube: &mut ParameterSpaceHypercube<LazyFrame>,
    collecting_lf:LazyFrame)->LazyFrame{

        let mu_val = cell.coschi.sqrt();
        spectral_grid.fill_corner_values_2d(mu_val, hypercube);
        let coords_in_param_space = vec![cell.t_eff,cell.log_g];
        let lf = hypercube.multilinear_interpolation(&coords_in_param_space).unwrap();
        let shifted_wavelengths = self.get_doppler_shifted_wavelengths(cell);
        let linear_lf = parse_intensity_grids::wavelength_interpolation(shifted_wavelengths, lf.clone());

        let flux_lf = linear_lf.clone().select(
            [col("pixel_id"),
            col("wavelength"),
            (col("mu_avg_s") * lit(cell.area)).alias("flux"),
            (col("mu_avg_c") * lit(cell.area)).alias("continuum")]
        );

        //let linear_df = flux_lf.clone().collect().unwrap();
        //println!("linear_df {:#?}",linear_df.head(Some(5)));
        FluxOfSpectra::add_into_current_data(collecting_lf.clone(),flux_lf.clone())
        //println!("self df {:#?}",self.flux_data.head(Some(5)));

}

fn add_into_current_data (total_flux:LazyFrame, cell_contribution:LazyFrame)->LazyFrame{
    let flux_lf = total_flux;
    let cell_lf = cell_contribution.clone();

    let flux_lf_to_add = flux_lf.clone().join(
        cell_lf,
        [col("pixel_id")],
        [col("pixel_id")],
        JoinArgs::new(JoinType::Inner)
    );

    let flux_added = flux_lf_to_add.clone().select([
        col("wavelength"),
        col("pixel_id"),
        col("time"),
        (col("flux") + col("flux_right")).alias("flux"),
        (col("continuum") + col("continuum_right")).alias("continuum")
    ]);

    flux_added
}

}

impl SpectralGrid{
    /// This function fills in the corner values of the [ParameterSpaceHypercube] with [LazyFrame]s of the [SpectralGrid]s 
    ///     
    /// ### Arguments:
    /// * `mu_val` - A [f64] variable that contains the Cos(θ) related to the point of view towards the observer
    /// * `hypercube` - An instance of a parameter space hypercube.
    fn fill_corner_values_2d(&mut self, mu_val:f64,hypercube:&mut ParameterSpaceHypercube<LazyFrame>){
        let (mu_index,fractional_distance_mu) = self.find_mu_index(mu_val);
        for i in 0..2usize{// Effective temperature
            for j in 0..2usize{// log gravity
                let corner_value_index = 2*i + j;
                hypercube.corner_values[corner_value_index] = 
                    avg_mu_lazyframe(self.grid_values[corner_value_index].clone().lazy().clone(),mu_index,fractional_distance_mu);
                }
            }
    } 
    
}

    /// This function takes a [LazyFrame] of the wavelength spectrum intensities for a given temperature, log_g, and in the future another coordinates (such as metalicity and chemical abundances) and returns an
    /// another [LazyFrame] with just 3 columns:
    /// 
    /// `|wavelength|mu_avg_s|mu_avg_c|`
    /// 
    /// With the purpose of interpolating only in wavelength at the very end. 
    fn avg_mu_lazyframe(lf:LazyFrame,mu_index:usize,fractional_distance:f64)->LazyFrame{
        let expr_final = if mu_index == 6usize || mu_index==0usize{
            vec![
                col("wavelength"),
                col(format!("mu{}_s",mu_index+1)).alias("mu_avg_s"),
                col(format!("mu{}_c",mu_index+1)).alias("mu_avg_c"),
            ]
        }else{
            let names = vec![
                format!("wavelength"),
                format!("mu{}_s",mu_index+1),
                format!("mu{}_s",mu_index+2),
                format!("mu{}_c",mu_index+1),
                format!("mu{}_c",mu_index+2),
            ];
            
            let cols:Vec<Expr> = names.iter().map(|x| col(x)).collect();

            vec![cols[0].clone(),//wavelength
                (cols[1].clone() * lit(fractional_distance) + cols[2].clone() * lit(1.0 - fractional_distance)).alias("mu_avg_s"),//mu_average_s
                (cols[3].clone() * lit(fractional_distance) + cols[4].clone() * lit(1.0 - fractional_distance)).alias("mu_avg_c")].clone()//mu_average_c
        };
        lf.clone().select(expr_final)
    }
