use super::*;
use temp_name_lib::interpolation::polars_lazyframe;
use temp_name_lib::interpolation::ParameterSpaceHypercube;

impl SpectralGrid{
    pub fn new_hypercube(& self,dimension:usize)->ParameterSpaceHypercube<LazyFrame>{
        //let dimension = 2usize;//T_eff,Log_g
        let mut cube = ParameterSpaceHypercube::<LazyFrame>::new(dimension);
        let (temps,log_g) = (self.t_eff.clone(),self.log_g.clone());
        match dimension{
            2usize => {cube.fill_coordinates(&vec![temps,log_g]);}
            _ => {panic!("there's no intensity grids for variable metalicity and chemical abundances.")}
        }
        cube
    }

    fn find_mu_index(&self, mu:f64)->(usize,f64){
        let mut index:usize =0;
        
        for (n,mu_val) in self.mu_values.iter().enumerate(){
            if mu<=*mu_val { index = n;
            break;}
        }
        let fractional_distance_mu = if index == 7usize {//because there's only 8 mu values
            0.0
        }else{
            (mu - self.mu_values[index-1])/(mu_values[index]-mu_values[index-1])
        };
        //println!("mu = {}; between {} and {}",mu,index,index+1);
        (index-1,fractional_distance_mu)
    }

    
    /// This function is used to store for the observed wavelength the indices of the wavelengths in [GridsData] that will be used for interpolation. 
    /// This function relies on the bisection algorithm to perform the query.
    pub fn extract_important_rows(&mut self,global_flux: &mut FluxOfSpectra){
        self.row_indices.fill(0);

        let mut counter = 0usize;
        for shifted_wavelength in global_flux.shifted_wavelength.iter() {
            self.row_indices[counter] = search_geq(&self.wavelengths, *shifted_wavelength)-1;
            self.row_indices[counter+1] = self.row_indices[counter]+1;
            counter += 2usize;
        }
    }

}
impl FluxOfSpectra{


pub fn collect_flux_from_cell_version2(& mut self,
    cell: & SurfaceCell,
    spectral_grid: &mut SpectralGrid,
    hypercube: &mut ParameterSpaceHypercube<LazyFrame>){

        spectral_grid.fill_corner_values_2d(wavelenght_index, mu_val, hypercube);

    }



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
    pub fn collect_flux_from_cell(& mut self, cell: & SurfaceCell, spectral_grid: &mut SpectralGrid, hypercube:& mut ParameterSpaceHypercube){
        spectral_grid.extract_important_rows(self);   

        let mu_index = if hypercube.fractional_coordinates.len()==4{
            let index= spectral_grid.return_mu_index(cell.coschi.sqrt());
            //fill coordinates of the hypercube in the parameter space
            //mu value
            hypercube.fractional_coordinates[3][0..=1]
                .copy_from_slice(&spectral_grid.mu_values[index..=index+1]);
            index
        }
        else{ 0usize};

        for (n,wavelength) in self.shifted_wavelength.iter().enumerate(){
            //------------------------------------------------
            //------Get coordinates in parameter space--------
            //------------------------------------------------

            let coordinate_point = if hypercube.fractional_coordinates.len()==4
                {vec![cell.t_eff,cell.log_g, *wavelength,cell.coschi.sqrt()]}
                else{vec![cell.t_eff,cell.log_g,*wavelength]};


            //----------------------------------------
            //------------------fill hypercube--------
            //----------------------------------------
            //fill wavelength coordinate
            let wavelength_index = spectral_grid.row_indices[2*n];
            hypercube.fractional_coordinates[2][0..=1]
                .copy_from_slice(&spectral_grid.wavelengths[wavelength_index..=wavelength_index+1]);
            
            spectral_grid.fill_corner_values(wavelength_index, mu_index, hypercube);
            self.flux[n] += hypercube.multilinear_interpolation(&coordinate_point).unwrap() * cell.area;

            //fill vertices values continuum
            spectral_grid.fill_corner_values(wavelength_index, mu_index+7, hypercube);
            self.continuum[n] += hypercube.multilinear_interpolation(&coordinate_point).unwrap() * cell.area;
        }
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
                hypercube.corner_values[corner_values_index] = 
                    avg_mu_lazyframe(self.grid_values[corner_value_index].lazy().clone(),mu_index,fractional_distance_mu);
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
        if mu_index == 7usize {
            let names = vec![
                format!("wavelength"),
                format!("mu{}_s",mu_index),
                format!("mu{}_s",mu_index+1),
                format!("mu{}_c",mu_index),
                format!("mu{}_c",mu_index+1),
            ];
            
            let cols:Vec<Expr> = names.iter().map(|x| col(x)).collect();

            let expr_final:Vec<Expr> = vec![cols[0].clone(),//wavelength
                (cols[1].clone() * lit(fractional_distance) + cols[2].clone() * lit(1.0 - fractional_distance)).alias("mu_avg_s"),//mu_average_s
                (cols[3].clone() * lit(fractional_distance) + cols[4].clone() * lit(1.0 - fractional_distance)).alias("mu_avg_c")].clone();//mu_average_s
            lf.clone().select(expr_final)
        }else{
            let expr_final = vec![
                col("wavelength"),
                col(format!("mu{}_s",mu_index)).alias("mu_avg_s"),
                col(format!("mu{}_c",mu_index)).alias("mu_avg_c"),
            ];
            lf.clone().select(expr_final)
        }
    }

    fn fill_corner_values(&mut self,wavelength_index:usize,mu_index:usize,hypercube:&mut ParameterSpaceHypercube){
        match hypercube.fractional_coordinates.len(){
            4usize=>{self.fill_corner_values_4d(wavelength_index, mu_index, hypercube);}
            3usize=>{self.fill_corner_values_3d(wavelength_index, hypercube);}
            2usize=>{self.fill_corner_values_2d(wavelength_index, mu_val,hypercube)}
            _=>{panic!("never thought of this case; in extract intensity fluxes the dimension of the hypercube in the parameter space is {}",{hypercube.fractional_coordinates.len()})}
        }
    }
}