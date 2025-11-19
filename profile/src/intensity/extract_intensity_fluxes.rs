use super::*;
use temp_name_lib::interpolation::ParameterSpaceHypercube;
impl SpectralGrid{
    pub fn new_hypercube(& self)->ParameterSpaceHypercube{
        let dimension = 4usize;//T_eff,Log_g,mu,lambda
        let mut cube = ParameterSpaceHypercube::new(dimension);
        let (temps,log_g) = (self.t_eff.clone(),self.log_g.clone());
        let wavelength:[f64;2]=[0.0,0.0];
        let mu_vals:[f64;2]=[0.0,0.0];
        cube.fill_coordinates(& vec![temps,log_g,wavelength,mu_vals]).unwrap();
        cube
    }

    fn return_mu_index(&self, mu:f64)->usize{
        let mut index:usize =0;
        
        for (n,mu_val) in self.mu_values.iter().enumerate(){
            if mu<=*mu_val { index = n;
            break;}
        }
        //println!("mu = {}; between {} and {}",mu,index,index+1);
        index-1
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
        let mu_index= spectral_grid.return_mu_index(cell.coschi.sqrt());

        //fill coordinates of the hypercube in the parameter space
        //mu value
        hypercube.fractional_coordinates[3][0..=1]
            .copy_from_slice(&spectral_grid.mu_values[mu_index..=mu_index+1]);

        for (n,wavelength) in self.shifted_wavelength.iter().enumerate(){
            //------------------------------------------------
            //------Get coordinates in parameter space--------
            //------------------------------------------------

            let coordinate_point = vec![cell.t_eff,cell.log_g, *wavelength,cell.coschi.sqrt()];


            //----------------------------------------
            //------------------fill hypercube--------
            //----------------------------------------
            //fill wavelength coordinate
            let wavelength_index = spectral_grid.row_indices[2*n];
            hypercube.fractional_coordinates[2][0..=1]
                .copy_from_slice(&spectral_grid.wavelengths[wavelength_index..=wavelength_index+1]);
            
            // Fill vertices values specific intensities
            for i in 0..2usize{// effective temperature
                for j in 0..2usize{// log gravity
                    let grid_number = 2*i+j;
                    for k in 0..2usize{//wavelength
                        for l in 0..2usize{//mu value
                            let corner_value_index = l+2*k+4*j+8*i;
                            hypercube.corner_values[corner_value_index]=spectral_grid.grid_values[[grid_number,wavelength_index+k,mu_index+l]].clone();
                        }
                    }
                }
            }
            self.flux[n] += hypercube.multilinear_interpolation(&coordinate_point).unwrap() * cell.area;

            //fill vertices values continuum
            for i in 0..2usize{// effective temperature
                for j in 0..2usize{// log gravity
                    let grid_number = 2*i+j;
                    for k in 0..2usize{//wavelength
                        for l in 0..2usize{//mu value
                            let corner_value_index = l+2*k+4*j+8*i;
                            hypercube.corner_values[corner_value_index]=spectral_grid.grid_values[[grid_number,wavelength_index+k,mu_index+l+7]].clone();
                        }
                    }
                }
            }

            self.continuum[n] += hypercube.multilinear_interpolation(&coordinate_point).unwrap() * cell.area;
        }
    }


}

