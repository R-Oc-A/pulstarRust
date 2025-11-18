use super::*;


impl SpectralGrid {
    /// This function  computes the specific intensities on a grid using it's limb darkening coefficients
    /// ### Arguments: 
    /// `cell` - A &[SurfaceCell] that has the projection angle with the line of sight `coschi`.
    /// ### Returns: 
    /// This function doesn't return anything but it updates the flux and continuum vectors contained in the [GridsData] instance. 
    pub fn compute_flux_en_continuum(&mut self,cell: &SurfaceCell){
        let mu = cell.coschi.sqrt();

        let bcoef = 1.0 - mu; 
        let ccoef = 1.0 - cell.coschi;
        let dcoef = 1.0 - mu.powi(3);



        /*for i in self.grids_indices.iter(){
            //for (j,amplitude) in self.flux[*i].iter_mut().enumerate(){
            for  j in 0..self.grid_wavelengths.len(){
                self.flux[*i][j] = self.limb_coef_flux[*i][0][j]
                    +bcoef*self.limb_coef_flux[*i][1][j]
                    +ccoef*self.limb_coef_flux[*i][2][j]
                    +dcoef*self.limb_coef_flux[*i][3][j];

                self.continuum[*i][j] = self.limb_coef_cont[*i][0][j]
                    +bcoef*self.limb_coef_cont[*i][1][j]
                    +ccoef*self.limb_coef_cont[*i][2][j]
                    +dcoef*self.limb_coef_cont[*i][3][j];
            }
        }*/
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


