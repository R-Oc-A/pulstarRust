use std::ops::Add;

use temp_name_lib::interpolation::ParameterSpaceHypercube;

use super::*;

mod parse_famias_grids;
///This structure contains all of the variables defined in FAMIAS to produce a Gaussian  profile
pub struct GaussianProfile{
    ///Intensity of the flux (amplitude of the gaussian depending on the surface cell)
    fl_in_ul:Vec<f64>,
    //There seems to be a flIn in famias which is defined as the norm of a vector times cos(χ) so it must be the proyection
    ///Total sum of fl_in_ul
    continuum:f64,
    ///Gaussian profile
    y_gauss:Vec<f64>,
    //Wavelengths,
    wavelength:Vec<f64>,
    ///Equivalent width
    eq_w:f64,
    ///Alpha W
    alpha_w:f64,
    ///This coefficient I don't know quite sure what it means
    sigmag_sqrtpi_sqrt2:f64,
    ///Neither this one. 
    sigmag_sqrt2_pow2:f64,
    ///Zero point shift of velocity
    zero_point_shift:f64,
    ///Central Wavelength
    central_wavelength:f64,
    ///Limb darkening coefficients
    limb:LimbDarkeningCoefficients,
    //Star temperature
    t_eff:f64,
    //Star Logg
    log_g:f64,
    ///Output Dataframe
    output:DataFrame,
}

///This structure contains the 4 limb darkening coefficients. 
#[derive(Clone, Copy)]
pub struct LimbDarkeningCoefficients([f64;4]);

impl GaussianProfile{
    pub fn new(number_of_surface_cells:usize,
        sampling_wavelengths:&[f64],central_wavelength:f64)->Self{
            let continuum = 1.0;
            let eq_w =1.0;
            let alpha_w =1.0;
            let sigmag_sqrtpi_sqrt2=1.0;
            let sigmag_sqrt2_pow2 = 1.0;
            let fl_in_ul:Vec<f64> = vec![0.0;number_of_surface_cells];
            let y_gauss:Vec<f64> = vec![0.0;sampling_wavelengths.len()];
            let zero_point_shift = 0.0;
            let limb = LimbDarkeningCoefficients([0.0;4]);
            let star_temperature = 10.0;
            let Star_logg =3.8;
            let output:DataFrame = DataFrame::empty();
            GaussianProfile { fl_in_ul, continuum, 
                y_gauss, wavelength:sampling_wavelengths.clone().to_vec(),
                eq_w, alpha_w,
                sigmag_sqrtpi_sqrt2,
                sigmag_sqrt2_pow2,
                zero_point_shift,
                central_wavelength,
                limb,t_eff:star_temperature,
                log_g:Star_logg,output}
        }
    
    ///This formula is taken from Joris de Ridder Thesis. 
    fn compute_gaussian_amplitude(& mut self, shifted_wavelength:&[f64],fl_in_ul:f64,d_temperature:f64){
        //rewrite gaussianprofile by adding cell's contrubution
        let w_eintr = self.get_w_eintr(d_temperature);
        self.y_gauss = 
        self.y_gauss.iter()
        .enumerate()
        .map(|(index,intensity)|
        {   intensity + 
            fl_in_ul * (
                1.0 -w_eintr * self.sigmag_sqrtpi_sqrt2 *
                (-(self.central_wavelength - shifted_wavelength[index]).powi(2)*self.sigmag_sqrt2_pow2).exp()
            )
        }
        ).collect();
    }
    fn compute_fl_in_ul(&mut self,surface_area:f64,mu:f64,cell_index:usize){
        self.fl_in_ul[cell_index]=surface_area *(
            1.0 
            - self.limb.0[0] * (1.0 - mu)
            + self.limb.0[1] * (1.0 - mu.powi(2))
            + self.limb.0[2] * (1.0 - mu.powf(3.5))
            + self.limb.0[3] * (1.0 - mu.powi(4))
        );
    }

    fn renormalize_fl_in_ul(&mut self){
        self.continuum = self.fl_in_ul.iter().fold(0.0,|acc,x|x+acc);
        if self.continuum > 0.0{
            self.fl_in_ul= self.fl_in_ul.iter().map(|x| x/self.continuum).collect();
        }
    }

    fn update_limb_darkening_coefficients(& mut self,t_eff:f64,log_g:f64,hypercube2d:& mut ParameterSpaceHypercube<LimbDarkeningCoefficients>){
        let coords = [t_eff,log_g];
        self.limb = hypercube2d.multilinear_interpolation(&coords).unwrap();
    }

    fn get_doppler_shifted_wavelengths(&self,relative_doppler_shift:f64)->Vec<f64>{
        self.wavelength.iter().map(|lambda| lambda * relative_doppler_shift).collect()
    }

    fn get_w_eintr(&self,d_temperature:f64)->f64{
        self.eq_w*(1.0 + self.alpha_w*d_temperature)
    }


    pub fn integrate(& mut self, surface_cells:&[SurfaceCell],
        hypercube2d:& mut ParameterSpaceHypercube<LimbDarkeningCoefficients>){
        self.update_limb_darkening_coefficients(self.t_eff, self.log_g, hypercube2d);
        self.y_gauss = vec![0.0;self.y_gauss.len()];

        for (index,cell) in surface_cells.iter().enumerate(){
            //self.update_limb_darkening_coefficients(cell.t_eff, cell.log_g, hypercube2d);
            self.compute_fl_in_ul(cell.area, cell.coschi.sqrt(), index);
        }
        self.renormalize_fl_in_ul();
        for (index,cell) in surface_cells.iter().enumerate(){
            let d_temperature = (cell.t_eff/self.t_eff)-1.0;
            let shifted_wavelength = self.get_doppler_shifted_wavelengths(cell.rel_dlamb);
            self.compute_gaussian_amplitude(&shifted_wavelength, self.fl_in_ul[index],d_temperature);
        }   
    }
}

impl Add for LimbDarkeningCoefficients{
    type Output = Self;
    
    fn add(self, rhs: Self) -> Self{
        let mut new_coeffs:Self=LimbDarkeningCoefficients([0.0;4]);
        for (index,coef) in new_coeffs.0.iter_mut().enumerate(){
            *coef = self.0[index] + rhs.0[index];
        }
        new_coeffs
    }
}

impl std::ops::Mul<f64> for LimbDarkeningCoefficients{
    type Output=Self;
    fn mul(self,rhs:f64)->Self{
        let mut new_coeffs:Self = LimbDarkeningCoefficients([0.0;4]);
        for (index,coef) in new_coeffs.0.iter_mut().enumerate(){
            *coef = self.0[index] * rhs
        }
        new_coeffs
    }
}


impl temp_name_lib::interpolation::LinearlyInterpolatable for LimbDarkeningCoefficients{
    fn linear_interpolation(left:& Self,right: & Self, fractional_distance:f64)->Self{
        *left * fractional_distance + *right * (1.0-fractional_distance)
    }    
}
