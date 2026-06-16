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
                y_gauss, wavelength:sampling_wavelengths.to_vec(),
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

impl LimbDarkeningCoefficients{

    fn new_parameter_space_cube(central_wavelength:f64,t_eff:f64,log_g:f64,df:&DataFrame)->ParameterSpaceHypercube<Self>{
        let mut new_cube = ParameterSpaceHypercube::<Self>::new(2);
        
        
        let teff= extract_column_as_vectorf64("Teff", df);
        let logg= extract_column_as_vectorf64("Teff", df);
        let l_coeffs1= extract_column_as_vectorf64("Teff", df);
        let teff= extract_column_as_vectorf64("Teff", df);
        new_cube
    }


}

mod parse_famis_grid{
    use super::*;
    const STROM_FILTER_CENTRAL_WAVELENGTH:[f64;4]=[3500.0, 4110.0, 4670.0, 5470.0];//in kelvin


    pub fn open_famias_grid(path:&str)->DataFrame{
        let schema:Vec<Field> = vec![
            Field::new("Teff".into(),DataType::Float64),
            Field::new("logg".into(),DataType::Float64),
            Field::new("M/H".into(),DataType::Float64),
            Field::new("u_a1".into(),DataType::Float64),
            Field::new("u_a2".into(),DataType::Float64),
            Field::new("u_a3".into(),DataType::Float64),
            Field::new("u_a4".into(),DataType::Float64),
            Field::new("v_a1".into(),DataType::Float64),
            Field::new("v_a2".into(),DataType::Float64),
            Field::new("v_a3".into(),DataType::Float64),
            Field::new("v_a4".into(),DataType::Float64),
            Field::new("b_a1".into(),DataType::Float64),
            Field::new("b_a2".into(),DataType::Float64),
            Field::new("b_a3".into(),DataType::Float64),
            Field::new("b_a4".into(),DataType::Float64),
            Field::new("y_a1".into(),DataType::Float64),
            Field::new("y_a2".into(),DataType::Float64),
            Field::new("y_a3".into(),DataType::Float64),
            Field::new("y_a4".into(),DataType::Float64),
        ];
        let path = format!("{}",path);
        let df = LazyCsvReader::new(path)
        .with_has_header(true)
        .with_separator(b' ')
        .with_schema(Some(Arc::new(
            Schema::from_iter(schema))))
        .finish().unwrap()
        .collect().unwrap();
        df
    }

    fn trim_df(central_wavelength:f64,lf:LazyFrame)->DataFrame{
        let filter_wl:Vec<f64> = Vec::from(STROM_FILTER_CENTRAL_WAVELENGTH.clone());
        if central_wavelength<STROM_FILTER_CENTRAL_WAVELENGTH[0]{panic!("central wavelenght is out of bounds, it should be between {} and {} Angstroms",STROM_FILTER_CENTRAL_WAVELENGTH[0],STROM_FILTER_CENTRAL_WAVELENGTH[3])};
        if central_wavelength>STROM_FILTER_CENTRAL_WAVELENGTH[3] {panic!("central wavelenght is out of bounds, it should be between {} and {} Angstroms",STROM_FILTER_CENTRAL_WAVELENGTH[0],STROM_FILTER_CENTRAL_WAVELENGTH[3])};
        let (left_index,_)=filter_wl.iter().enumerate()
        .fold((0usize,STROM_FILTER_CENTRAL_WAVELENGTH[0]),
        |(index_acc,acc),(index,filter_wavelength)|
        {if *filter_wavelength <= central_wavelength {(index,*filter_wavelength)}
        else{(index_acc,acc)}});
        let (namel,namer):(char,char) = match left_index{
            0usize=>{('u','b')}
            1usize=>{('b','v')}
            2usize=>{('v','y')}
            3usize=>{('y','y')}
            _=>{(' ',' ')}
        };
        let mut col_names:Vec<Expr> = Vec::with_capacity(10);
        col_names.push(col("Teff"));
        col_names.push(col("logg"));
        for i in 1..=4{
            col_names.push(col(format!("{}_a{}",namel,i)))
        }
        for i in 1..=4{
            col_names.push(col(format!("{}_a{}",namer,i)))
        }
        lf.clone().select(
            col_names
        ).collect().unwrap()
    }
    fn trim_teff_logg(t_eff:f64,log_g:f64,df:&DataFrame)->DataFrame{
        let teffs = extract_column_as_vectorf64("Teff", df);
        let min_teff = teffs.iter().fold(teffs[0],|acc,t|{if *t<t_eff {*t}else{acc}});
        let max_teff =  teffs.iter().fold(teffs[0],|acc,t|{if acc>=t_eff {acc}else{*t}});
        let ddf = df.clone().lazy().filter(col("Teff").eq(lit(min_teff)).or(col("Teff").eq(lit(max_teff))))
        .sort(["Teff","logg"],Default::default())
        .collect().unwrap();
        
        let loggs = extract_column_as_vectorf64("logg", &ddf);
        let min_logg = loggs.iter().fold(loggs[0],|acc,lg|{ if *lg<log_g{*lg}else{acc}});
        let max_logg = loggs.iter().fold(loggs[0],|acc,lg|{ if acc<log_g{acc}else{*lg}});
        ddf.clone().lazy().filter(
            col("logg").eq(lit(min_logg)).and(col("logg").eq(lit(max_logg)))
        ).sort(["Teff","logg"],Default::default()).collect().unwrap()
    }

}