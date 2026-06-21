use std::ops::Add;

use temp_name_lib::interpolation::ParameterSpaceHypercube;

use crate::famias_profiles::parse_famias_grid::{give_left_right_wavelengths, parse_lib_coefs};

use super::*;

mod parse_famias_grids;
mod init_gaussian_profile;
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

pub fn gaussian_profile_mkr(toml_string:&str,star_df:DataFrame){
   //---------------------------------------- 
   //------Parsing profile_input.toml--------
   //----------------------------------------
   // |--> Check that the toml file exists
   // |--> Check if the Profile_input.toml is well written.
   // |--> Check if the Intensity Grid files exist.
   // |--> Initialize the profile parameters.
    //let profile_config = ProfileConfig::read_from_toml(toml_string);
    let mut profile_gauss:GaussianProfile=init_gaussian_profile::init_profile(toml_string);
   //---------------------------------------- 
   //----Parsing rasterized_star.parquet-----
   //----------------------------------------
    let lf = star_df.lazy();
    let tf = lf.clone().select([col("time").unique(),]).collect().unwrap();
    let extract_time_series = tf.column("time").unwrap();
    let time_points:Vec<f64> = extract_time_series.f64().unwrap().into_iter().flatten().collect();


    //init output dataframe
    let mut intensity_collection = crate::utils::IntensityFlux::new(time_points.len());


    // load limb darkening coefficients
    
   let (
        mut spectral_grid,
        mut hypercube2d
    )= loading_intensity_grids(lf.clone(), & profile_config);
    //----------------------------------------------------------------
    //-------------- Collect fluxes for each time point  -------------
    //----------------------------------------------------------------
    for pulsation_phase in time_points.iter() {
    
        intensity_collection.append_fluxes(fluxes.integrate(
            lf.clone(),
            *pulsation_phase,
            & mut spectral_grid,
            & mut hypercube2d));
        println!("done computing flux");

        println!("finished collecting fluxes {}",pulsation_phase);

        //append_into_df
    }
    
    //write output into parquet file
    intensity_collection.collect_into_single_df()

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
        let coords = [t_eff,log_g,self.central_wavelength];
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

        let (teffs,loggs,
            leftfilter_a1,leftfilter_a2,leftfilter_a3,leftfilter_a4,
            rightfilter_a1,rightfilter_a2,rightfilter_a3,rightfilter_a4) = parse_lib_coefs(central_wavelength, t_eff, log_g);
        
        let (left_wl,right_wl)=give_left_right_wavelengths(central_wavelength);

        let coords1 = [teffs[0],teffs[2]];
        let coords2 = [loggs[0],loggs[2]];
        let coords3 = [left_wl,right_wl];
        
        let mut vertices_data:Vec<LimbDarkeningCoefficients>=Vec::new();

        for i in 0usize..2{
            for j in 0usize..2{
                let index = 2usize.pow(i as u32)+j;
                for k in 0usize..2{
                    let limb = match k{
                        0usize=>{
                            LimbDarkeningCoefficients([leftfilter_a1[index],leftfilter_a2[index],leftfilter_a3[index],leftfilter_a4[index]])
                        }
                        1usize=>{
                            LimbDarkeningCoefficients([rightfilter_a1[index],rightfilter_a2[index],rightfilter_a3[index],rightfilter_a4[index]])

                        }
                        _=>{
                            LimbDarkeningCoefficients([leftfilter_a1[index],leftfilter_a2[index],leftfilter_a3[index],leftfilter_a4[index]])
                        }
                    };
                    vertices_data.push(limb)
                }
            }
        }

        new_cube.fill_coordinates(&[coords1,coords2,coords3]);
        new_cube.fill_vertices_data(&vertices_data);

        new_cube
    }


}

mod parse_famias_grid{
    use super::*;
    const STROM_FILTER_CENTRAL_WAVELENGTH:[f64;4]=[3500.0, 4110.0, 4670.0, 5470.0];//in kelvin

    pub fn give_left_right_wavelengths(central_wavelength:f64)->(f64,f64){
        let filter_wl:Vec<f64> = Vec::from(STROM_FILTER_CENTRAL_WAVELENGTH.clone());
        let (left_index,_)=filter_wl.iter().enumerate()
        .fold((0usize,STROM_FILTER_CENTRAL_WAVELENGTH[0]),
        |(index_acc,acc),(index,filter_wavelength)|
        {if *filter_wavelength <= central_wavelength {(index,*filter_wavelength)}
        else{(index_acc,acc)}});
        (STROM_FILTER_CENTRAL_WAVELENGTH[left_index],STROM_FILTER_CENTRAL_WAVELENGTH[left_index+1])
    }
    fn get_column_names(central_wavelength:f64)->Vec<String>{
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
        let mut col_names:Vec<String> = Vec::with_capacity(10);
        col_names.push(format!("Teff"));
        col_names.push(format!("logg"));
        for i in 1..=4{
            col_names.push(format!("{}_a{}",namel,i))
        }
        for i in 1..=4{
            col_names.push(format!("{}_a{}",namer,i))
        }
        col_names

    }

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

    pub fn parse_lib_coefs(central_wavelength:f64,t_eff:f64,log_g:f64)->
    (Vec<f64>,//teff
    Vec<f64>,//logg
    Vec<f64>,//leftfilter_a1
    Vec<f64>,//leftfilter_a2
    Vec<f64>,//leftfilter_a3
    Vec<f64>,//leftfilter_a4
    Vec<f64>,//rightfilter_a1
    Vec<f64>,//rightfilter_a2
    Vec<f64>,//rightfilter_a3
    Vec<f64>)//rightfilter_a4
    {   let column_names = get_column_names(central_wavelength);
        let columns:Vec<Expr> = column_names.clone().into_iter().map(|x|col(x)).collect();
        let path = format!("./grids/FAMIAS_grids/limbcoef_MHp00.stromgren");
        let df_grids = open_famias_grid(&path);
        let trimmed1= df_grids.clone().lazy().select(columns).collect().unwrap();
        let trimmed2 = trim_teff_logg(t_eff, log_g, &trimmed1.clone());

        let teffs = extract_column_as_vectorf64(&column_names[0], &trimmed2);
        let loggs = extract_column_as_vectorf64(&column_names[1], &trimmed2);
        let leftfilter_a1 = extract_column_as_vectorf64(&column_names[2], &trimmed2);
        let leftfilter_a2 = extract_column_as_vectorf64(&column_names[3], &trimmed2);
        let leftfilter_a3 = extract_column_as_vectorf64(&column_names[4], &trimmed2);
        let leftfilter_a4 = extract_column_as_vectorf64(&column_names[5], &trimmed2);
        let rightfilter_a1 = extract_column_as_vectorf64(&column_names[6], &trimmed2);
        let rightfilter_a2 = extract_column_as_vectorf64(&column_names[7], &trimmed2);
        let rightfilter_a3 = extract_column_as_vectorf64(&column_names[8], &trimmed2);
        let rightfilter_a4 = extract_column_as_vectorf64(&column_names[9], &trimmed2);

        (teffs,loggs,
        leftfilter_a1,leftfilter_a2,leftfilter_a3,leftfilter_a4,
        rightfilter_a1,rightfilter_a2,rightfilter_a3,rightfilter_a4,
        )
    }
}

