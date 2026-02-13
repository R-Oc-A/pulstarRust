use pyo3::prelude::*;
/*use pulstar::*;
use profile::*;
use temp_name_lib::*;
use pyo3::exceptions::PyValueError;
*/

mod make_pulstar_input;


/// A Python module implemented in Rust.
#[pymodule]
mod pulstar_py {
    use pyo3::prelude::*;

    #[pyclass]
    #[derive(FromPyObject)]
    struct StarParameters{
        mass:f64,
    }
/*
    #[pyclass]
    struct ModeParameters{
        l:f64
    }

    #[pyclass]
    struct ProfileParameters{
        wavelength:f64
    }
*/

    #[pymethods]
    impl StarParameters{
        #[new]
        fn new(mass:f64)->PyResult<Self>{
            Ok(
                Self{mass:mass}
            )
        }
        #[setter]
        fn set_mass(&mut self, mass:f64)->PyResult<()>{
            if mass<=0.0 {
                Err(pyo3::exceptions::PyValueError::new_err("Please enter an adequate mass for your star."))
            } else {
                self.mass = mass;
                Ok(())
            }
        }
        #[getter]
        fn get_mass(&self)->PyResult<f64>{
            Ok(self.mass)
        }
    }

    /// Formats the sum of two numbers as string.
    #[pyfunction]
    fn get_intensity_flux(
        star_parameters: StarParameters, 
/*        mode_parameters: ModeParameters,
        profile_parameters: ProfileParameters,
*/    ) -> PyResult<String> {
        Ok(format!("{} is the star mass",star_parameters.mass).to_string())
    }

    #[pyfunction]
    fn propulse(path_to_profile_input:String,path_to_pulstar_input:String)->PyResult<()>{
        match pulstar::pulstar_mkr::pulstar_main(&path_to_pulstar_input){
            Some(star_df)=>{
                profile::profile_mkr::profile_main(&path_to_profile_input,star_df);
            }
            None => {println!("unable to create rasterized star")}
        };
        Ok(())        
    }

}

// create pulstar_input.toml (both file and string)
/* */


// create profile_input.toml (both file and string)
// run pulstar from pulstar_py-> done
// run profile from pulstar_py->Done