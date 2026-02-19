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
    struct StarData{
        mass:f64,
        radius:f64,
        effective_temperature:f64,
        vsini:f64,
        inclination_angle:f64,
    }

    #[pyclass]
    #[derive(FromPyObject)]
    struct PulsationMode{
        l:u16,
        m:i16,
        rel_dr:f64,
        k:f64,
        frequency:f64,
        phase_offset:f64,
        rel_dtemp:f64,
        phase_rel_dtemp:f64,
        rel_dg:f64,
        phase_rel_dg:f64,
    }

    #[pyclass]
    #[derive(FromPyObject)]
    struct ModeCollection{
        modes:Vec<PulsationMode>
    }

    
    #[pyclass]
    #[derive(FromPyObject)]
    struct TimeTypeUniform{
        start:f64,
        end:f64,
        step:f64,
    }

    #[pyclass]
    #[derive(FromPyObject)]
    struct TimeTypeExplicit{
        collection:Vec<f64>,
    }

    #[pyclass]
    enum TimeType {
        Uniform,
        Explicit
    }


    #[pyclass]
    struct ProfileParameters{
        wavelength:f64
    }

///Estoy aqui
    #[pymethods]
    impl StarData{
        #[new]
        #[pyo3(signature=(mass = 10.0, radius = 6.93, effective_temperature=15000.0, vsini=20.0,inclination_angle=45.0))]
        fn new(mass:f64,radius:f64,effective_temperature:f64,vsini:f64,inclination_angle:f64)->PyResult<Self>{
            Ok(
                Self{mass:mass,
                    radius:radius,
                    effective_temperature:effective_temperature,
                    vsini:vsini,
                    inclination_angle:inclination_angle,
                }
            )
        }
    }

    #[pymethods]
    impl PulsationMode{
        #[new]
        #[pyo3(signature=(l=4,m=1,rel_dr=0.024,k=0.05,frequency=6.74,phase_offset=0.0,rel_dtemp=2.62,phase_rel_dtemp=180.0,rel_dg=10.0,phase_rel_dg=34.0))]
        fn new(
            l:u16,m:i16,
            rel_dr:f64,k:f64,
            frequency:f64,phase_offset:f64,
            rel_dtemp:f64,phase_rel_dtemp:f64,
            rel_dg:f64,phase_rel_dg:f64)->PyResult<Self>{
                Ok(Self{
                    l:l,
                    m:m,
                    rel_dr:rel_dr,
                    k:k,
                    frequency:frequency,
                    phase_offset:phase_offset,
                    rel_dtemp:rel_dtemp,
                    phase_rel_dtemp:phase_rel_dtemp,
                    rel_dg:rel_dg,
                    phase_rel_dg:phase_rel_dg,
                })
            }
    }

    #[pymethods]
    impl TimeTypeUniform{
        #[new]
        #[pyo3(signature=(start=0.0,end=0.1,step=0.01))]
        fn new(start:f64,end:f64,step:f64)->PyResult<Self>{
            Ok(Self{start:start,
            end:end,
            step:step})
        }
    }

    #[pymethods]
    impl TimeTypeExplicit{
        #[new]
        #[pyo3(signature=(collection=vec![0.01,0.3,0.32]))]
        fn new(collection:Vec<f64>)->PyResult<Self>{
            Ok(Self{collection:collection})
        }
    }


    #[pyfunction]
    fn propulse(profile_input:String,pulstar_input:String)->PyResult<()>{
        match pulstar::pulstar_mkr::pulstar_main(&pulstar_input){
            Some(star_df)=>{
                profile::profile_mkr::profile_main(&profile_input,star_df);
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