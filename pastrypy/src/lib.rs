use pyo3::prelude::*;
use pulstar::pulstar_mkr;
use profile::profile_mkr;
use marching_step_triangulation;
/// A Python module implemented in Rust.
#[pymodule]
mod pastrypy {
    use polars::frame::DataFrame;
    use pyo3::prelude::*;
    use pyo3_polars::PyDataFrame;
    use crate::pulstar_mkr;
    use crate::profile_mkr;
    use crate::marching_step_triangulation;
    

    #[pyfunction]
    fn propulse(profile_input:&str,pulstar_input:&str)->PyResult<PyDataFrame>{


        println!("---------------------------");
        println!("---------------------------");
        let profile_input_rs=profile_input.replace("\\n", "\n");
        let pulstar_input_rs=pulstar_input.replace("\n",&format!("\n"));
        //println!("{}",pulstar_input_rs);

        println!("---------------------------");
        println!("---------------------------");

//        println!("{}",profile_input_rs);

        let df =match pulstar_mkr::pulstar_main(&pulstar_input_rs){
            Some(star_df)=>{
                profile_mkr::profile_main(&profile_input_rs,star_df)
            }
            None => {panic!("Unable to create rasterized star")}};
        

        let pydf = PyDataFrame(df);

        Ok(pydf)

    }        

    #[pyfunction]
    fn pulstar(pulstar_input:&str)->PyResult<PyDataFrame>{
       println!("--------------------");
       println!("--------------------");
       let pulstar_input_rs=pulstar_input.replace("\n",&format!("\n"));
       let df= match pulstar_mkr::pulstar_main(&pulstar_input_rs){
        Some(star_df)=>{star_df},
        None => {panic!("unable to produce a data frame. Please check if your input was well written.")}
       };

       Ok(PyDataFrame(df))
    }

    #[pyfunction]
    fn profile(profile_input:&str,
        star_df: PyDataFrame)-> PyResult<PyDataFrame>{
            
            println!("----------------------------------------");
            println!("----------------------------------------");
            let profile_input_rs = profile_input.replace("\n",&format!("\n"));
            println!("{:?}",profile_input_rs);
            let star_df_rust: DataFrame = star_df.into();

            let df = profile_mkr::profile_main(&profile_input_rs,star_df_rust);


            Ok(PyDataFrame(df))
    }

    #[pyfunction]
    fn profile_gauss(gaussian_profile_input:&str,
        star_df:PyDataFrame)-> PyResult<PyDataFrame>{
            println!("----------------------------------------");
            println!("----------------------------------------");

            let profile_input_rs = gaussian_profile_input.replace("\n", &format!("\n"));
            let star_df_rust:DataFrame = star_df.into();
            println!("{:?}",star_df_rust.head(Some(5)));
            let df = profile_mkr::profile_gauss_main(&profile_input_rs,star_df_rust);

            Ok(PyDataFrame(df))
        }


    #[pyfunction]
    fn sphere_triangulation_points()
    -> PyResult<PyDataFrame>{
        let mut tetra = marching_step_triangulation::Tetrahedrization::default();
        let points=PyDataFrame(tetra.extract_points());
        Ok(points)
    }

    #[pyfunction]
    fn sphere_triangulation_triangles()
    -> PyResult<PyDataFrame>{
        let mut tetra = marching_step_triangulation::Tetrahedrization::default();
        let points=PyDataFrame(tetra.extract_triangles());
        Ok(points)
    }
    
    #[pyfunction]
    fn roche_deformed_triangulation_points(w:f64,length:f64)
    -> PyResult<PyDataFrame>{
        let delta_t =length;
        let mut tetra = marching_step_triangulation::working_examples::roche_model(w, delta_t)
        .expect("unable to create triangularization");
        let points=PyDataFrame(tetra.extract_points());
        Ok(points)
    }

    #[pyfunction]
    fn roche_deformed_triangulation_triangles(w:f64,length:f64)
    -> PyResult<PyDataFrame>{
        let delta_t =length;
        let mut tetra = marching_step_triangulation::working_examples::roche_model(w, delta_t)
        .expect("unable to create a triangularization");
        let points=PyDataFrame(tetra.extract_triangles());
        Ok(points)
    }

    

}



