use crate::type_def::VectorComponents;
use crate::utils::MathErrors;


///Given the components of both the surface normal vector and the unit vector
///pointing towards the observer, this function gives the cosine of the angle
///between those two vectors. 
/// 
/// The z-axis coincides with the rotation axis.
pub fn cos_chi(surface_normal: &VectorComponents, 
               k: &VectorComponents) -> Result<f64,MathErrors> {

    match surface_normal.base == k.base{

    true=> {
        let length = surface_normal.coords.norm();
        let machine_precision = 1.0e-15;
        
        match length > machine_precision{
            true => {

                Ok( surface_normal.coords.dot(&k.coords)
                /(surface_normal.coords.norm()*k.coords.norm()) )
            }
            false => {
                Err(MathErrors::VectorLengthZero)
            }

        }
    }
    false => {
        Err(MathErrors::DifferentVectorBase)
    }
 
    }
}

///Given the components of both the surface normal vector and the unit vector
///pointing towards the observer, this function gives the projected area
///of a surface cell
/// 
/// The z-axis coincides with the rotation axis.
pub fn project_area(surface_normal: &VectorComponents,
                    k: &VectorComponents)->Result<f64,MathErrors>{
    match surface_normal.base == k.base {
        true => { Ok(surface_normal.coords.dot(&k.coords))}
        false =>{Err(MathErrors::DifferentVectorBase)}
    }
}