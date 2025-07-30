use super::Coordinates;
use super::MathErrors;
use super::na;
use std::ops;

//----------------------------------------
//-----Overloading Operators--------------
//----------------------------------------
///This is operation overloading so that we can add vectors component by component on the same basis.
impl ops::Add for Coordinates{
    type Output = Result<Coordinates,MathErrors>;
    fn add(self,other: Coordinates) -> Result<Coordinates,MathErrors>{
        match self {
            Coordinates::Spherical(value) => {
                match other{
                    Coordinates::Spherical(value2)=>{
                        Ok(Coordinates::Spherical(value + value2))
                    }
                    Coordinates::Cartesian(_)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
            Coordinates::Cartesian(value) =>{
                match other{
                    Coordinates::Cartesian(value2)=>{
                        Ok(Coordinates::Spherical(value + value2))
                    }
                    Coordinates::Spherical(_)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
        }
    }
}

///This is operation overloading so that we can addasign (self += rhs) vectors component by component having the same basis.
impl ops::AddAssign<Coordinates> for Coordinates{
    fn add_assign(&mut self, rhs: Self) {
        let dummy = (*self + rhs).unwrap();
        *self = dummy;       
    }

}

/// This is operation overloading so that we can multiply by a scalar (from the left) as `scalar * vector`
impl ops::Mul<Coordinates> for f64{
    type Output = Coordinates;

    fn mul(self, rhs: Coordinates) -> Self::Output {
        match rhs{
            Coordinates::Cartesian(value)=>{ Coordinates::Cartesian(self * value)}
            Coordinates::Spherical(value)=>{ Coordinates::Spherical(self * value)}
        }
    }
}



//----------------------------------------
//-----Methods and implementations--------
//----------------------------------------
impl Coordinates {
    /// Projects a vector onto another one
    /// ### Arguments:
    /// * `other` - A borrowed reference variant of [Coordinates]
    /// ### Returns:
    /// This method returns A [Result] that has the following values binded to it:
    /// * `Ok(f64)` - where the binded value is the projection of the vector into `other`
    /// * `Err(DiferentVectorBase)` - This returns a Different Vector Base to the error to the caller function.
    pub fn project_vector(&self,other:&Self)->Result<f64,MathErrors>{
        match self{
            Coordinates::Cartesian(value)=>{
                match other {
                    Coordinates::Cartesian(value2)=>{
                        Ok(value.dot(value2))
                    }
                    Coordinates::Spherical(_)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
            Coordinates::Spherical(value)=>{
                match other {
                    Coordinates::Spherical(value2)=>{
                        Ok(value.dot(value2))
                    }
                    Coordinates::Cartesian(_)=>{
                        Err(MathErrors::DifferentVectorBase)
                    }
                }
            }
            }
        }
    /// This method produces a unit vector k in cartesian coordinates that points towards the observer, 
    /// in this frame of reference, the z-axis coincides with the rotation axis. 
    /// 
    /// ### Arguments: 
    /// * `inclination_angle` - The inclination angle in rads with respect to the observer. 
    /// ### Returns:
    /// * `Coordinate::Cartesian(unit_vector_k)` - The unit vector k in cartesian coordinates.
    pub fn unit_vector_k(inclination_angle:f64)->Self{
        Coordinates::Cartesian(na::Vector3::new(-inclination_angle.sin(),
            0.0,
            inclination_angle.cos() ))  
    }

    /// This method provides a straightforward way to compute the norm of a vector
    /// 
    /// ### Arguments:
    /// 
    /// ### Returns: 
    /// * `f64` - the L2 norm of a vector.
    pub fn vector_length(&self)->f64{
        match self {
            Coordinates::Cartesian(value)=>{value.norm()}
            Coordinates::Spherical(value)=>{value.norm()}
        }
    }

    /// This method obtains the transformation matrix useful for changing basis
    /// 
    /// ### Arguments
    /// * `theta_rad` - The colatitude angle in radians
    /// * `phi_rad` - The azimuthal angle in radians
    /// ### Returns
    /// * `t_a` -The transformation Matrix to change from spherical to cartesian coordinates and vice versa.
    fn transformation_matrix(&self,
        theta_rad:f64,
        phi_rad:f64)-> na::Matrix3<f64>{
            let sintheta = theta_rad.sin();
            let costheta = theta_rad.cos();
            let sinphi = phi_rad.sin();
            let cosphi = phi_rad.cos();


            let t_matrix = na::Matrix3::new(
                sintheta*costheta, costheta*cosphi, -sinphi,
                sintheta*sinphi, costheta*cosphi, cosphi,
                costheta, -sintheta, 0.0,
            );//spherical to cartesian matrix

            match self{
                Coordinates::Spherical(_) => {t_matrix}
                Coordinates::Cartesian(_) => {t_matrix.transpose()}
            }
        }
    /// This method transforms coordinates from spherical to cartesian and vice versa. 
    /// 
    /// ### Arguments:
    /// * `theta_rad` - The colatitude angle in radians.
    /// * `phi_rad` - The azimuthal angle in radians.
    /// 
    /// ### Returns:
    /// * ` ` - a new instance of [Coordinates] with different base.
    pub fn transform(&self,theta_rad:f64,phi_rad:f64)->Coordinates{
        let t_matrix = self.transformation_matrix(theta_rad, phi_rad);
        match self{
            Coordinates::Cartesian(value)=>{
                Coordinates::Spherical(t_matrix * value)}
            Coordinates::Spherical(value)=>{
                Coordinates::Cartesian(t_matrix * value)}
        }
    }
    
    /// This method obtains the r component of the spherical coordinates
    /// 
    /// ### Arguments: 
    /// * `self` - A [Coordinates] instance that should be Spherical variant.
    /// 
    /// ### Returns:
    /// * `Option` - This method returns a `Some(f64)` that has the binded the r component value, in case it was casted by a sperical variant, or a  `None` otherwise.
    pub fn r_component(&self)->Option<f64>{
        match self{
            Coordinates::Spherical(value)=>{
                Some(value[0])
            }
            _=>{None}
        }
    }
}