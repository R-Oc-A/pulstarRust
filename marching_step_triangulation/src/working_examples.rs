use crate::na::Vector3;
use crate::{Tetrahedrization, step0};
use temp_name_lib::utils::MathErrors;


impl Default for Tetrahedrization{
    /// The [default] for a tetrahedrization is a unit sphere
    fn default()->Self{
        let delta_t = 0.3;
        //sphere potential test Unit Sphere
        let potential = |point_coords:&Vector3<f64>|{
        point_coords.norm().powi(2)-1.0}; // radius 1
        
        let grad_potential = |point_coords:&Vector3<f64>|{
            2.0*point_coords
        };

        let starting_point:Vector3<f64> = Vector3::from([-1.0,0.0,-4.0e-2]);

        tetrahedrize(delta_t, starting_point, potential, grad_potential)
    }
}

/// Triangularization of a sphere of radius 1
pub fn sphere_radius1(delta_t:f64)->Result<Tetrahedrization,MathErrors>{
    if (delta_t < 0.1 || delta_t>0.5){ 
        println!("triangle side length is {} which is outside of [0.1,0.5]",delta_t);
        Err(MathErrors::ResolutionNotSupported)}
    else{
    let potential = |point_coords:&Vector3<f64>|{
    point_coords.norm().powi(2)-1.0}; // radius 1
    
    let grad_potential = |point_coords:&Vector3<f64>|{
        2.0*point_coords
    };

    let starting_point:Vector3<f64> = Vector3::from([-1.0,0.0,-4.0e-2]);

    Ok(tetrahedrize(delta_t, starting_point, potential, grad_potential))
    }    
}

/// Triangularization of a sphere of radius 4
pub fn sphere_radius4()-> Tetrahedrization{
    let delta_t =0.3;
    let potential = |point_coords:&Vector3<f64>|{
        point_coords.norm().powi(2)-4.0
    };
    let grad_potential = |point_coords:&Vector3<f64>|{
        2.0*point_coords
    };
    let starting_point:Vector3<f64> = Vector3::from([4.2,0.0,0.0]);
    tetrahedrize(delta_t, starting_point, potential, grad_potential)
}

/// Triangularization of the Roche aproximation of a flatenned sphere
/// This model is applied to rotating stars in the classical approximation where it is assumed the point mass approximation.
/// While not necessary, in this case we assume solid body like rotation 
/// The changing of factors is taken from the page 24 of chapter 2 of the book "Mechanical equilibrium of rotating stars" by Maeder et al. 2009
/// ### Arguments:
/// * `rotation_frequency` - a f64 value that can go up to 98% critical rotation
/// * `delta_t` - The approximate length of each triangle that covers the surface in normalized units
/// ### Returns: 
/// * This function returns a [Result] where the [Ok] variant contains a [Tetrahedrization] of a flattened star and the [Err] variant is presented when the requested resolution is out of bounds. 
pub fn roche_model(rotation_frequency:f64,delta_t:f64)->Result<Tetrahedrization,MathErrors>{
    if (delta_t < 0.1 || delta_t>0.5)|| rotation_frequency>0.98{ 
        println!("triangle side length is {}",delta_t);
        println!("rotation frequency is {} critical, which is above what's supported",rotation_frequency);
        Err(MathErrors::ResolutionNotSupported)}

    else{
    let w=(8.0f64/27.0f64).sqrt()*rotation_frequency;//
    let potential = move |point_coords:&Vector3<f64>|{
        let ww = w;
        1.0/point_coords.norm() 
        + 0.5 * ww.powi(2) *(point_coords.x.powi(2) + point_coords.y.powi(2))
        - 1.0
    };
    let grad_potential = move |point_coords:&Vector3<f64>|{
        let ww = w;
        -point_coords/(point_coords.norm().powi(3))
        +ww.powi(2) * (Vector3::<f64>::x()*point_coords.x
        +Vector3::<f64>::y()*point_coords.y)
    };
    let starting_point:Vector3<f64> = Vector3::from([0.01,0.01,1.1]);

    Ok(tetrahedrize(delta_t, starting_point, potential, grad_potential))
    }
}

pub fn roche_lobe()->Tetrahedrization{
    let delta_t = 0.3;
    let q:f64 =0.5;//1.5;//1.0;//0.5;
    let f:f64 = 4.0e-3;//3.75e-3;
    let delta:f64 = 1.04e1;//1.3e1;//8.0;
    let starting_point = Vector3::from([1.1,0.0,4.1]);
    let potential = move |point_coords:&Vector3<f64>|{
        let x = point_coords.x;
        let y = point_coords.y;
        let z = point_coords.z;
        let r = point_coords.norm();
        let mod_r = ((x-delta).powi(2) + y.powi(2) + z.powi(2)).sqrt();   

        -2.5e-1 + 1.0/r + q * (1.0/mod_r - x / delta.powi(2) )
        + 0.5 * (1.0 + q) * f.powi(2) *(x.powi(2) + y.powi(2))
    };

    let grad_potential = move |point_coords:&Vector3<f64>|{
        let x = point_coords.x;
        let y = point_coords.y;
        let z = point_coords.z;
        let r = point_coords.norm();
        let mod_r = ((x-delta).powi(2) + y.powi(2) + z.powi(2)).sqrt();

        let xx = -x/r.powi(3) - q * (1.0/mod_r.powi(3)*(x-delta)
        + 1.0/delta.powi(2) ) + (1.0 + q) * f.powi(2) * x;

        let yy = -y/r.powi(3) - q * y/mod_r.powi(3)
        + (1.0 + q) * f.powi(2) * y;

        let zz = -z/r.powi(3) - q * z/mod_r.powi(3);

        
        Vector3::<f64>::from([xx,yy,zz])
    };

    tetrahedrize(delta_t, starting_point, potential, grad_potential)

}



pub fn tetrahedrize<F,Gf>(delta_t:f64,starting_point:Vector3<f64>,potential:F,grad_potential:Gf)->Tetrahedrization
where 
    F:Fn(&Vector3<f64>)->f64 + 'static,
    Gf:Fn(&Vector3<f64>)->Vector3<f64> + 'static,
{    
    let mut tetra= Tetrahedrization::new(delta_t,potential,grad_potential);
    step0(&mut tetra, starting_point);
    tetra.step4();

    tetra
}