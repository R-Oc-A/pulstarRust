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

/// Triangularization of a flattened
pub fn roche_model(w:f64,delta_t:f64)->Result<Tetrahedrization,MathErrors>{
    if delta_t < 0.2 || delta_t>0.7{ Err(MathErrors::ResolutionNotSupported)}
    else{
    let ww=w;
    let potential = move |point_coords:&Vector3<f64>|{
        let rotation_freq = ww;
        1.0/point_coords.norm() 
        + 0.5 * rotation_freq *(point_coords.x.powi(2) + point_coords.y.powi(2))
        - 1.0
    };
    let grad_potential = move |point_coords:&Vector3<f64>|{
        let rotation_freq = w;
        -point_coords/(point_coords.norm().powi(3))
        +rotation_freq * (Vector3::<f64>::x()*point_coords.x
        +Vector3::<f64>::y()*point_coords.y)
    };
    let starting_point:Vector3<f64> = Vector3::from([0.0,0.0,1.1]);

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