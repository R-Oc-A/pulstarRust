use marching_step_triangulation::{
    Point,
    write_output::ExtractedTriangulation};

use temp_name_lib::utils::{MACHINE_PRECISION, MathErrors};
use crate::MeshConfig::TSphere;
use crate::na::Vector3;
//use crate::reference_frames::rotation_treatment::tar::TARCollection;
use crate::{PI, PulstarConfig, TARCollection};
use crate::reference_frames::*;
use crate::{observed_pulsation_velocity,project_vrot,local_surface_temperature_logg};


use crate::SurfaceCell;

#[derive(Clone,Debug)]
pub struct Triangles {
    pub triangles: ExtractedTriangulation,
    pub effective_temp_of_points:Vec<f64>,
    pub log_g_of_points:Vec<f64>,
    pub velocity_of_points:Vec<f64>,
    //some other quantity that I'll need to add. 
}

pub fn new_triangles(
    parameters:&PulstarConfig,
    triangulation:ExtractedTriangulation)->Triangles{
    let (log_g_of_points,
        effective_temp_of_points)=unperturbed_log_g_and_effective_temperature_for_triangles(parameters, &triangulation);
    
    let velocity_of_points:Vec<f64> = vec![0.0;effective_temp_of_points.len()];
    Triangles { triangles:triangulation, effective_temp_of_points, log_g_of_points, velocity_of_points }
}

//This function computes the unperturbed surface gravity and effective temperature for all of the points that discretize the surface of the star. 
pub fn unperturbed_log_g_and_effective_temperature_for_triangles(
    parameters:& PulstarConfig,
    generated_triangles:&ExtractedTriangulation
    )->(Vec<f64>,Vec<f64>){

        
        match &parameters.mesh{
            TSphere{triangle_length:_}=>{
                let log_g0 = 4.438 + parameters.star_data.mass.log10()
                    - 2.0 * parameters.star_data.radius.log10();
                
                let t_eff = parameters.star_data.effective_temperature;
                let number_of_points = generated_triangles.points.len();
                (vec![log_g0;number_of_points],vec![t_eff;number_of_points])
            }
            // Here I'll deal with the pain in the ass deformed sphere.
            _ =>{panic!("unexpected behaviour, this function should have been called only by meshing done via triangulation.")}
        }
}

pub fn cell_from_triangle(triangles:&Triangles,index:usize)->SurfaceCell{
    let triangulation = &triangles.triangles;
    let centroid = triangulation.centroid(index);
    let centroid = point_coords_to_spherical(centroid);
    let (_surface_normal,surface_area) = triangulation.area_and_surface_normal(index);
    
    let mut new_surface_cell = SurfaceCell::new(centroid.y,centroid.z);

    new_surface_cell.area=surface_area;
    new_surface_cell.triangle_index=Some(index);

    new_surface_cell
}

fn point_coords_to_spherical(cartesian_coords:Vector3<f64>)->Vector3<f64>{
    let x = cartesian_coords.x;
    let y = cartesian_coords.y;
    let z = cartesian_coords.z;

    let r = cartesian_coords.norm();
    let cos_theta = z/r;
    let theta = (cos_theta.acos()).to_degrees();
    
    let phi = if cos_theta.abs()<MACHINE_PRECISION{0.0}else{
        if x<MACHINE_PRECISION { 
            if y>0.0{0.5*PI}
            else{1.5 * PI}
        }else{
            let tan_phi = y/x;
            if x>0.0{
                if y>0.0{tan_phi.atan()}
                else{ tan_phi.atan()+PI}
            }else{
                if y<0.0{tan_phi.atan()+PI}
                else{tan_phi.atan()}
            }
        }        
    }.to_degrees();
    
    Vector3::from([r,theta,phi])
}

impl SurfaceCell{
    pub fn barycenter_surface_normal_area(&self,triangulation:&ExtractedTriangulation)->(Coordinates,Coordinates,f64){
        if let Some(triangle_index)= self.triangle_index{

            let mut centroid = triangulation.centroid(triangle_index);//
            let (mut surface_normal,area) = triangulation.area_and_surface_normal(triangle_index);

            if centroid.dot(&surface_normal)>0.0 {surface_normal=-surface_normal};//Make sure that the surface normal is positively oriented. For more complex geometries I might need to reconsider this approach.
            centroid = point_coords_to_spherical(centroid);
            (
                Coordinates::Spherical(centroid),
                Coordinates::Cartesian(surface_normal),
                area
            )

        }else{panic!("barycenter computation on a surface cell not computed via a triangulation")}
    }
}


impl Triangles {
    /*
    pub fn update_points_of_the_triangulation(& mut self){

    }
    */
    /// This function updates the point quanities, that is to say 
    /// * it moves the point by computing the lagrangian displacement
    /// * it computes the associated total velocity field 
    /// * it computes variations of the effective temperature, and logarithm of surface gravity. 
    /// ### Returns:
    /// * This function returns a [Result] where [Ok] (()) if everything went well or [Err] if something went wrong in the callable methods. 
    pub fn update_point_quantities(
        & mut self,
        parameters:& PulstarConfig,
        k:&Coordinates,
        tar_collections:&[Option<TARCollection>],
        index:usize)-> Result<(),MathErrors>{
        
        let coordinates = point_coords_to_spherical(self.triangles.points[index].coords);
        let theta = coordinates.y;//This weird way of expressing it is because the nalgebra crate extracts the second and third coordinate of a vector as y and z.
        let phi = coordinates.z;
        
        displace_point(
            parameters, 
            & mut self.triangles.points[index],
            theta,
            phi,
            tar_collections)?;
        let g0 = 1.0e1f64.powf(self.log_g_of_points[index]);
        let temperature_0 = self.effective_temp_of_points[index];
        let v_tot = observed_pulsation_velocity(parameters, theta, phi,k,tar_collections).unwrap()+project_vrot(parameters, theta, phi, k);
        let (t_eff,log_g) = local_surface_temperature_logg(parameters, theta,phi, g0, temperature_0, tar_collections);
        self.log_g_of_points[index] = log_g;
        self.effective_temp_of_points[index] = t_eff;
        self.velocity_of_points[index]=v_tot;
        Ok(())
    }   

}

/// This function computes the lagrangian displacement in each point of the surface and moves the coordinates of the point. 
fn displace_point(
    parameters:& PulstarConfig,
    point:& mut Point,
    theta:f64,
    phi:f64,
    //mode:&PulsationMode,
    tar_collections:&[Option<TARCollection>]
    )->Result<(),MathErrors>{
    // compute theta and phi 
    let mut total_p_ds = Coordinates::Spherical(Vector3::new(0.0,0.0,0.0));
    // Compute total displacement
    for (index,mode) in parameters.mode_data.iter().enumerate(){
        let radial_amplitude = ampl_r(mode);
        let tangential_amplitude = ampl_t(mode);
        let pulsation_displacement = displacement(
            mode,
            theta,
            phi,
            radial_amplitude,
            tangential_amplitude,
            mode.get_spin_parameter(parameters),
            &tar_collections[index]
        )?;
        total_p_ds += pulsation_displacement;
    }
    //transform displacement into cartesian
    total_p_ds = total_p_ds.transform(theta, phi);

    // displace each coordinate
    if let Coordinates::Cartesian(displacement_coords) = total_p_ds {
        point.coords = point.coords + displacement_coords;
        Ok(())
    }else{Err(MathErrors::DifferentVectorBase)}
}


