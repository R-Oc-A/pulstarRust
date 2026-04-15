
use crate::reference_frames::rotation_treatment::tar::TARCollection;

use super::{*,na};
use temp_name_lib::math_module::{spherical_harmonics::norm_factor::ylmnorm};





/// This enum binds the components of a 3D vector depending on the type of coordinates that are used
/// 
/// as of now there're the following:
/// * 'Spherical' - The component order of this variant is `(r,θ,φ)` 
/// * 'Cartesian' - The component order of this variant is `(x,y,z)`
#[derive(Debug,PartialEq,Clone,Copy)]
pub enum Coordinates{
    Spherical(na::Vector3<f64>),//<-[Ricardo:]  I'll be using the nalgebra crate as it's well suited for linear algebra operation on small fixed size arrays
    Cartesian(na::Vector3<f64>),
}

/// In this module there's a collection of methods and implementations useful to handle the coordinate type of coordinates. 
/// So far we can ad vectors, multiply by a scalar from the left, transform from one basis to the other.
/// get the unit vector pointing along the observer
/// project from one vector onto the other 
/// get the spherical r component of a vector
mod implementations_for_coordinates;

/// Compute the Lagrangian displacement vector in spherical coordinates
/// 
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// * `theta` - The colatitude coordinate θ in radians
/// * 'dtheta' - difference between two colatitud coordinates between two neighbouring [SurfaceCell]s, it has the same units as θ and it assumes a regular mesh of the sphere in the colatitude coordinate.
/// * `phi`   - azimuthal coordinate  in radians
/// * `radial_amplitude`     - amplitude in the radial direction 
/// * `tangential_amplitude` - amplitude in the tangential direction 
/// * `tar_functions` - An [Option] enum that has the following variants:
///     - [Some] variant that has binded a reference to a [TARCollection]
///     - [None] in case tar functions are not needed.
/// ### Returns:
/// This function can return an [Ok] or [Err] variants of [Result] that will have the following values binded to them:
/// * `Ok(Coordinates::Spherical)` - an Ok  variant that has binded the spherical components of the displacement vector in the`r,θ,φ` order.
/// * `Err(DivisionByZero)` - an Err variant that has binded the error produced if the colatitude  coordinate (theta) is too small.
pub fn displacement(
    mode: &PulsationMode,
    theta:f64,
    dtheta:f64,
    phi:f64,
    radial_amplitude:f64,
    tangential_amplitude:f64,
    tar_functions:&Option<TARCollection>)->Result<Coordinates,MathErrors>{
            match mode.rotation_effects{
                RotationRegime::NonRotating =>{rotation_treatment::non_rotating::non_rotating_displacement(
                    mode,
                    theta.sin(),
                    theta.cos(),
                    phi,
                    radial_amplitude,
                    tangential_amplitude,
                )},
                RotationRegime::PerturbativeCoriolis =>{rotation_treatment::non_rotating::non_rotating_displacement(
                    mode,
                    theta.sin(),
                    theta.cos(),
                    phi,
                    radial_amplitude,
                    tangential_amplitude,
                )},
                RotationRegime::Tar =>{
                    if let Some(tar_collection) = tar_functions{
                        rotation_treatment::tar::tar_displacement(
                        mode,
                        theta,
                        dtheta,
                        phi,
                        radial_amplitude,
                        tangential_amplitude,
                        tar_collection)}
                    else{panic!("hough's functions were not properly loaded")}
                },
                RotationRegime::CentrifugalDeformation =>{rotation_treatment::non_rotating::non_rotating_displacement(
                    mode,
                    theta.sin(),
                    theta.cos(),
                    phi,
                    radial_amplitude,
                    tangential_amplitude,
                )},
            }
    }

/// This module contains the functions to calculate the derivatives of the lagrangian displacement vector over 
/// the surface of a star using spherical coordinates.
mod displacement_derivatives;

/// This module contains the functions used to implement rotational effects into the modelling of the pulsation equations. 
/// Currently these are the included mechanisms:
/// * `Non rotating` - in a non rotating case treatment of rotation, the rotation effects are only included in a frequency shift. 
/// * `TAR` - The traditional approximation of rotation includes the effects by disregarding the horizontal components of the coriolis force.\
/// * `Perturbative coriolis` - 
/// * `Centrifugal deformation`
pub mod rotation_treatment;

/// This function computes the spherical components of the surface normal vector on a reference frame where the z-axis 
/// coincides with the rotation axis. A surface normal is a vector which stands locally perpendicular to the surfaces and
/// has a length equal to the area of the local surface cell. 
/// 
/// *WARNING* This function does NOT multiply each component with 
/// R_0^2 dθdφ. Users must do this themselves. 
/// 
/// ### Arguments:
/// * `parameters` - The data contained in [PulstarConfig], here you find the parameters that describe the pulsation modes and the star.
/// * `theta` - The colatitude angle θ in rads, must not be too small in order to avoid the poles.
/// * 'dtheta' - difference between two colatitud coordinates between two neighbouring [SurfaceCell]s, it has the same units as θ and it assumes a regular mesh of the sphere in the colatitude coordinate.
/// * `phi` - The azimuthal angle in rads
/// * `tar_functions` - An [Option] enum that has the following variants:
///     - [Some] variant that has binded a reference to a [TARCollection]
///     - [None] in case tar functions are not needed.
/// ### Returns: 
/// This function returns [Ok] or [Err] variants of [Result]
/// * Ok(`Coordinates::Spherical(surface_normal_coords)`) - Where surface_normal_coords  is a [na::Vector3]
/// * Err(DivisionByZero) - Where the error is pased to the calling function if the colatitude angle θ is too small. 
pub fn surface_normal(
parameters: &PulstarConfig,
theta: f64,
dtheta: f64,
phi: f64,
tar_collections:&[Option<TARCollection>],
)->Result<Coordinates,MathErrors>{
    let sintheta = theta.sin();
    let costheta = theta.cos();

    let mut total_p_ds = Coordinates::Spherical(na::Vector3::new(0.0,0.0,0.0));
	let mut total_dev1=0.0;
	let mut total_dev2=0.0;
	let mut total_dev3=0.0;
	let mut total_dev4=0.0;

    for (index,mode) in parameters.mode_data.iter().enumerate(){

        let radial_amplitude = ampl_r(mode);
        let tangential_amplitude = ampl_t(mode);
        
        let pulsation_displacement = displacement(
            mode,
            theta,
            dtheta, 
            phi, 
            radial_amplitude, 
            tangential_amplitude,
            &tar_collections[index])?;
        
        let drdtheta = displacement_derivatives::d_dr_rdtheta(
            mode,
            theta, 
            dtheta, 
            phi,
            &tar_collections[index]);
        
        let drdphi = displacement_derivatives::d_dr_rdphi(
            mode,
            theta,
            dtheta,
            phi,
            &tar_collections[index]);
        
        let dtdtheta = displacement_derivatives::d_dtheta_dtheta(
            mode,
            theta,
            dtheta,
            phi,
            &tar_collections[index]);
        
        let dpdphi = displacement_derivatives::d_dphi_dphi(
            mode,
            theta,
            dtheta,
            phi,
            &tar_collections[index])?;//<- the ? is necesary to pas to the calling function if the colatitude angle theta is too close to the poles.
        
        
        total_p_ds += pulsation_displacement;

        total_dev1 += drdtheta;
        total_dev2 += drdphi;
        total_dev3 += dtdtheta;
        total_dev4 += dpdphi;
    } 
    let mut surface_normal_coords = na::Vector3::new(0.0,0.0,0.0);
    if let Coordinates::Spherical(coords) = total_p_ds{//<- This is an idiomatic way to extract values of an enum variant in Rust
        let r_hat = 1.0;
        let theta_hat = coords[1]-total_dev1;
        let phi_hat = sintheta * coords[2] - total_dev2/sintheta;

        let normal = na::Vector3::new(r_hat,theta_hat,phi_hat);
        let length = sintheta * (1.0 + 2.0 * coords[0]
            + costheta/sintheta * coords[1]
            + total_dev3 + total_dev4);
        
        surface_normal_coords = normal * length;

    }
    Ok(Coordinates::Spherical(surface_normal_coords))

}


///Given the components of both the surface normal vector and the unit vector
///pointing towards the observer, this function gives the cosine of the angle
///between those two vectors. 
/// 
/// The z-axis coincides with the rotation axis.
/// 
/// ### Arguments:
/// * `surface_normal` - a spherical [Coordinates] vector normal to a surface cell that has as length the (normalized) area of the cell. See [surface_normal]
/// * `k` - a unit vector pointing towards the observer on a frame of reference where the z-axis coincides with the rotation axis.
/// * `theta` - the colatitude angle on the star in radians. This angle should not be too small. 
/// * `phi` - the azimuthal angle on the star in radians.
/// 
/// ### Returns:
/// * `cos_chi` - a `f64` value that is the cosine of the angle between `surface_normal` and `k`
pub fn cos_chi(
    surface_normal:&Coordinates,
    k: &Coordinates,
    theta: f64,
    phi:f64)->f64{
    match k {
        Coordinates::Cartesian(_)=>{
            let k_spherical = k.transform(theta, phi);
            k_spherical.project_vector(surface_normal).unwrap()
            /(k_spherical.vector_length()*surface_normal.vector_length())
        }
        Coordinates::Spherical(_)=>{
            k.project_vector(surface_normal).unwrap()
            /(k.vector_length()*surface_normal.vector_length())
        }
    }
}


/// This function calculates the amplitude of the relative radial displacement multiplied by the normalization factor `Y_l^m`
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// ### Returns:
/// * `radial_amplitude` - A `f64` value that contains the amplitude of relative radial displacement (thus without units) caused by the pulsations of a given mode. 
pub fn ampl_r(mode:&PulsationMode)->f64{
    match mode.rotation_effects{
        RotationRegime::NonRotating => {mode.rel_dr * ylmnorm(mode.l, mode.m)},
        RotationRegime::CentrifugalDeformation => {mode.rel_dr 
            * ylmnorm(mode.l,mode.m)},
        RotationRegime::PerturbativeCoriolis => {
            mode.rel_dr * ylmnorm(mode.l, mode.m)
        },
        RotationRegime::Tar => {mode.rel_dr}
    }
}

/// This function calculates the amplitude of the relative tangential displacement multiplied by the normalization factor `Y_l^m`
/// ### Arguments:
/// * `mode` - This is a struct that contains the parameters of a pulsation mode in the star. See [crate::PulstarConfig]
/// ### Returns:
/// * `tangential_amplitude` - amplitude in the tangential direction times the normalization factor  'Y_l^m' (see [temp_name_lib::math_module::spherical_harmonics::norm_factors])
pub fn ampl_t(mode:&PulsationMode)->f64{
    match mode.rotation_effects{
        RotationRegime::NonRotating => {mode.rel_dr * ylmnorm(mode.l, mode.m)*mode.k}
        RotationRegime::PerturbativeCoriolis => {mode.rel_dr * ylmnorm(mode.l, mode.m)*mode.k}
        RotationRegime::CentrifugalDeformation=> {mode.rel_dr * ylmnorm(mode.l, mode.m)*mode.k}
        RotationRegime::Tar => {mode.rel_dr*mode.k}
    }
}