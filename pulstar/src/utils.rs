//use std::any::Any;
//use std::hint::assert_unchecked;
use temp_name_lib::type_def::{Eigenfunctions,Config};
#[derive(Debug,PartialEq)]
pub struct PulstarConfig{
    pub time_pts_nmbr:u16,
    pub mode_config:Config,
    pub freqcycli:Vec<f64>,
    pub temperature_config:Vec<Eigenfunctions>,
    pub gravity_config:Vec<Eigenfunctions>,
    pub star_config:StarInfo,
    pub is_time_dependent:bool,
    pub suppress_pulse:bool,
    pub print_amplitude:bool,
}
#[derive(Debug,PartialEq)]
pub struct StarInfo{
    pub mass:f64,
    pub radius:f64,
    pub effective_temperature:f64,
    pub rotation_velocity:f64,
    pub inclination_angle:i16,
}

mod parse_value;
mod parse_string;
pub mod parse_file;

pub mod print_info;

#[derive(Debug,Clone,Copy)]
enum InputKind{
    F64(f64),
    U16(u16),
    I16(i16),
    BOOL(bool),
}

enum InputLines{
    //number of time points. Must be a u16
    TimePoints([InputKind;1]),
    
    //number of modes, entry must be a u16
    NumberOfModes([InputKind;1]),
    
    //l, m, amplitude r/r0, k, freq, phase offset
    //(u16,i16,f64,f64,f64,f64);
    ModeInfo([InputKind;6]),

    //Equatorial Rotation velocity and Inclination angle
    // f64,u16
    RotationInclination([InputKind;2]),
    
    //Factor ΔT/T_0  phase difference ΔT/T0
    //(f64,f64)
    ModeTemperature([InputKind;2]),

    //Factor Δg/g_0  phase difference Δg/g0
    //(f64,f64)
    ModeGravity([InputKind;2]),

    // Mass/Mass_sun, Radius/Radius_Sun, T_eff
    //(f64,f64,f64)
    StarInfo([InputKind;3]),

    //bool Expresses whether it will be time dependent (true) or not(false)
    IsTimeDependent([InputKind;1]),

    //bool Supress the pulsational velocity
    SuppressPulseVel([InputKind;1]),

    //bool Compute and print the max velocity amplitude and relative disp vector
    PrintMaxVelAmplitude([InputKind;1]),
}

fn search_geq(vector:&[f64],key:f64)-> usize {
    //non empty vector
    let first_element =vector.get(0);
    if let Some(first_value) = first_element{
        if key <= *first_value {0usize}
        else{
            let size = vector.len();
            let mut top = size - 1;
            let mut bottom = 0;
            let mut middle = bottom + (top - bottom) / 2;
            if key > vector[top]{panic!("key value is too large")}

            while bottom<top{
                if vector[middle] <key{
                    bottom = middle + 1;
                } else {
                    top = middle;
                }
                middle = bottom + (top - bottom)/2;
            }
            top
        }
    }
    else{ panic!("empty vector!")}
}