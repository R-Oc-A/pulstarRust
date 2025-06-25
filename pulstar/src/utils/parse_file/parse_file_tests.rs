use super::*;
use std::{fs, vec};
#[test]
fn open_file_and_get_first_line(){
    let path = String::from("dummy_file_line.txt");
    match fs::write(&path,"First Line"){
        Ok(_)=>{println!("Created dummy file")}
        Err(e)=>{eprintln!("couldn't create file, error {}",e);
            return;}
    }
    let mut ssstring =get_line(&path).unwrap().map(|l| l.expect("error while reading"));
    assert_eq!(ssstring.next(),Some(String::from("First Line")));
    match fs::remove_file(path){
        Ok(_)=>{println!("Removing dummy file")}
        Err(e)=>{eprintln!("couldn't remove dummy file due to {}",e)}
    }
}

#[test]
fn open_file_and_reach_to_the_end(){
    let path = String::from("dummy_file_lines.txt");
    match fs::write(&path,"First line 
    2
    3
    4
    5
    Last line"){
        Ok(_)=>{println!("Dummy file created.")}
        Err(e)=>{eprintln!("couldn't create dummy file, error {} ",e);
                        return;}
    }
    let mut counter = 1;

    for item in get_line(&path).unwrap().map(|l| l.expect("error while reading")){
        if counter == 6 {
            assert_eq!(item.trim_start(),"Last line")
        }
        counter += 1;
    }
    match fs::remove_file(path){
        Ok(_)=>{println!("Removing dummy file.")}
        Err(e)=>{eprintln!("couldn't remove dummy file due to {}",e)}
    }
}

#[test]
fn open_and_read_input_file(){
    let path = String::from("dummy_input_file.txt");
    match fs::write(&path, "
# This is an input file for the program PULSTAR
# 
# Some rules: . between 2 input lines you can always put comment lines with `#'
#             . you cannot put comments after a number on the same line
#
# 
# Number of time phase points in time points file
# (time points should be expressed in days)
# -----------------------------------------------
     15


# Total number of modes
# ---------------------
     1

# Make a table with mode information. For each mode:
#
# l     m     ampl. delta r/r0    K value      freq.(c/d)    phase offset [0,1]
#--    ---    ----------------    -------      ----------    ------------------
  4     1       0.024              0.05           6.74            0.00


# equatorial rotational velocity (km/s)     inclination angle (integer degrees)
# -------------------------------------     -----------------------------------
             20.0                                        45


# Make a table with mode information. For each mode:
#
# factor delta T/T_0         phase difference delta T/T_0 (degrees, float)
# ------------------         ---------------------------------------------
      2.62                                        180.0



# Make a table with mode information. For each mode:
# 
# factor delta g/g_0         phase difference delta g/g_0 (degrees, float)
# ------------------         ----------------------------------------------
      10.0                                        34.0


# Mass/Mass_sun  Radius/Radius_sun  effective temperature (K)
# -------------  -----------------  -------------------------
    10.0            6.93                   22642.0


# A time dependent surface normal (0: no; 1: yes)
# -----------------------------------------------
                       1

# Artificially suppress the pulsational velocity field (0: no; 1: yes)
# --------------------------------------------------------------------
                       0


# Compute and Print the maximum length of velocity vector and relative
# displacement vector. (0: no; 1: yes)
# --------------------------------------------------------------------
                       1

    
"){
    Ok(_)=>{println!("Dummy input file created.")}
    Err(e)=>{eprintln!("couldn't create dummy input file")}
}

let mode_pars = Config{
    n_modes:1u16,
    l: vec![4u16],
    m:vec![1i16],
    rel_deltar:vec![0.024],
    k:vec![0.05],
    phase:vec![0.00],
};

let freqcycli = vec![6.74];
let mode_temperature = vec![Eigenfunctions{
ampl:2.62,
phasedif:180.0,
}];
let mode_gravity = vec![Eigenfunctions{
    ampl:10.0,
    phasedif:34.0,
}];

let new_star= StarInfo{
    mass:10.0,
    radius: 6.93,
    effective_temperature:22642.0,
    rotation_velocity: 20.0,
    inclination_angle: 45i16
};

let time_points = 15u16;
let is_time_dependent:bool=true;
let suppress_pulse:bool = false;
let print_amplitude:bool = true;


let build_config = PulstarConfig{
    time_pts_nmbr: time_points,
    mode_config: mode_pars,
    star_config: new_star,
    freqcycli: freqcycli,
    temperature_config: mode_temperature,
    gravity_config:mode_gravity,
    is_time_dependent:is_time_dependent,
    suppress_pulse:suppress_pulse,
    print_amplitude:print_amplitude,
};

assert_eq!(parse_from_file(&path),build_config);

match fs::remove_file(path){
    Ok(_)=>{println!("Dummy input file removed.")}
    Err(e)=>{eprintln!("couldn't remove file due to error {}",e);
        return;}
}
}

//should panic
#[test]
#[should_panic]
fn open_corrupt_input_file_incomplete(){

    let path = String::from("dummy_input_file_incomplete.txt");
    match fs::write(&path, "
# This is an input file for the program PULSTAR
# 
# Some rules: . between 2 input lines you can always put comment lines with `#'
#             . you cannot put comments after a number on the same line
#
# 
# Number of time phase points in time points file
# (time points should be expressed in days)
# -----------------------------------------------
     15

"){
    Ok(_)=>{println!("Dummy input file created.")}
    Err(e)=>{eprintln!("couldn't create dummy input file")}
}

let build_config =parse_from_file(&path);

match fs::remove_file(path){
    Ok(_)=>{println!("Dummy input file removed.")}
    Err(e)=>{eprintln!("couldn't remove file due to error {}",e);
        return;}
}
}

#[test]
#[should_panic]
fn open_corrupt_file_mismatch_with_nmodes(){

    let path = String::from("dummy_input_file_mismatch.txt");
    match fs::write(&path, "
# This is an input file for the program PULSTAR
# 
# Some rules: . between 2 input lines you can always put comment lines with `#'
#             . you cannot put comments after a number on the same line
#
# 
# Number of time phase points in time points file
# (time points should be expressed in days)
# -----------------------------------------------
     15


# Total number of modes
# ---------------------
     2 

# Make a table with mode information. For each mode:
#
# l     m     ampl. delta r/r0    K value      freq.(c/d)    phase offset [0,1]
#--    ---    ----------------    -------      ----------    ------------------
  4     1       0.024              0.05           6.74            0.00

# equatorial rotational velocity (km/s)     inclination angle (integer degrees)
# -------------------------------------     -----------------------------------
             20.0                                        45


# Make a table with mode information. For each mode:
#
# factor delta T/T_0         phase difference delta T/T_0 (degrees, float)
# ------------------         ---------------------------------------------
      2.62                                        180.0



# Make a table with mode information. For each mode:
# 
# factor delta g/g_0         phase difference delta g/g_0 (degrees, float)
# ------------------         ----------------------------------------------
      10.0                                        34.0


# Mass/Mass_sun  Radius/Radius_sun  effective temperature (K)
# -------------  -----------------  -------------------------
    10.0            6.93                   22642.0


# A time dependent surface normal (0: no; 1: yes)
# -----------------------------------------------
                       1

# Artificially suppress the pulsational velocity field (0: no; 1: yes)
# --------------------------------------------------------------------
                       0


# Compute and Print the maximum length of velocity vector and relative
# displacement vector. (0: no; 1: yes)
# --------------------------------------------------------------------
                       1

    
"){
    Ok(_)=>{println!("Dummy input file created.")}
    Err(e)=>{eprintln!("couldn't create dummy input file")}
}

let build_config = parse_from_file(&path); 

match fs::remove_file(path){
    Ok(_)=>{println!("Dummy input file removed.")}
    Err(e)=>{eprintln!("couldn't remove file due to error {}",e);
        return;}
}
}

#[test]
#[should_panic]
fn open_corrupt_file_wrongvalue(){

    let path = String::from("dummy_input_file_wrong_val.txt");
    match fs::write(&path, "
# This is an input file for the program PULSTAR
# 
# Some rules: . between 2 input lines you can always put comment lines with `#'
#             . you cannot put comments after a number on the same line
#
# 
# Number of time phase points in time points file
# (time points should be expressed in days)
# -----------------------------------------------
     15


# Total number of modes
# ---------------------
     1

# Make a table with mode information. For each mode:
#
# l     m     ampl. delta r/r0    K value      freq.(c/d)    phase offset [0,1]
#--    ---    ----------------    -------      ----------    ------------------
  4     1       0              0.05           6.74            0.00


# equatorial rotational velocity (km/s)     inclination angle (integer degrees)
# -------------------------------------     -----------------------------------
             20.0                                        45.0               


# Make a table with mode information. For each mode:
#
# factor delta T/T_0         phase difference delta T/T_0 (degrees, float)
# ------------------         ---------------------------------------------
      2.62                                        180.0



# Make a table with mode information. For each mode:
# 
# factor delta g/g_0         phase difference delta g/g_0 (degrees, float)
# ------------------         ----------------------------------------------
      10.0                                        34.0


# Mass/Mass_sun  Radius/Radius_sun  effective temperature (K)
# -------------  -----------------  -------------------------
    10.0            6.93                   22642.0


# A time dependent surface normal (0: no; 1: yes)
# -----------------------------------------------
                       1

# Artificially suppress the pulsational velocity field (0: no; 1: yes)
# --------------------------------------------------------------------
                       0


# Compute and Print the maximum length of velocity vector and relative
# displacement vector. (0: no; 1: yes)
# --------------------------------------------------------------------
                       1

    
"){
    Ok(_)=>{println!("Dummy input file created.")}
    Err(e)=>{eprintln!("couldn't create dummy input file")}
}

let build_config = parse_from_file(&path);

match fs::remove_file(path){
    Ok(_)=>{println!("Dummy input file removed.")}
    Err(e)=>{eprintln!("couldn't remove file due to error {}",e);
        return;}
}
}