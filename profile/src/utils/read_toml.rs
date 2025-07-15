use serde::{Serialize,Deserialize};
use std::fs;
use std::process::exit;
use toml;


#[derive(Deserialize,Debug,PartialEq)]
struct ConfigParams{
    wavelength_range: WavelengthRange,
    pulsation_velocity: PulsationVelocity,
    intensity_grids_id: IntensityGrids,
}

#[derive(Deserialize,Debug,PartialEq)]
struct WavelengthRange{
start: f64,
end: f64,
step: f64,
}

#[derive(Deserialize,Debug,PartialEq)]
struct PulsationVelocity{
    max_velocity: f64,
}

#[derive(Deserialize,Debug,PartialEq)]
struct IntensityGrids{
    temperatures:Vec<f64>,
    log_gravity: Vec<f64>,
    filenames:Vec<String>,
} 


fn load_config_params (path_to_file: &str)->ConfigParams{
    let contents = match fs::read_to_string(path_to_file){
        Ok(c)=>c,
        Err(_) => { panic!("Could not read file {}",path_to_file)}
    };
    let params = match toml::from_str(&contents){
        Ok(d) => d,
        Err(_) => { panic!("Unable to load data from {}",contents)}
    }; 

    params 
}

#[cfg(test)]
#[test]
fn compare_toml(){
    let dummy_data:ConfigParams = toml::from_str(
"# This will be the typical structure of \n
# the input file for the input file for the profile code.\n 
\n
# Set the wave length range. The units are in nm \n 
[wavelength_range]\n
start = 412.0\n
end = 413.0\n
step = 0.01\n
\n
\n
# Upper limit of velocity pulsation. Units are in km/s.\n
# This velocity is used to compute what wave_lenght block will be read\n
# from the intensity grids.\n
[pulsation_velocity]\n
max_velocity = 1.0e2\n
\n
# Information of the intensity grid files\n
[intensity_grids_id]\n
\n
# Temperatures are in Kelvin\n
temperatures = [21000.0,21000.0,24000.0,24000.0]\n
# Logarithm of surface gravity\n
log_gravity = [3.5,4.5,3.5,4.5]\n
# Path to  intensity grid files\n
filenames = [\n
\"grids/21000g35.txt\",\n
\"grids/21000g45.txt\",\n
\"grids/24000g35.txt\",\n
\"grids/24000g45.txt\"] "
    ).expect("Shouldn't happen");

let lambda_range = WavelengthRange{
    start:412.0,
    end:413.0,
    step: 0.01,
};
let pulse = PulsationVelocity{
    max_velocity : 1.0e2,
};
let temps = vec![21000.0,21000.0,24000.0,24000.0];
let grids = IntensityGrids{
    temperatures: temps,
    log_gravity: vec![3.5,4.5,3.5,4.5],
    filenames: vec![String::from("grids/21000g35.txt"),
    String::from("grids/21000g45.txt"),
    String::from("grids/24000g35.txt"),
    String::from("grids/24000g45.txt")]
};
let params = ConfigParams{
    wavelength_range: lambda_range,
    pulsation_velocity: pulse,
    intensity_grids_id: grids
};

assert_eq!(params,dummy_data)

}
