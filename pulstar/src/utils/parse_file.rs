use std::io::{self, BufRead, Lines};
use std::fs::{self,File};
use std::path::Path;
use super::{Config,Eigenfunctions};
use super::parse_string::parse_from_string;
use super::{InputKind,InputLines};
use super::{StarInfo,PulstarConfig};

///Returns a result with Ok(iterator containing a line);
fn get_line<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn parse_from_file(file_name:&str)->PulstarConfig{    
    //This is in case reading file from buffer needs more time to set up.
    let parameter_file_contents = fs::read_to_string(file_name)
        .expect("Should have been able to read the file");

    //----List of parameters-----
    let mut counter:u16 =1;
    let mut TimePoints:u16=0;
    let mut n_modes:u16=0;
    let mut freqcycli:Vec<f64>=Vec::new();
    let mut l: Vec<u16>=Vec::new();
    let mut m: Vec<i16>=Vec::new();
    let mut rel_deltar:Vec<f64>=Vec::new();
    let mut k:Vec<f64>=Vec::new();
    let mut phase:Vec<f64>=Vec::new();
    let mut mass:f64=0.0;
    let mut radius:f64=0.0;
    let mut inclination_angle:i16=0;
    let mut rotation_velocity:f64=0.0;
    let mut effective_temperature:f64 = 0.0;
    let mut ModeTemperature:Vec<Eigenfunctions>=Vec::new();
    let mut ModeGravity:Vec<Eigenfunctions>=Vec::new();
    let mut IsTimeDependent:bool=true;
    let mut suppress_pulse:bool = false;
    let mut print_amplitude:bool = false;

    let mut n_modes_ctr:u16=0;
    //for s in parameter_file_contents.lines(){
    for s in get_line(file_name).unwrap().map(|l| l.expect("error while reading") ){
        match counter{
        1u16=>{
        if let Some(item)=parse_from_string(&s,counter){
        match item{
            InputLines::TimePoints(value)=>{
                match value[0]{
                    InputKind::U16(value2)=>{
                        TimePoints = value2;
                        counter +=1;
                    }
                    _=>panic!("this should not occur!error in number of timepoints"),
                }
            }
            _=>{panic!("It was expecting number of time points")}
            };
        }
        }
        2u16=>{
        if let Some(item)=parse_from_string(&s,counter){
        match item{
            InputLines::NumberOfModes(value)=>{
                match value[0]{
                    InputKind::U16(value2)=>{
                        n_modes = value2;
                        counter +=1;
                    }
                    _=>{panic!("this should not occur! error on format of number of modes");}
                }
            }
            _=>{panic!("It was expecting number of time points")}
            };
        }
        }
        3u16=>{
        if let Some(item) = parse_from_string(&s, counter){
        match item{
            InputLines::ModeInfo(value)=>{
                match value[0]{
                    InputKind::U16(value2)=>{
                        l.push(value2); 
                    }
                    _=>{panic!("was expecting a u16 value for l and failed")}
                }
                match value[1]{
                    InputKind::I16(value2)=>{
                        m.push(value2); 
                    }
                    _=>{panic!("was expecting a i16 value for m and failed")}
                }
                match value[2]{
                    InputKind::F64(value2)=>{
                        rel_deltar.push(value2); 
                    }
                    _=>{panic!("was expecting a f64 value for delr/r0 and failed")}
                }
                match value[3]{
                    InputKind::F64(value2)=>{
                        k.push(value2);
                    }
                    _=>{panic!("was expecting a f64 value for K value")}
                }
                match value[4]{
                    InputKind::F64(value2)=>{
                        freqcycli.push(value2);
                    }
                    _=>{panic!("Was expecting a f64 value for freqcycli")}
                }
                match value[5]{
                    InputKind::F64(value2)=>{
                        if value2>0.0||value2<=1.0{
                        phase.push(value2);
                        n_modes_ctr +=1;}
                        else{panic!("phase offset of mode {} is not within bounds [0,1]",n_modes_ctr)}
                    }
                    _=>{panic!("Was expecting a f64 value for phase offset")}       
                }
                if n_modes_ctr == n_modes{
                    n_modes_ctr=0;
                    counter +=1;
                }
            }
            _=>{panic!("error parsing mode info l,m,amplitude,k,...")}
        } 
        }
        }
        4u16=>{
            if let Some(item)=parse_from_string(&s, counter){
                match item{
                    InputLines::RotationInclination(value)=>{
                        match value[0]{
                            InputKind::F64(value2)=>{
                                rotation_velocity= value2;
                            }
                            _=>{panic!("Was expecting f64 value for rotational velocity")}
                        }
                        match value[1]{
                            InputKind::I16(value2)=>{
                                inclination_angle = value2;
                            }
                            _=>{panic!("Was expecting an integer value for the inclination angle")}
                        }
                        counter +=1;
                    }
                    _=>{panic!("Something went wrong on rotational velocity and inclination angle")}
                }
            }
        }
        5u16=>{
            if let Some(item)=parse_from_string(&s, counter){
                match item{
                    InputLines::ModeTemperature(value)=>{
                        let mut amplitude=0.0;
                        let mut phase_dif= 0.0;
                        match value[0]{
                            InputKind::F64(value2)=>{
                                amplitude = value2;
                            }
                            _=>{panic!("Was expecting a f64 value for T/T0")}
                        }
                        match value[1]{
                            InputKind::F64(value2)=>{
                                phase_dif =value2;
                            }
                            _=>{panic!("Was expecting a f64 value for phase difference T/T0")}
                        }
                        ModeTemperature.push(Eigenfunctions { ampl: amplitude, phasedif: phase_dif});
                        n_modes_ctr += 1;
                        if n_modes_ctr == n_modes{
                            n_modes_ctr=0;
                            counter += 1;
                        }
                    }
                    _=>{panic!("Something went wrong on temperature values for mode {}",n_modes_ctr)}
                }
            }
        }
        6u16=>{
            if let Some(item)=parse_from_string(&s, counter){
                match item {
                    InputLines::ModeGravity(value)=>{
                        let mut amplitude =0.0;
                        let mut phase_dif= 0.0;
                        match value[0]{
                            InputKind::F64(value2)=>{
                                amplitude = value2;
                            }
                            _=>{panic!("Was expecting a f64 value for g/g0")}
                        }
                        match value[1]{
                            InputKind::F64(value2)=>{
                                phase_dif =value2;
                            }
                            _=>{panic!("Was expecting a f64 value for phase difference g/g0")}
                        }
                        ModeGravity.push(Eigenfunctions { ampl: amplitude, phasedif: phase_dif});
                        n_modes_ctr += 1;
                        if n_modes_ctr == n_modes{
                            n_modes_ctr=0;
                            counter += 1;
                        }
                    }
                    _=>{panic!("Something went wrong on gravity values for mode {}",n_modes_ctr)}
                }
            }
        }
        7u16=>{
            if let Some(item)=parse_from_string(&s, counter){
                match item {
                    InputLines::StarInfo(value)=>{
                        match value[0]{
                            InputKind::F64(value2)=>{mass = value2;}
                            _=>{panic!("error, expecting a f64 value for the star's mass but got something else")}
                        }
                        match value[1]{
                            InputKind::F64(value2)=>{radius = value2;}
                            _=>{panic!("error, expecting a f64 value for the star's radius but got something else")}
                        }
                        match value[2]{
                            InputKind::F64(value2)=>{effective_temperature = value2;}
                            _=>{panic!("error, expecting a f64 value for the star's but got something else")}
                        }
                    counter += 1;    
                    }
                    _=>{panic!("Something went wrong on reading star's information on mass, Radius and/or effective temperature.")}
                }
            }
        }
        8u16=>{
            if let Some(item) = parse_from_string(&s, counter){
                match item {
                    InputLines::IsTimeDependent(value)=>{
                        match value[0]{
                            InputKind::U16(value2)=>{
                                match value2{
                                    1 => {IsTimeDependent = true;}
                                    0 => {IsTimeDependent = false;}
                                    _ => {panic!("Value for flag is not correct (either 0 or 1)")}
                                }
                                counter += 1;
                            }
                            _=>{panic!("error, expecting a u16 value for time dependent flag")}
                        }
                    }
                    _=>{panic!("Something went wrong on time dependent flag for surface normal")}
                }
            }
        }
        9u16=>{
            if let Some(item) = parse_from_string(&s, counter){
                match item {
                    InputLines::SuppressPulseVel(value)=>{
                        match value[0]{
                            InputKind::U16(value2)=>{
                                match value2{
                                    1 => {suppress_pulse = true;}
                                    0 => {suppress_pulse = false;}
                                    _ => {panic!("Value for flag is not correct (either 0 or 1)")}
                                }
                                counter += 1;
                            }
                            _=>{panic!("error, expecting a u16 value for suppress pulse flag")}
                        }
                    }
                    _=>{panic!("Something went wrong on suppressing pulsation flag")}
                }
            }
        }
        10u16=>{
            if let Some(item) = parse_from_string(&s, counter){
                match item {
                    InputLines::PrintMaxVelAmplitude(value)=>{
                        match value[0]{
                            InputKind::U16(value2)=>{
                                match value2{
                                    1 => {print_amplitude= true;}
                                    0 => {print_amplitude= false;}
                                    _ => {panic!("Value for flag is not correct (either 0 or 1)")}
                                }
                                counter += 1;
                            }
                            _=>{panic!("error, expecting a u16 value for print flag")}
                        }
                    }
                    _=>{panic!("Something went wrong on print amplitude flag for surface normal")}
                }
            }
        }

        _=>{if counter == 11 {break;} else {panic!("undefined behavior!")}}        
        }// match counter

    }//for line iteration
    
    
    let mode_parameters= Config{
        n_modes: n_modes,
        l: l,
        m: m,
        rel_deltar:rel_deltar,
        k:k,
        phase:phase,
    };//config
    
    let new_star=StarInfo{
        mass: mass,
        radius: radius,
        effective_temperature: effective_temperature,
        rotation_velocity: rotation_velocity,
        inclination_angle: inclination_angle,
    };

    PulstarConfig { 
         time_pts_nmbr: TimePoints,
         mode_config: mode_parameters,
         freqcycli: freqcycli, 
         temperature_config: ModeTemperature,
         gravity_config: ModeGravity, 
         star_config: new_star,
         is_time_dependent: IsTimeDependent,
         suppress_pulse: suppress_pulse,
         print_amplitude: print_amplitude }

}

#[cfg(test)]
mod parse_file_tests;