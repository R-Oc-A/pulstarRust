//use std::any::Any;
//use std::hint::assert_unchecked;
use std::io::{self, BufRead, Lines};
use std::fs::{self,File};
use std::path::Path;
use temp_name_lib::type_def::Config;
use temp_name_lib::type_def::Eigenfunctions;

pub struct PulstarConfig{
    time_pts_nmbr:u16,
    mode_config:Config,
    freqcycli:Vec<f64>,
    temperature_config:Vec<Eigenfunctions>,
    gravity_config:Vec<Eigenfunctions>,
    star_config:StarInfo,
    is_time_dependent:bool,
    suppress_pulse:bool,
    print_amplitude:bool,
}
pub struct StarInfo{
    mass:f64,
    radius:f64,
    effective_temperature:f64,
    rotation_velocity:f64,
    inclination_angle:i16,
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


///Returns a result with Ok(iterator containing a line);
fn get_line<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

//This function fills a tupple of inputs
fn parse_from_string(s: &str, input_line_id:u16 )->Option<InputLines>{
    //skip if # or end of line
    let trimmedline = s.trim();
    if trimmedline.is_empty()||trimmedline.starts_with('#'){
        return None;//skip 
    }else{
        let input_vec:Vec<InputKind> =  trimmedline
        .split_whitespace()
        .map(|s| parse_from_cell(s))
        .collect();
        match input_line_id{
            1u16=>{
                //check if well defined
                if input_vec.len()!= 1{ panic!("Corrupt input file, error in number of time points field")}
                else{
                    match input_vec[0]{
                        InputKind::U16(value)=>{
                            Some(InputLines::time_points([InputKind::U16(value)]))
                        }
                        _=>{panic!("Corrupt input file, in number of time points field\n
                                    expected a u16 value and received something else")}
                    }
                }
            }

            2u16=>{
                //check if well defined
                if input_vec.len()!= 1{ panic!("Corrupt input file, error in number of modes field")}
                else{
                    match input_vec[0]{
                        InputKind::U16(value)=>{
                            Some(InputLines::number_of_modes([InputKind::U16(value)]))
                        }
                        _=>{panic!("Corrupt input file, in number of modes field\n
                                    expected a u16 value and received something else")}
                    }
                }
            }

            3u16=>{
                //Check if well defined
                if input_vec.len()!=6usize{panic!("Corrupt input file, error in mode info")}
                else{
                    let mut aux_vec:Vec<InputKind>=Vec::new();
                    match input_vec[0]{
                        InputKind::U16(value)=>{ aux_vec.push(InputKind::U16(value));}
                        _=>{panic!("Corrupt input file, error in mode info expected u16 for l")}
                    }
                    match input_vec[1]{
                        InputKind::I16(value)=>{ aux_vec.push(InputKind::I16(value));}
                        InputKind::U16(value)=>{
                            if value > (std::i16::MAX as u16){panic!("Corrupt input file, error in mode info i16 out of bounds")}
                            else {aux_vec.push(InputKind::I16(value as i16))}
                        }
                        _=>{panic!("Corrupt input file, error in mode info expected i16 for m")}
                    }
                    match input_vec[2]{
                        InputKind::F64(value)=>{ aux_vec.push(InputKind::F64(value));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude r/r0")}
                    }
                    match input_vec[3]{
                        InputKind::F64(value)=>{ aux_vec.push(InputKind::F64(value));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude K")}
                    }
                    match input_vec[4]{
                        InputKind::F64(value)=>{ aux_vec.push(InputKind::F64(value));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude freq.(c/d)")}
                    }
                    match input_vec[5]{
                        InputKind::F64(value)=>{ 
                            if value>1.0||value<0.0 {panic!("Corrupt input file,phase offset value outside of [0,1] bounds")}
                            aux_vec.push(InputKind::F64(value));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude phase offset")}
                    }
                    Some(InputLines::mode_info([aux_vec[0],
                                        aux_vec[1],
                                        aux_vec[2],
                                        aux_vec[3],
                                        aux_vec[4],
                                        aux_vec[5]]))
                }
            }
            
            4u16=>{
                //Check if well defined
                if input_vec.len()!=2 {panic!("Corrupt input file, error in rotation and inclination fields")}
                let mut aux_vec: Vec<InputKind> = Vec::new();
                match input_vec[0]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for equatorial rotational velocity")}
                }
                match input_vec[1]{
                    InputKind::I16(value)=>{
                        aux_vec.push(InputKind::I16(value));
                    }
                    InputKind::U16(value)=>{
                        if value > (std::i16::MAX as u16) {panic!("Corrupt input file, error inclination angle out of bounds, expected i16")}
                        else {aux_vec.push(InputKind::I16(value as i16))}
                    }
                    _=>{panic!("Corrupt input file, expected i16 value for inclination angle")}
                }
                    Some(InputLines::rotation_inclination([aux_vec[0],
                                                        aux_vec[1]]))
            }

            5u16=>{
                //Check if well defined
                if input_vec.len()!=2usize {panic!("Corrupt input file, error in temperature fields")}
                let mut aux_vec:Vec<InputKind>=Vec::new();
                match input_vec[0]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for factor delta T/T0")}
                }
                match input_vec[1]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for phase difference delta T/T0")}
                }
                Some(InputLines::mode_temperature([aux_vec[0],
                                                    aux_vec[1]]))
            }

            6u16=>{
                //Check if well defined
                if input_vec.len()!=2usize {panic!("Corrupt input file, error in gravity fields")}
                let mut aux_vec:Vec<InputKind>=Vec::new();
                match input_vec[0]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for factor delta g/g0")}
                }
                match input_vec[1]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for phase difference delta g/g0")}
                }
                
                Some(InputLines::mode_gravity([aux_vec[0],
                                                aux_vec[1]]))
            }
            
            7u16=>{
                //Check if well defined
                if input_vec.len()!= 3usize{panic!("Corrupt input file, error in Star information field")}
                let mut aux_vec:Vec<InputKind>=Vec::new();
                match input_vec[0]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for Mass/Mass_sun")}
                }
                match input_vec[1]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for Radius/Radius_sun")}
                }
                match input_vec[2]{
                    InputKind::F64(value)=>{
                        aux_vec.push(InputKind::F64(value));
                    }
                    _=>{panic!("Corrupt input file, expected f64 value for effective temperature")}
                }
                Some(InputLines::star_info([aux_vec[0],
                                            aux_vec[1],
                                            aux_vec[2]]))
            }

            8u16=>{
                //Check if well defined
                if input_vec.len()!=1usize {panic!("Corrupt input file, error in time independence flag field")}
                let mut aux_vec:Vec<InputKind>=Vec::new();
                match input_vec[0]{
                    InputKind::U16(value)=>{
                        if value > 1 {panic!("Corrupt input file, flag for time independence can only be 0 or 1")}
                        aux_vec.push(InputKind::U16(value));
                    }
                    _=>{panic!("Corrupt input file, expected u16 for time independence flag")}
                }
                Some(InputLines::is_time_dependent([aux_vec[0]]))
            }

            9u16=>{
                //Check if well defined
                if input_vec.len()!=1usize {panic!("Corrupt input file, error in flag for suppressing pulse velocity")}
                let mut aux_vec:Vec<InputKind>=Vec::new();
                match input_vec[0]{
                    InputKind::U16(value)=>{
                        if value > 1 {panic!("Corrupt input file, flag for suppressing pulse velocity can only be 0 or 1")}
                        aux_vec.push(InputKind::U16(value));
                    }
                    _=>{panic!("Corrupt input file, expected u16 for suppressing pulse velocity flag")}
                }
                Some(InputLines::suppress_pulse_vel([aux_vec[0]]))
            }
            
            10u16=>{
                //Check if well defined
                if input_vec.len()!=1usize {panic!("Corrupt input file, error in flag for print amplitude of velocity")}
                let mut aux_vec:Vec<InputKind>=Vec::new();
                match input_vec[0]{
                    InputKind::U16(value)=>{
                        if value > 1 {panic!("Corrupt input file, flag for print amplitude of velocity can only be 0 or 1")}
                        aux_vec.push(InputKind::U16(value));
                    }
                    _=>{panic!("Corrupt input file, expected u16 for print amplitude of velocity flag")}
                }
                
                Some(InputLines::print_max_vel_amplitude([aux_vec[0]]))
            }

            _=>{panic!("There's not an input line with that id number")}
            }
            }
}

#[derive(Debug,Clone,Copy)]
enum InputKind{
    F64(f64),
    U16(u16),
    I16(i16),
    USIZE(usize),
    BOOL(bool),
}

enum InputLines{
    //number of time points. Must be a u16
    time_points([InputKind;1]),
    
    //number of modes, entry must be a u16
    number_of_modes([InputKind;1]),
    
    //l, m, amplitude r/r0, k, freq, phase offset
    //(u16,i16,f64,f64,f64,f64);
    mode_info([InputKind;6]),

    //Equatorial Rotation velocity and Inclination angle
    // f64,u16
    rotation_inclination([InputKind;2]),
    
    //Factor ΔT/T_0  phase difference ΔT/T0
    //(f64,f64)
    mode_temperature([InputKind;2]),

    //Factor Δg/g_0  phase difference Δg/g0
    //(f64,f64)
    mode_gravity([InputKind;2]),

    // Mass/Mass_sun, Radius/Radius_Sun, T_eff
    //(f64,f64,f64)
    star_info([InputKind;3]),

    //bool Expresses whether it will be time dependent (true) or not(false)
    is_time_dependent([InputKind;1]),

    //bool Supress the pulsational velocity
    suppress_pulse_vel([InputKind;1]),

    //bool Compute and print the max velocity amplitude and relative disp vector
    print_max_vel_amplitude([InputKind;1]),
}

fn parse_from_cell (s: &str)->InputKind{
    // Try bool
    if let Ok(b) = s.parse::<bool>(){
        return InputKind::BOOL(b);
    }

    // Try f64
    if s.contains('.')||s.contains('e')||s.contains('E'){
        if let Ok(f) = s.parse::<f64>(){
            return InputKind::F64(f);
        }
    }

    // Try U16
    if let Ok(u_16) = s.parse::<u16>(){
        return InputKind::U16(u_16);
    }

    // Try i16
    if let Ok(i_16) = s.parse::<i16>(){
        return InputKind::I16(i_16);
    }
    
    if let Ok(u_size) = s.parse::<usize>(){
        return InputKind::USIZE(u_size);
    }

    panic!("Failed to parse {} as a recognized primitive value",s);
}

pub fn parse_from_file(file_name:&str)->PulstarConfig{    
    //This is in case reading file from buffer needs more time to set up.
    let parameter_file_contents = fs::read_to_string(file_name)
        .expect("Should have been able to read the file");

    //----List of parameters-----
    let mut counter:u16 =1;
    let mut time_points:u16=0;
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
    let mut mode_temperature:Vec<Eigenfunctions>=Vec::new();
    let mut mode_gravity:Vec<Eigenfunctions>=Vec::new();
    let mut is_time_dependent:bool=true;
    let mut suppress_pulse:bool = false;
    let mut print_amplitude:bool = false;

    let mut n_modes_ctr:u16=0;
    //for s in parameter_file_contents.lines(){
    for s in get_line(file_name).unwrap().map(|l| l.expect("error while reading") ){
        match counter{
        1u16=>{
        if let Some(item)=parse_from_string(&s,counter){
        match item{
            InputLines::time_points(value)=>{
                match value[0]{
                    InputKind::U16(value2)=>{
                        time_points = value2;
                        counter +=1;
                    }
                    _=>panic!("this should not occur!"),
                }
            }
            _=>{panic!("It was expecting number of time points")}
            };
        }
        }
        2u16=>{
        if let Some(item)=parse_from_string(&s,counter){
        match item{
            InputLines::time_points(value)=>{
                match value[0]{
                    InputKind::U16(value2)=>{
                        n_modes = value2;
                        counter +=1;
                    }
                    _=>{panic!("this should not occur! error on format of number of timepoints");}
                }
            }
            _=>{panic!("It was expecting number of time points")}
            };
        }
        }
        3u16=>{
        if let Some(item) = parse_from_string(&s, counter){
        match item{
            InputLines::mode_info(value)=>{
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
                    InputLines::rotation_inclination(value)=>{
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
                    InputLines::mode_temperature(value)=>{
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
                        mode_temperature.push(Eigenfunctions { ampl: amplitude, phasedif: phase_dif});
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
                    InputLines::mode_gravity(value)=>{
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
                        mode_gravity.push(Eigenfunctions { ampl: amplitude, phasedif: phase_dif});
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
                    InputLines::star_info(value)=>{
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
                    InputLines::is_time_dependent(value)=>{
                        match value[0]{
                            InputKind::U16(value2)=>{
                                match value2{
                                    1 => {is_time_dependent = true;}
                                    0 => {is_time_dependent = false;}
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
                    InputLines::suppress_pulse_vel(value)=>{
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
                    InputLines::print_max_vel_amplitude(value)=>{
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
         time_pts_nmbr: time_points,
         mode_config: mode_parameters,
         freqcycli: freqcycli, 
         temperature_config: mode_temperature,
         gravity_config: mode_gravity, 
         star_config: new_star,
         is_time_dependent: is_time_dependent,
         suppress_pulse: suppress_pulse,
         print_amplitude: print_amplitude }

}

                    
#[cfg(test)]
mod tests{
    use super::*;

    use crate::utils::parse_from_cell;

    #[test]
    fn parsing_cell_value_float(){
        let s1 = String::from("1.0");
        match parse_from_cell(&s1){
            InputKind::F64(value)=>{
                assert_eq!(value,1.0);
            }
            _=>{panic!("not a floating number")}
        }
    }

    #[test]
    fn parsing_cell_value_bool(){
        let s1=String::from("true");
        let s2 = String::from("false");
        match parse_from_cell(&s1){
            InputKind::BOOL(value)=>{assert_eq!(true,value);}
            _=>{panic!("not correct bool")}
        }
        match parse_from_cell(&s2){
            InputKind::BOOL(value)=>{assert_eq!(false,value);}
            _=>{panic!("not correct bool")}
        }
    }

    #[test]
    fn parsing_cell_value_u16(){
        let s1 = String::from("1");
        match parse_from_cell(&s1){
            InputKind::U16(value)=>{assert_eq!(1u16,value);}
            _=>{panic!("not u16")}
        } 
    }

    #[test]
    fn parsing_cell_value_i16(){
        let s1 = String::from("1");
        let s2 = String::from("-1");
        /*match parse_from_cell(&s1){
            InputKind::I16(value)=>{assert_eq!(1i16,value);}
            _=>{panic!("not i16")}
        }*/ //This test fails, then if i16 is expected, the u16 must be casted as i16 when storing it in the inputs.
        match parse_from_cell(&s2){
            InputKind::I16(value)=>{assert_eq!(-1i16,value);}
            _=>{panic!("not i16")}
        }
    }


}
