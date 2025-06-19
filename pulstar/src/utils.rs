use std::any::Any;
use std::io::{self, BufRead};
use std::fs::File;
use std::path::Path;


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
                        InputKind::U16(value){ aux_vec.push(InputKind::U16((value)));}
                        _=>{panic!("Corrupt input file, error in mode info expected u16 for l")}
                    }
                    match input_vec[1]{
                        InputKind::I16(value){ aux_vec.push(InputKind::I16((value)));}
                        _=>{panic!("Corrupt input file, error in mode info expected i16 for m")}
                    }
                    match input_vec[2]{
                        InputKind::F64(value){ aux_vec.push(InputKind::F64((value)));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude r/r0")}
                    }
                    match input_vec[3]{
                        InputKind::F64(value){ aux_vec.push(InputKind::F64((value)));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude K")}
                    }
                    match input_vec[4]{
                        InputKind::F64(value){ aux_vec.push(InputKind::F64((value)));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude freq.(c/d)")}
                    }
                    match input_vec[5]{
                        InputKind::F64(value){ 
                            if value>1.0||value<0 {panic!("Corrupt input file,phase offset value outside of [0,1] bounds")}
                            aux_vec.push(InputKind::F64((value)));}
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
                    _=>{panic!("Corrupt input file, expected i16 value for inclination angle")}
                }
                    Some(InputLines::rotation_inclination([aux_vec[0],
                                                        aux_vec[1]]))
            }

            5u16=>{
                //Check if well defined
                if input_vec.len()!=2usize {panic!("Corrupt input file, error in temperature fields")}
                let aux_vec:Vec<InputKind>=Vec::new();
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
                let aux_vec:Vec<InputKind>=Vec::new();
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
                let aux_vec:Vec<InputKind>=Vec::new();
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
                Some(InputLines::is_time_dependent([input_vec[0]]))
            }

            9u16=>{
                //Check if well defined
                //-snip
                Some(InputLines::suppress_pulse_vel([input_vec[0]]))
            }
            
            10u16=>{
                //Check if well defined
                //-snip
                Some(InputLines::print_max_vel_amplitude([input_vec[0]]))
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
//cfg(Tests)
//Test that gets you a key value.

