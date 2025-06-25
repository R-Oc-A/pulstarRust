use super::parse_value::parse_from_cell;
use crate::utils::{InputLines,InputKind};
use temp_name_lib::type_def;
pub fn parse_from_string(s: &str, input_line_id:u16 )->Option<InputLines>{
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
                            if value < type_def::MAX_N_TIMES {
                            Some(InputLines::TimePoints([InputKind::U16(value)]))
                            }
                            else{panic!("number of timepoints must be less than {}",type_def::MAX_N_TIMES)}
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
                            if value < type_def::MAX_N_MODES{
                            Some(InputLines::NumberOfModes([InputKind::U16(value)]))}
                            else{panic!("number of modes must be less than {}",type_def::MAX_N_MODES)}
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
                    let mut l_check:u16 = 0;
                    match input_vec[0]{
                        InputKind::U16(value)=>{ 
                            l_check = value;
                            aux_vec.push(InputKind::U16(value));
                        }
                        _=>{panic!("Corrupt input file, error in mode info expected u16 for l.")}
                    }
                    match input_vec[1]{
                        InputKind::I16(value)=>{ 
                            if value < l_check as i16 {aux_vec.push(InputKind::I16(value));}
                            else{panic!("Corrupt input file, error in mode info, m value greater than l.")}
                        }
                        InputKind::U16(value)=>{
                            if value > (std::i16::MAX as u16){panic!("Corrupt input file, error in mode info i16 out of bounds")}
                            else if value > l_check {panic!("Corrupt input file, error in mode info, m value greater than l.")}
                            else {aux_vec.push(InputKind::I16(value as i16))}
                        }
                        _=>{panic!("Corrupt input file, error in mode info expected i16 for m")}
                    }
                    match input_vec[2]{
                        InputKind::F64(value)=>{ aux_vec.push(InputKind::F64(value));}
                        _=>{panic!("Corrupt input file, error in mode info expected f64 for amplitude r/r0")}
                    }
                    match input_vec[3]{
                        InputKind::F64(value)=>{ 
                            if l_check ==0 || value != 0.0{ panic!("Corrupt input file, error in mode info k value !=0 for l=0.")}
                            else{aux_vec.push(InputKind::F64(value));}
                        }
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
                    Some(InputLines::ModeInfo([aux_vec[0],
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
                    Some(InputLines::RotationInclination([aux_vec[0],
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
                Some(InputLines::ModeTemperature([aux_vec[0],
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
                
                Some(InputLines::ModeGravity([aux_vec[0],
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
                Some(InputLines::StarInfo([aux_vec[0],
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
                Some(InputLines::IsTimeDependent([aux_vec[0]]))
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
                Some(InputLines::SuppressPulseVel([aux_vec[0]]))
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
                
                Some(InputLines::PrintMaxVelAmplitude([aux_vec[0]]))
            }

            _=>{panic!("There's not an input line with that id number")}
            }
            }
}


#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    #[ignore = "allready passed"]
    fn parsing_time_points(){
        let s = String::from("     15");
        if let Some(value) = parse_from_string(&s, 1){
            match value{
                InputLines::TimePoints(value2)=>{
                    match value2[0]{
                        InputKind::U16(value3)=>{assert_eq!(value3,15u16,
                            "it went almost completely allright except that it received a '{}' instead of 15",value3)}
                        _=>{panic!("read an input but it wasn't 15u16")}
                    }
                }
                _=>{panic!("received an incorrect input line id")}
            }
        }else{panic!("It couldn't even get an input line")};
    }
    
    #[test]
    #[ignore = "allready passed"]
    fn parsing_number_of_modes(){
        let s = String::from("     1");
        if let Some(value) = parse_from_string(&s, 2){
            match value{
                InputLines::NumberOfModes(value2)=>{
                    match value2[0]{
                        InputKind::U16(value3)=>{assert_eq!(value3,1u16,
                            "it went almost completely allright except that it received a '{}' instead of 1",value3)}
                        _=>{panic!("read an input but it wasn't 1u16, rather a {:#?}",value2[0])}
                    }
                }
                _=>{panic!("received an incorrect input line id")}
            }
        }else{panic!("It couldn't even get an input line")};
    }

    #[test]
    #[ignore = "allready passed"]
    fn parsing_mode_info(){
        let s = String::from("  4     1       0.024              0.05           6.74            0.00");
        if let Some(value) = parse_from_string(&s, 3){
            match value{
                InputLines::ModeInfo(value2)=>{
                    match value2[0]{
                        InputKind::U16(value3)=>{assert_eq!(value3,4u16,
                            "it went almost completely fine except that it received a '{}' instead of  4", value3)}
                        _=>{panic!("read an input for l but it wasn't a u16, rather a {:#?}",value2[0])}
                    }
                    match value2[1]{
                        InputKind::I16(value3)=>{assert_eq!(value3,1i16,
                            "it went almost completely fine except that it received a '{}' instead of a 1",value3)}
                        _=>{panic!("read an input for m but it wasn't a i16, rather a {:#?}",value2[0])}
                    }
                    match value2[2]{
                        InputKind::F64(value3)=>{assert_eq!(value3,0.024,
                            "it went almost completely fine except that it received a '{}' instead of a 0.024",value3)}
                        _=>{panic!("read an input for ampl delta r/r0 but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[3]{
                        InputKind::F64(value3)=>{assert_eq!(value3,0.05,
                            "it went almost completely fine except that it received a '{}' instead of a 0.05",value3)}
                        _=>{panic!("read an input for K value but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[4]{
                        InputKind::F64(value3)=>{assert_eq!(value3,6.74,
                            "it went almost completely fine except that it received a '{}' instead of a 6.74",value3)}
                        _=>{panic!("read an input for frequency cycle (c/d) but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[5]{
                        InputKind::F64(value3)=>{assert_eq!(value3,0.00,
                            "it went almost completely fine except that it received a '{}' instead of a 0.00",value3)}
                        _=>{panic!("read an input for phase offset but it wasn't a i16, rather a {:#?}",value2[0])}
                    }
                }
                _=>{panic!("received an incorrect input line id")}
            }
        }else{
            panic!("It couldn't even get an input line")
        }
    }


    #[test]
    #[ignore = "allready passed"]
    fn parsing_rotvel_inclination(){
        let s = String::from("             20.0                                        45");
        if let Some(value) = parse_from_string(&s, 4){
            match value{
                InputLines::RotationInclination(value2)=>{
                    match value2[0]{
                        InputKind::F64(value3)=>{assert_eq!(value3,20.0,
                            "It received a '{}' instead of a 20.0",value3)}
                        _=>{panic!("read an input for rotation but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[1]{
                        InputKind::I16(value3)=>{assert_eq!(value3,45i16,
                            "It received a '{}' instead of a 45",value3)}
                        _=>{panic!("read an input for inclination angle but it wasn't a i16, rather a {:#?}",value2[0])}
                    }
                }
                _=>{panic!("Received an incorrect input line id")}
            }
        }else { panic!("It couldn't even get an input line")}
    }

    #[test]
    #[ignore = "allready passed"]
    fn parsing_mode_temperature(){
        let s = String::from("      2.62                                        180.0");
        if  let Some(value) = parse_from_string(&s, 5){
            match value{
                InputLines::ModeTemperature(value2)=>{
                    match value2[0]{
                        InputKind::F64(value3)=>{assert_eq!(value3,2.62,
                            "It received a '{}' instead of a 2.62",value3)}
                        _=>{panic!("read an input for factor delta T/T0 but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[1]{
                        InputKind::F64(value3)=>{assert_eq!(value3,180.0,
                            "It received a '{}' instead of a 180.0",value3)}
                        _=>{panic!("read an input for phase difference delta T/T0 but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                }
                _=>{panic!("Received an incorrect input line id")}
            }
        }else { panic!("It couldn't even get an input line")}
    }
    
    #[test]
    #[ignore = "allready passed"]
    fn parsing_mode_gravity(){
        let s = String::from("      10.0                                        34.0");
        if  let Some(value) = parse_from_string(&s, 6){
            match value{
                InputLines::ModeGravity(value2)=>{
                    match value2[0]{
                        InputKind::F64(value3)=>{assert_eq!(value3,10.0,
                            "It received a '{}' instead of a 10.0",value3)}
                        _=>{panic!("read an input for factor delta g/g0 but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[1]{
                        InputKind::F64(value3)=>{assert_eq!(value3,34.0,
                            "It received a '{}' instead of a 34.0",value3)}
                        _=>{panic!("read an input for phase difference delta g/g0 but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                }
                _=>{panic!("Received an incorrect input line id")}
            }
        }else { panic!("It couldn't even get an input line")}
    }

    #[test]
    #[ignore = "allready passed"]
    fn parsing_star_info(){
        let s = String::from("    10.0            6.93                   22642.0");
        if let Some(value) = parse_from_string(&s, 7){
            match value{
                InputLines::StarInfo(value2)=>{
                    match value2[0]{
                        InputKind::F64(value3)=>{assert_eq!(value3,10.0,
                            "It received a '{}' instead of a 10.0",value3)}
                        _=>{panic!("read an input for the star mass but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[1]{
                        InputKind::F64(value3)=>{assert_eq!(value3,6.93,
                            "It received a '{}' instead of a 6.93",value3)}
                        _=>{panic!("read an input for Radius/Radius_sun but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                    match value2[2]{
                        InputKind::F64(value3)=>{assert_eq!(value3,22642.0,
                            "It received a '{}' instead of a 22642.0",value3)}
                        _=>{panic!("read an input for effective temperature but it wasn't a f64, rather a {:#?}",value2[0])}
                    }
                }
                _=>{panic!("Received an incorrect input line id")}
            }
        }else { panic!("It couldn't even get an input line")}
    }
    
    #[test]
    #[ignore = "allready passed"]
    fn parsing_is_time_dependent(){
        let s = String::from("                       1");
        if let Some(value) = parse_from_string(&s, 8){
            match value{
                InputLines::IsTimeDependent(value2)=>{
                    match value2[0]{
                        InputKind::U16(value3)=>{assert_eq!(value3,1u16,
                            "it went almost completely allright except that it received a '{}' instead of 1",value3)}
                        _=>{panic!("read an input but it wasn't u16")}
                    }
                }
                _=>{panic!("received an incorrect input line id")}
            }
        }else{panic!("It couldn't even get an input line")};
    }
    
    #[test]
    #[ignore = "allready passed"]
    fn parsing_suppress_pulse(){
        let s = String::from("                       0");
        if let Some(value) = parse_from_string(&s, 9){
            match value{
                InputLines::SuppressPulseVel(value2)=>{
                    match value2[0]{
                        InputKind::U16(value3)=>{assert_eq!(value3,0u16,
                            "it went almost completely allright except that it received a '{}' instead of 0",value3)}
                        _=>{panic!("read an input but it wasn't u16")}
                    }
                }
                _=>{panic!("received an incorrect input line id")}
            }
        }else{panic!("It couldn't even get an input line")};
    }
    
    #[test]
    #[ignore = "allready passed"]
    fn parsing_print_vel(){
        let s = String::from("                       0");
        if let Some(value) = parse_from_string(&s, 10){
            match value{
                InputLines::PrintMaxVelAmplitude(value2)=>{
                    match value2[0]{
                        InputKind::U16(value3)=>{assert_eq!(value3,0u16,
                            "it went almost completely allright except that it received a '{}' instead of 0",value3)}
                        _=>{panic!("read an input but it wasn't u16")}
                    }
                }
                _=>{panic!("received an incorrect input line id")}
            }
        }else{panic!("It couldn't even get an input line")};
    }
}