use super::InputKind;
pub fn parse_from_cell (s: &str)->InputKind{
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

    panic!("Failed to parse {} as a recognized primitive value",s);
}



#[cfg(test)]
mod tests{

    use super::{parse_from_cell,InputKind};

    #[test]
    #[ignore = "allready passed"]
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
    #[ignore = "allready passed"]
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
    #[ignore = "allready passed"]
    fn parsing_cell_value_u16(){
        let s1 = String::from("1");
        match parse_from_cell(&s1){
            InputKind::U16(value)=>{assert_eq!(1u16,value);}
            _=>{panic!("not u16")}
        } 
    }

    #[test]
    #[ignore = "allready passed"]
    fn parsing_cell_value_i16(){
        let s2 = String::from("-1");
        match parse_from_cell(&s2){
            InputKind::I16(value)=>{assert_eq!(-1i16,value);}
            _=>{panic!("not i16")}
        }
    }


}