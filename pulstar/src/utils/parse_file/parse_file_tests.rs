use super::*;
use std::{fs,io};
#[test]
fn open_file_and_get_first_line(){
    let path = String::from("dummy_file_lines.txt");
    match fs::write(&path,"First Line"){
        Ok(_)=>{println!("Created dummy file")}
        Err(e)=>{return;}
    }
    let mut ssstring =get_line(&path).unwrap().map(|l| l.expect("error while reading"));
    assert_eq!(ssstring.next(),Some(String::from("First Line")));
    match fs::remove_file(path){
        Ok(_)=>{println!("Removing dummy file")}
        Err(e)=>{eprintln!("couldn't remove dummy file due to {}",e)}
    }
}