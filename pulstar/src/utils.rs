use cdshealpix::nested::map::fits::error;
use temp_name_lib::utils::MathErrors;

use crate::{ParsingFromToml, PulsationMode, PulstarConfig};

pub mod write_grid_data;

pub mod print_info;

mod parse_input_file;

impl ParsingFromToml for PulstarConfig{
    /// This function is used to fill the parameters required for the pulstar program to run out of the toml configuration file.
    /// #### Arguments:
    /// * `path_to_file` - this is a string that indicates the path to the `profile_input.toml` file
    /// #### Returns:
    /// * new instance of the profile config structure.
    fn read_from_toml(path_to_file:&str)->Self {
        let input_parameters = 
        parse_input_file::InputParameters::read_from_toml(path_to_file);
        let mode_data = parse_input_file::PulsationModeNoPhases
            ::get_initial_phases(input_parameters.mode_data);
        check_valid_lm_values_all_modes(& mode_data).expect("There's an ill defined pulsation mode. Please verify your quantities.");
        Self { mode_data: mode_data,
		star_data: input_parameters.star_data,
		time_points: input_parameters.time_points,
		mesh: input_parameters.mesh}

    }
}


fn check_for_valid_lm_values(mode: &PulsationMode)-> Result<(),MathErrors>{
    let error_string = format!(
        "frequency = {},\n l = {}, \n m = {}\n rotation regime : {:#?}",mode.frequency,mode.l,mode.m,mode.rotation_effects
    );
    if mode.l < mode.m.abs() as u16{
                println!("Error in mode \n{}",error_string);
                println!("azimuthal order is such that |m|>l");
                Err(MathErrors::OutOfBounds)
    }else{
        match mode.rotation_effects{
            crate::RotationRegime::CentrifugalDeformation{coefficient_expansion:values}=>{
                if mode.l==0||mode.l==1{
                    if values[2]!=0.0{//Check if its index=2 or index=1
                        println!("Error in mode \n{}",error_string);
                        println!("coefficient expansion not possible, l-2<0, please set this coefficient as zero");
                        return Err(MathErrors::OutOfBounds)
                    }
                }else{
                    if mode.l-2 < mode.m.abs() as u16{
                        if values[2]!=0.0{
                            println!("Error in mode \n{}",error_string);
                            println!("Coefficient expansion not possible, azimuthal order is such that m>l-2");
                            println!("Please choose an azimuthal order that complies with the expansion");
                            println!("or set the last coefficient to zero");
                            return Err(MathErrors::OutOfBounds)}
                    }
                }
                Ok(())
            }
            _=>{Ok(())}
        }
    }
}

fn check_valid_lm_values_all_modes(modes:&[PulsationMode])->Result<(),MathErrors>{
    for mode in modes.iter(){
        let _ = check_for_valid_lm_values(mode)?;
    }
    Ok(())
}