use crate::{ParsingFromToml, PulstarConfig};

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
        Self { mode_data: parse_input_file::PulsationModeNoPhases::
            get_initial_phases(input_parameters.mode_data),
		star_data: input_parameters.star_data,
		time_points: input_parameters.time_points,
		mesh: input_parameters.mesh}

    }
}