pub use std::f64::consts::PI as PI;

/// G constant in m^3 s^{-2} kg^{-1}
pub const GRAVCONSTANT:f64 = 6.67259e-11;       // SI-units (m^3/s^2/kg) */
/// Mass of the sun in Kg
pub const MASSSUN:f64 = 1.9891e30;              // SI-units (kg)         */
/// Solar radius in m
pub const RADIUSSUN:f64 = 6.9599e8;             // SI-units (m)          */
/// Speed of light in m/s
pub const CLIGHT:f64 = 299792458.0;             // SI-units (m/s)        */

//pub const PI:f64 = 3.14159265358979;
/// Conversion constant from degrees to radians
pub const DEG2RAD:f64 = PI/180.0;             // Conversion from degrees to radians */
/// Conversion constant from radians to degrees
pub const RAD2DEG:f64 = 180.0/PI;             // Conversion from radians to degrees */
/// Number of seconds in a day.
pub const SEC_IN_DAY:f64 = 86400.0;             // Number of seconds in a day         */
/// Conversion constant from cycles per day to radians per second
pub const CYCLI2RAD:f64 = 2.0*PI/SEC_IN_DAY;  // Conversion from cycli/day to rad/s */


pub const N_FLUX_POINTS:u16=10000; // number of points in one flux profile
pub const MAX_N_TIMES:u16 = 3000; //max number of time points 
