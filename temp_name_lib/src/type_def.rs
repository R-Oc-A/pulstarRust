use std::f64::consts::PI as PI;

pub const GRAVCONSTANT:f64 = 6.67259e-11;       // SI-units (m^3/s^2/kg) */
pub const MASSSUN:f64 = 1.9891e30;              // SI-units (kg)         */
pub const RADIUSSUN:f64 = 6.9599e8;             // SI-units (m)          */
pub const CLIGHT:f64 = 299792458.0;             // SI-units (m/s)        */

//pub const PI:f64 = 3.14159265358979;
pub const DEG2RAD:f64 = PI/180.0;             // Conversion from degrees to radians */
pub const RAD2DEG:f64 = 180.0/PI;             // Conversion from radians to degrees */
pub const SEC_IN_DAY:f64 = 86400.0;             // Number of seconds in a day         */
pub const CYCLI2RAD:f64 = 2.0*PI/SEC_IN_DAY;  // Conversion from cycli/day to rad/s */


pub const N_FLUX_POINTS:u16=10000; // number of points in one flux profile
pub const MAX_N_TIMES:u16 = 3000; //max number of time points 
