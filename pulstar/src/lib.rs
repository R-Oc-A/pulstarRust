
const GRAVCONSTANT:f64 = 6.67259e-11;       // SI-units (m^3/s^2/kg) */
const MASSSUN:f64 = 1.9891e30;              // SI-units (kg)         */
const RADIUSSUN:f64 = 6.9599e8;             // SI-units (m)          */
const CLIGHT:f64 = 299792458.0;             // SI-units (m/s)        */

const PI:f64 = 3.14159265358979;
const DEG2RAD:f64 = PI/180.0;             // Conversion from degrees to radians */
const RAD2DEG:f64 = 180.0/PI;             // Conversion from radians to degrees */
const SEC_IN_DAY:f64 = 86400.0;             // Number of seconds in a day         */
const CYCLI2RAD:f64 = 2.0*PI/SEC_IN_DAY;  // Conversion from cycli/day to rad/s */

pub mod utils;