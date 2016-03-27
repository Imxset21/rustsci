/// Special function common definitions
use libc::{c_int, c_double};

pub const GSL_PREC_DOUBLE: i32 = 0i32;
pub const GSL_PREC_SINGLE: i32 = 1i32;
pub const GSL_PREC_APPROX: i32 = 2i32;

#[repr(C)]
pub struct gsl_sf_result_struct
{
    pub val: c_double,
    pub err: c_double
}
#[allow(non_camel_case_types)]
pub type gsl_sf_result = gsl_sf_result_struct;

#[repr(C)]
pub struct gsl_sf_result_e10_struct
{
    pub val: c_double,
    pub err: c_double,
    pub e10: c_int
}
#[allow(non_camel_case_types)]
pub type gsl_sf_result_e10 = gsl_sf_result_e10_struct;
