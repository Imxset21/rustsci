/// The Airy functions Ai(x) and Bi(x) are defined by the integral representations:
/// Ai(x) = (1/\pi) \int_0^\infty \cos((1/3) t^3 + xt) dt
/// Bi(x) = (1/\pi) \int_0^\infty (e^(-(1/3) t^3 + xt) + \sin((1/3) t^3 + xt)) dt
/// For further information see Abramowitz & Stegun, Section 10.4.

use libc::{c_int, c_uint, c_double};
use gsl_sf;

#[link(name = "gsl")]
extern
{
    /// Compute the Airy function Ai(x) with an accuracy specified by mode.
    fn gsl_sf_airy_Ai_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Compute the Airy function Bi(x) with an accuracy specified by mode.
    fn gsl_sf_airy_Bi_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Compute a scaled version of the Airy function S_A(x) Ai(x).
    /// For x>0 the scaling factor S_A(x) is exp(+(2/3) x^(3/2)), 1 for x<0. 
    fn gsl_sf_airy_Ai_scaled_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Compute a scaled version of the Airy function S_B(x) Bi(x). For x>0 the
    /// scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
    fn gsl_sf_airy_Bi_scaled_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the Airy function derivative Ai'(x) with an
    /// accuracy specified by mode.
    fn gsl_sf_airy_Ai_deriv_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the Airy function derivative Bi'(x) with an
    /// accuracy specified by mode.
    fn gsl_sf_airy_Bi_deriv_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the scaled Airy function derivative S_A(x)
    /// Ai'(x). For x>0 the scaling factor S_A(x) is \exp(+(2/3) x^(3/2)), and
    /// is 1 for x<0.
    fn gsl_sf_airy_Ai_deriv_scaled_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the scaled Airy function derivative S_B(x)
    /// Bi'(x). For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is
    /// 1 for x<0.
    fn gsl_sf_airy_Bi_deriv_scaled_e(
        x: c_double,
        mode: c_int,
        result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the location of the s-th zero of the Airy
    /// function Ai(x).
    fn gsl_sf_airy_zero_Ai_e(s: c_uint, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the location of the s-th zero of the Airy
    /// function Bi(x).
    fn gsl_sf_airy_zero_Bi_e(s: c_uint, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the location of the s-th zero of the Airy
    /// function derivative Ai'(x).
    fn gsl_sf_airy_zero_Ai_deriv_e(
        s: c_uint, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// These routines compute the location of the s-th zero of the Airy
    /// function derivative Bi'(x).
    fn gsl_sf_airy_zero_Bi_deriv_e(
        s: c_uint, result: *mut gsl_sf::gsl_sf_result) -> c_int;
}

/// Compute the Airy function Ai(x), returning the value and estimated error
pub fn airy_ai(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Ai_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the Airy function Bi(x), returning the value and estimated error
pub fn airy_bi(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Bi_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute a scaled version of the Airy function S_A(x) Ai(x).
/// For x>0 the scaling factor S_A(x) is exp(+(2/3) x^(3/2)), 1 for x<0. 
pub fn airy_ai_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Ai_scaled_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute a scaled version of the Airy function S_B(x) Bi(x).
/// For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
pub fn airy_bi_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Bi_scaled_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the Airy function derivative Ai'(x)
pub fn airy_ai_deriv(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Ai_deriv_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy derivative function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the Airy function derivative Bi'(x)
pub fn airy_bi_deriv(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Bi_deriv_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy derivative function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the scaled Airy function derivative S_A(x) Ai'(x). For x>0 the
/// scaling factor S_A(x) is exp(+(2/3) x^(3/2)), and is 1 for x<0.
pub fn airy_ai_deriv_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Ai_deriv_scaled_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy derivative function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the scaled Airy function derivative S_B(x) Bi'(x). For x>0 the
/// scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
pub fn airy_bi_deriv_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_Bi_deriv_scaled_e(x, gsl_sf::GSL_PREC_DOUBLE, &mut s) != 0
        {
            panic!("Airy derivative function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the location of the n-th zero of the Airy function Ai(x).
pub fn airy_zero_ai(n: u32) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_zero_Ai_e(n, &mut s) != 0
        {
            panic!("Airy zero function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the location of the n-th zero of the Airy function Bi(x).
pub fn airy_zero_bi(n: u32) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_zero_Bi_e(n, &mut s) != 0
        {
            panic!("Airy zero function calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Compute the location of the n-th zero of the Airy function Ai(x).
pub fn airy_zero_ai_deriv(n: u32) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_zero_Ai_deriv_e(n, &mut s) != 0
        {
            panic!("Airy zero function calculation failed");
        }
    }
    (s.val, s.err)
}

/// Compute the location of the n-th zero of the Airy function Bi(x).
pub fn airy_zero_bi_deriv(n: u32) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_airy_zero_Bi_deriv_e(n, &mut s) != 0
        {
            panic!("Airy zero function calculation failed");
        }
    }
    (s.val, s.err)    
}
