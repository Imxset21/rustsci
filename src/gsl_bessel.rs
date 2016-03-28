// The routines described in this section compute the Cylindrical Bessel
// functions J_n(x), Y_n(x), Modified cylindrical Bessel functions I_n(x),
// K_n(x), Spherical Bessel functions j_l(x), y_l(x), and Modified Spherical
// Bessel functions i_l(x), k_l(x). For more information see Abramowitz &
// Stegun, Chapters 9 and 10.

use libc::{c_int, c_uint, c_double, size_t};
use gsl_sf;
use array::Array;
use array::Order;

#[link(name = "gsl")]
extern
{
    /// Regular Bessel Function J_0(x)
    fn gsl_sf_bessel_J0_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular Bessel Function J_1(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_J1_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular Bessel Function J_n(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_Jn_e(
        n: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular Bessel Function J_n(x),  nmin <= n <= nmax
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_Jn_array(
        nmin: c_int,
        nmax: c_int,
        x: c_double,
        result_array: *mut c_double) -> c_int;
    
    /// Irregular Bessel function Y_0(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_Y0_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular Bessel function Y_1(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_Y1_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular Bessel function Y_n(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_Yn_e(
        n: c_int,x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular Bessel function Y_n(x), nmin <= n <= nmax
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_Yn_array(
        nmin: c_int,
        nmax: c_int,
        x: c_double,
        result_array: *mut c_double) -> c_int;

    /// Regular modified Bessel function I_0(x)
    /// exceptions: GSL_EOVRFLW
    fn gsl_sf_bessel_I0_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular modified Bessel function I_1(x)
    /// exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_I1_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular modified Bessel function I_n(x)
    /// exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_In_e(
        n: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular modified Bessel function  I_n(x) for n=nmin,...,nmax
    /// nmin >=0, nmax >= nmin
    /// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_In_array(
        nmin: c_int,
        nmax: c_int,
        x: c_double,
        result_array: *mut c_double) -> c_int;

    /// Scaled regular modified Bessel function
    ///  exp(-|x|) I_0(x)
    /// exceptions: none
    fn gsl_sf_bessel_I0_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Scaled regular modified Bessel function
    ///  exp(-|x|) I_1(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_I1_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Scaled regular modified Bessel function
    ///  exp(-|x|) I_n(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_In_scaled_e(
        n: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Scaled regular modified Bessel function
    ///  exp(-|x|) I_n(x)  for n=nmin,...,nmax
    /// nmin >=0, nmax >= nmin
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_In_scaled_array(
        nmin: c_int,
        nmax: c_int,
        x: c_double,
        result_array: *mut c_double) -> c_int;

    /// Irregular modified Bessel function K_0(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_K0_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular modified Bessel function K_1(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_K1_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular modified Bessel function K_n(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_Kn_e(
        n: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular modified Bessel function  K_n(x)  for n=nmin,...,nmax
    /// x > 0.0, nmin >=0, nmax >= nmin
    /// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
    fn gsl_sf_bessel_Kn_array(
        nmin: c_int,
        nmax: c_int,
        x: c_double,
        result_array: *mut c_double) -> c_int;

    /// Scaled irregular modified Bessel function
    ///  exp(x) K_0(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM
    fn gsl_sf_bessel_K0_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Scaled irregular modified Bessel function
    ///  exp(x) K_1(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_K1_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int; 

    /// Scaled irregular modified Bessel function
    ///  exp(x) K_n(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_Kn_scaled_e(
        n: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Scaled irregular modified Bessel function  exp(x) K_n(x)  for n=nmin,...,nmax
    /// x > 0.0, nmin >=0, nmax >= nmin
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_Kn_scaled_array(
        nmin: c_int,
        nmax: c_int,
        x: c_double,
        result_array: *mut c_double) -> c_int;

    /// Regular spherical Bessel function j_0(x) = sin(x)/x
    fn gsl_sf_bessel_j0_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular spherical Bessel function j_1(x) = (sin(x)/x - cos(x))/x
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_j1_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular spherical Bessel function j_2(x) = ((3/x^2 - 1)sin(x) - 3cos(x)/x)/x
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_j2_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular spherical Bessel function j_l(x)
    /// l >= 0, x >= 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_jl_e(
        l: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular spherical Bessel function j_l(x) for l=0,1,...,lmax
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_jl_array(
        lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;

    /// Regular spherical Bessel function j_l(x) for l=0,1,...,lmax
    /// Uses Steed's method.
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_jl_steed_array(
        lmax: c_int, x: c_double, jl_x_array: *mut c_double) -> c_int;

    /// Irregular spherical Bessel function y_0(x)
    fn gsl_sf_bessel_y0_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular spherical Bessel function y_1(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_y1_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular spherical Bessel function y_2(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_y2_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular spherical Bessel function y_l(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_yl_e(
        l: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular spherical Bessel function y_l(x) for l=0,1,...,lmax
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_yl_array(
        lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;

    /// Regular scaled modified spherical Bessel function
    /// Exp[-|x|] i_0(x)
    /// exceptions: none
    fn gsl_sf_bessel_i0_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular scaled modified spherical Bessel function
    /// Exp[-|x|] i_1(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_i1_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular scaled modified spherical Bessel function
    /// Exp[-|x|] i_2(x)
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_i2_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular scaled modified spherical Bessel functions
    /// Exp[-|x|] i_l(x)
    /// i_l(x) = Sqrt[Pi/(2x)] BesselI[l+1/2,x]
    /// l >= 0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_il_scaled_e(
        l: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular scaled modified spherical Bessel functions
    /// Exp[-|x|] i_l(x)
    /// for l=0,1,...,lmax
    /// exceptions: GSL_EUNDRFLW
    fn gsl_sf_bessel_il_scaled_array(
        lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;

    /// Irregular scaled modified spherical Bessel function
    /// Exp[x] k_0(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_k0_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular modified spherical Bessel function
    /// Exp[x] k_1(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
    fn gsl_sf_bessel_k1_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular modified spherical Bessel function
    /// Exp[x] k_2(x)
    /// x > 0.0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
    fn gsl_sf_bessel_k2_scaled_e(
        x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular modified spherical Bessel function
    /// Exp[x] k_l[x]
    /// k_l(x) = Sqrt[Pi/(2x)] BesselK[l+1/2,x]
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_kl_scaled_e(
        l: c_int, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular scaled modified spherical Bessel function
    /// Exp[x] k_l(x)
    /// for l=0,1,...,lmax
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_kl_scaled_array(
        lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;

    /// Regular cylindrical Bessel function J_nu(x)
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_Jnu_e(
        nu: c_double, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Irregular cylindrical Bessel function Y_nu(x)
    fn gsl_sf_bessel_Ynu_e(
        nu: c_double, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Regular cylindrical Bessel function J_nu(x)
    /// evaluated at a series of x values. The array
    /// contains the x values. They are assumed to be
    /// strictly ordered and positive. The array is
    /// over-written with the values of J_nu(x_i).
    /// exceptions: GSL_EDOM, GSL_EINVAL
    fn gsl_sf_bessel_sequence_Jnu_e(
        nu: c_double, mode: c_int, size: size_t, v: *mut c_double) -> c_int;

    /// Scaled modified cylindrical Bessel functions
    /// Exp[-|x|] BesselI[nu, x]
    /// x >= 0, nu >= 0
    /// exceptions: GSL_EDOM
    fn gsl_sf_bessel_Inu_scaled_e(
        nu: c_double, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Modified cylindrical Bessel functions
    /// BesselI[nu, x]
    /// x >= 0, nu >= 0
    /// exceptions: GSL_EDOM, GSL_EOVRFLW
    fn gsl_sf_bessel_Inu_e(
        nu: c_double, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Scaled modified cylindrical Bessel functions
    /// Exp[+|x|] BesselK[nu, x]
    /// x > 0, nu >= 0
    /// exceptions: GSL_EDOM
    fn gsl_sf_bessel_Knu_scaled_e10_e(
        nu: c_double,
        x: c_double,
        result: *mut gsl_sf::gsl_sf_result_e10) -> c_int;

    /// Modified cylindrical Bessel functions
    /// BesselK[nu, x]
    /// x > 0, nu >= 0
    /// exceptions: GSL_EDOM, GSL_EUNDRFLW
    fn gsl_sf_bessel_Knu_e(
        nu: c_double, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// Logarithm of modified cylindrical Bessel functions.
    /// Log[BesselK[nu, x]]
    /// x > 0, nu >= 0
    /// exceptions: GSL_EDOM
    fn gsl_sf_bessel_lnKnu_e(
        nu: c_double, x: c_double, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// s'th positive zero of the Bessel function J_0(x).
    fn gsl_sf_bessel_zero_J0_e(
        s: c_uint, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// s'th positive zero of the Bessel function J_1(x).
    fn gsl_sf_bessel_zero_J1_e(
        s: c_uint, result: *mut gsl_sf::gsl_sf_result) -> c_int;

    /// s'th positive zero of the Bessel function J_nu(x).
    fn gsl_sf_bessel_zero_Jnu_e(
        nu: c_double, s: c_uint, result: *mut gsl_sf::gsl_sf_result) -> c_int;
}

/// Regular Bessel Function J_0(x)
pub fn gslbessel_j0r(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_J0_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular Bessel Function J_1(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_j1r(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_J1_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular Bessel Function J_n(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_jnr(n: i32, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Jn_e(n, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular Bessel Function J_n(x),  nmin <= n <= nmax
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_jnr_array(nmin: i32, nmax: i32, x: f64) -> Array<f64>
{
    if nmin > nmax
    {
        panic!("Invalid nmin, nmax values: {} vs {}", nmin, nmax);
    }
    let size = (nmax - nmin) as usize;
    let mut result_array = Array::<f64>::new_filled(0f64, size, Order::Row);
    unsafe {
        if gsl_sf_bessel_Jn_array(nmin, nmax, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");            
        }
    }
    result_array
}

/// Irregular Bessel function Y_0(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_y0r(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Y0_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular Bessel function Y_1(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_y1r(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Y1_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular Bessel function Y_n(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_ynr(n: i32, x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Yn_e(n, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular Bessel function Y_n(x), nmin <= n <= nmax
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_ynr_array(nmin: i32, nmax: i32, x: f64) -> Array<f64>
{
    if nmin > nmax
    {
        panic!("Invalid nmin, nmax values: {} vs {}", nmin, nmax);
    }
    let size = (nmax - nmin) as usize;
    let mut result_arr = Array::<f64>::new_filled(0f64, size, Order::Row);
    unsafe {
        if gsl_sf_bessel_Yn_array(nmin, nmax, x, result_arr.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_arr
}

/// Regular modified Bessel function I_0(x)
/// exceptions: GSL_EOVRFLW
pub fn gslbessel_i0r(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_I0_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular modified Bessel function I_1(x)
/// exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_i1r(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_I1_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular modified Bessel function I_n(x)
/// exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_inr(n: i32, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_In_e(n, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular modified Bessel function  I_n(x) for n=nmin,...,nmax
/// nmin >=0, nmax >= nmin
/// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_inr_array(nmin: i32, nmax: i32, x: f64) -> Array<f64>
{
    if nmin > nmax
    {
        panic!("Invalid nmin, nmax values: {} vs {}", nmin, nmax);
    }
    let size = (nmax - nmin) as usize;
    let mut result_arr = Array::<f64>::new_filled(0f64, size, Order::Row);
    unsafe {
        if gsl_sf_bessel_In_array(nmin, nmax, x, result_arr.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_arr
}

/// Scaled regular modified Bessel function
///  exp(-|x|) I_0(x)
/// exceptions: none
pub fn gslbessel_i0r_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_I0_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Scaled regular modified Bessel function
///  exp(-|x|) I_1(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_i1r_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_I1_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Scaled regular modified Bessel function
///  exp(-|x|) I_n(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_inr_scaled(n: i32, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_In_scaled_e(n, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Scaled regular modified Bessel function
///  exp(-|x|) I_n(x)  for n=nmin,...,nmax
/// nmin >=0, nmax >= nmin
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_inr_scaled_array(nmin: i32, nmax: i32, x: f64) -> Array<f64>
{
    if nmin > nmax
    {
        panic!("Invalid nmin, nmax values: {} vs {}", nmin, nmax);
    }
    let size = (nmax - nmin) as usize;
    let mut result_array = Array::<f64>::new_filled(0f64, size, Order::Row);
    unsafe {
        if gsl_sf_bessel_In_scaled_array(nmin, nmax, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_array
}

/////////////////////////
// Irregular Functions //
/////////////////////////

/// Irregular modified Bessel function K_0(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_k0i(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_K0_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular modified Bessel function K_1(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_k1i(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_K1_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Irregular modified Bessel function K_n(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_kni(n: i32, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Kn_e(n, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular modified Bessel function  K_n(x)  for n=nmin,...,nmax
/// x > 0.0, nmin >=0, nmax >= nmin
/// exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
pub fn gslbessel_nn_array(nmin: i32, nmax: i32, x: f64) -> Array<f64>
{
    if nmin > nmax
    {
        panic!("Invalid nmin, nmax values: {} vs {}", nmin, nmax);
    }
    let size = (nmax - nmin) as usize;
    let mut result_array = Array::<f64>::new_filled(0f64, size, Order::Row);
    unsafe {
        if gsl_sf_bessel_Kn_array(nmin, nmax, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_array
}

/// Scaled irregular modified Bessel function
///  exp(x) K_0(x)
/// x > 0.0
/// exceptions: GSL_EDOM
pub fn gslbessel_k0i_scaled(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_K0_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Scaled irregular modified Bessel function
///  exp(x) K_1(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_k1i_scaled(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_K1_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Scaled irregular modified Bessel function
///  exp(x) K_n(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_kni_scaled(n: i32, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Kn_scaled_e(n, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Scaled irregular modified Bessel function  exp(x) K_n(x)  for n=nmin,...,nmax
/// x > 0.0, nmin >=0, nmax >= nmin
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_kni_scaled_array(nmin: i32, nmax: i32, x: f64) -> Array<f64>
{
    if nmin > nmax
    {
        panic!("Invalid nmin, nmax values: {} vs {}", nmin, nmax);
    }
    let size = (nmax - nmin) as usize;
    let mut result_array = Array::<f64>::new_filled(0f64, size, Order::Row);
    unsafe {
        if gsl_sf_bessel_Kn_scaled_array(nmin, nmax, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_array
}

/// Regular spherical Bessel function j_0(x) = sin(x)/x
pub fn gslbessel_j0i(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_j0_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular spherical Bessel function j_1(x) = (sin(x)/x - cos(x))/x
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_j1i(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_j1_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Regular spherical Bessel function j_2(x) = ((3/x^2 - 1)sin(x) - 3cos(x)/x)/x
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_j2i(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_j2_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular spherical Bessel function j_l(x)
/// l >= 0, x >= 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_jli(l: i32, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Kn_scaled_e(l, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Regular spherical Bessel function j_l(x) for l=0,1,...,lmax
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_jl_array(lmax: usize, x: f64) -> Array<f64>
{
    let mut result_array = Array::<f64>::new_filled(0f64, lmax, Order::Row);
    unsafe {
        if gsl_sf_bessel_jl_array(lmax as i32, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_array
}

/// Regular spherical Bessel function j_l(x) for l=0,1,...,lmax
/// Uses Steed's method.
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_jlr_steed_array(lmax: usize, x: f64) -> Array<f64>
{
    let mut result_array = Array::<f64>::new_filled(0f64, lmax, Order::Row);
    unsafe {
        if gsl_sf_bessel_jl_steed_array(lmax as i32, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_array
}

/// Irregular spherical Bessel function y_0(x)
pub fn gslbessel_y0i(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_y0_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Irregular spherical Bessel function y_1(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_y1i(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_y1_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)        
}

/// Irregular spherical Bessel function y_2(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_y2i(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_y2_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular spherical Bessel function y_l(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_yli(l: i32, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_yl_e(l, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular spherical Bessel function y_l(x) for l=0,1,...,lmax
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_yli_array(lmax: usize, x: f64) -> Array<f64>
{
    let mut result_array = Array::<f64>::new_filled(0f64, lmax, Order::Row);
    unsafe {
        if gsl_sf_bessel_jl_steed_array(lmax as i32, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_array
}

/// Regular scaled modified spherical Bessel function
/// Exp[-|x|] i_0(x)
/// exceptions: none
pub fn gslbessel_i0_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_i0_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular scaled modified spherical Bessel function
/// Exp[-|x|] i_1(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_i1_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_i1_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Regular scaled modified spherical Bessel function
/// Exp[-|x|] i_2(x)
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_i2_scaled(x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_i2_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Regular scaled modified spherical Bessel functions
/// Exp[-|x|] i_l(x)
/// i_l(x) = Sqrt[Pi/(2x)] BesselI[l+1/2,x]
/// l >= 0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_il_scaled(l: i32, x: f64) -> (f64, f64)
{
    if l < 0
    {
        panic!("l must be positive");
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_il_scaled_e(l, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular scaled modified spherical Bessel functions
/// Exp[-|x|] i_l(x)
/// for l=0,1,...,lmax
/// exceptions: GSL_EUNDRFLW
pub fn gslbessel_il_scaled_array(lmax: usize, x: f64) -> Array<f64>
{
    let mut result_array = Array::<f64>::new_filled(0f64, lmax, Order::Row);
    unsafe {
        if gsl_sf_bessel_il_scaled_array(lmax as i32, x, result_array.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_array
}

/// Irregular scaled modified spherical Bessel function
/// Exp[x] k_0(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_k0_scaled(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_k0_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular modified spherical Bessel function
/// Exp[x] k_1(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
pub fn gslbessel_k1_scaled(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_k1_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular modified spherical Bessel function
/// Exp[x] k_2(x)
/// x > 0.0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
pub fn gslbessel_k2_scaled(x: f64) -> (f64, f64)
{
    if x <= 0.0
    {
        panic!("Invalid value for x: {}", x);
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_k2_scaled_e(x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular modified spherical Bessel function
/// Exp[x] k_l[x]
/// k_l(x) = Sqrt[Pi/(2x)] BesselK[l+1/2,x]
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_kl_scaled(l: i32, x: f64) -> (f64, f64)
{
    if l < 0
    {
        panic!("l must be positive");
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_kl_scaled_e(l, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular scaled modified spherical Bessel function
/// Exp[x] k_l(x)
/// for l=0,1,...,lmax
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_kl_scaled_array(lmax: usize, x: f64) -> Array<f64>
{
    let mut result_arr = Array::<f64>::new_filled(0f64, lmax, Order::Row);
    unsafe {
        if gsl_sf_bessel_kl_scaled_array(lmax as c_int, x, result_arr.as_mut_ptr()) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    result_arr
}


/// Regular cylindrical Bessel function J_nu(x)
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_jnur(nu: f64, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Jnu_e(nu, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Irregular cylindrical Bessel function Y_nu(x)
pub fn gslbessel_ynui(nu: f64, x: f64) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Ynu_e(nu, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Regular cylindrical Bessel function J_nu(x)
/// evaluated at a series of x values. The array
/// contains the x values. They are assumed to be
/// strictly ordered and positive. The array is
/// over-written with the values of J_nu(x_i).
/// exceptions: GSL_EDOM, GSL_EINVAL
pub fn gslbessel_sequence_jnur(nu: f64, v: &Array<f64>) -> Array<f64>
{
    let mut _v = v.clone();

    unsafe {
        if 0 != gsl_sf_bessel_sequence_Jnu_e(
            nu, gsl_sf::GSL_PREC_DOUBLE, v.len(), _v.as_mut_ptr())
        {
            panic!("Bessel calulation failed");
        }
    }

    _v
}

/// Scaled modified cylindrical Bessel functions
/// Exp[-|x|] BesselI[nu, x]
/// x >= 0, nu >= 0
/// exceptions: GSL_EDOM
pub fn gslbessel_inu_scaled(nu: f64, x: f64) -> (f64, f64)
{
    if x < 0f64 || nu < 0f64
    {
        panic!("nu, x out of range");
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Inu_scaled_e(nu, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Modified cylindrical Bessel functions
/// BesselI[nu, x]
/// x >= 0, nu >= 0
/// exceptions: GSL_EDOM, GSL_EOVRFLW
pub fn gslbessel_inu(nu: f64, x: f64) -> (f64, f64)
{
    if x < 0f64 || nu < 0f64
    {
        panic!("nu, x out of range");
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Inu_scaled_e(nu, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Scaled modified cylindrical Bessel functions
/// Exp[+|x|] BesselK[nu, x]
/// x > 0, nu >= 0
/// exceptions: GSL_EDOM
pub fn gslbessel_knu_scaled(nu: f64, x: f64) -> (f64, f64)
{
    if x < 0f64 || nu < 0f64
    {
        panic!("nu, x out of range");
    }
    let mut s = gsl_sf::gsl_sf_result_e10_struct{val: 0f64, err: 0f64, e10: 0};
    unsafe {
        if gsl_sf_bessel_Knu_scaled_e10_e(nu, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// Modified cylindrical Bessel functions
/// BesselK[nu, x]
/// x > 0, nu >= 0
/// exceptions: GSL_EDOM, GSL_EUNDRFLW
pub fn gslbessel_knu(nu: f64, x: f64) -> (f64, f64)
{
    if x <= 0f64 || nu < 0f64
    {
        panic!("nu, x out of range");
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_Knu_e(nu, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// Logarithm of modified cylindrical Bessel functions.
/// Log[BesselK[nu, x]]
/// x > 0, nu >= 0
/// exceptions: GSL_EDOM
pub fn gslbessel_lnknu(nu: f64, x: f64) -> (f64, f64)
{
    if x <= 0f64 || nu < 0f64
    {
        panic!("nu, x out of range");
    }
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_lnKnu_e(nu, x, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// n'th positive zero of the Bessel function J_0(x).
pub fn gslbessel_zero_j0(n: u32) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_zero_J0_e(n, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}

/// n'th positive zero of the Bessel function J_1(x).
pub fn gslbessel_zero_j1(n: u32) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_zero_J1_e(n, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)    
}

/// n'th positive zero of the Bessel function J_nu(x).
pub fn gslbessel_zero_jnu(nu: f64, n: u32) -> (f64, f64)
{
    let mut s = gsl_sf::gsl_sf_result_struct{val: 0f64, err: 0f64};
    unsafe {
        if gsl_sf_bessel_zero_Jnu_e(nu, n, &mut s) != 0
        {
            panic!("Bessel calculation failed");
        }
    }
    (s.val, s.err)
}
