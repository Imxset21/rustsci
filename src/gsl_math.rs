/// Basic GSL mathematical function alternatives for Rust's standard library
use libc::{c_int, c_double, c_uint};

#[link(name = "gsl")]
extern
{
    /// This function computes the value of log(1 + x) in a way that is accurate for
    /// small x.  It provides an alternative to the BSD math function log1p(x).
    fn gsl_log1p(x: c_double) -> c_double;

    /// This function computes the value of exp(x) − 1 in a way that is accurate
    /// for small x. It provides an alternative to the BSD math function expm1(x).
    fn gsl_expm1(x: c_double) -> c_double;

    /// This function computes the value of sqrt{x^2 + y^2} in a way that
    /// avoids overflow. It provides an alternative to the BSD math function
    /// hypot(x,y).
    fn gsl_hypot(x: c_double, y: c_double) -> c_double;

    /// This function computes the value of sqrt{x^2 + y^2 + z^2} in a way that
    /// avoids overflow.
    fn gsl_hypot3(x: c_double, y: c_double, z: c_double) -> c_double;

    /// This function computes the value of arccosh(x). It provides an
    /// alternative to the standard math function acosh(x).
    fn gsl_acosh(x: c_double) -> c_double;

    /// This function computes the value of arcsinh(x). It provides an
    /// alternative to the standard math function asinh(x).
    fn gsl_asinh(x: c_double) -> c_double;

    /// This function computes the value of arctanh(x). It provides an
    /// alternative to the standard math function atanh(x).
    fn gsl_atanh(x: c_double) -> c_double;

    /// This function computes the value of x * 2^e. It provides an alternative
    /// to the standard math function ldexp(x,e).
    fn gsl_ldexp(x: c_double) -> c_double;

    /// This function splits the number x into its normalized fraction f and
    /// exponent e, such that x = f * 2^e and 0.5 <= f < 1. The function returns
    /// f and stores the exponent in e. If x is zero, both f and e are set to
    /// zero. This function provides an alternative to the standard math
    /// function frexp(x, e).
    fn gsl_frexp (x: c_double, e: *mut c_int) -> c_double;

    /// Efficient power function
    fn gsl_pow_int(x: c_double, n: c_int) -> c_double;
    fn gsl_pow_uint(x: c_double, n: c_uint) -> c_double;
    fn gsl_pow_2(x: c_double) -> c_double;
    fn gsl_pow_3(x: c_double) -> c_double;
    fn gsl_pow_4(x: c_double) -> c_double;
    fn gsl_pow_5(x: c_double) -> c_double;
    fn gsl_pow_6(x: c_double) -> c_double;
    fn gsl_pow_7(x: c_double) -> c_double;
    fn gsl_pow_8(x: c_double) -> c_double;
    fn gsl_pow_9(x: c_double) -> c_double;

    /// This function determines whether x and y are approximately equal to a
    /// relative accuracy epsilon.  The relative accuracy is measured using an
    /// interval of size 2*delta, where delta = 2^k * epsilon and k is the
    /// maximum base-2 exponent of x and y as computed by the function frexp.
    /// If x and y lie within this interval, they are considered approximately
    /// equal and the function returns 0. Otherwise if x < y, the function
    /// returns -1, or if x > y, the function returns +1.  Note that x and y are
    /// compared to relative accuracy, so this function is not suitable for
    /// testing whether a value is approximately zero.  The implementation is
    /// based on the package fcmp by T.C. Belding.
    fn gsl_fcmp(x: c_double, y: c_double, epsilon: c_double) -> c_int;
}

/// Computes the value of log(1 + x) in a way that is accurate for small x.
pub fn gslmath_log1p(x: f64) -> f64
{
    unsafe {
        gsl_log1p(x)
    }
}

/// Computes the value of exp(x) − 1 in a way that is accurate for small x.
pub fn gslmath_expm1(x: f64) -> f64
{
    unsafe {
        gsl_expm1(x)
    }
}

/// Computes the value of sqrt{x^2 + y^2} in a way that avoids overflow
pub fn gslmath_hypot(x: f64, y: f64) -> f64
{
    unsafe {
        gsl_hypot(x, y)
    }
}

/// Computes the value of sqrt{x^2 + y^2 + z^2} in a way that avoids overflow.
pub fn gslmath_hypot3(x: f64, y: f64, z: f64) -> f64
{
    unsafe {
        gsl_hypot3(x, y, z)
    }
}

/// Computes the value of arccosh(x).
pub fn gslmath_acosh(x: f64) -> f64
{
    unsafe {
        gsl_acosh(x)
    }
}

/// Computes the value of arcsinh(x)
pub fn gslmath_asinh(x: f64) -> f64
{
    unsafe {
        gsl_asinh(x)
    }
}

/// Computes the value of arctanh(x).
pub fn gslmath_atanh(x: f64) -> f64
{
    unsafe {
        gsl_atanh(x)
    }
}

/// Computes the value of x * 2^e.
pub fn gslmath_ldexp(x: f64) -> f64
{
    unsafe {
        gsl_ldexp(x)
    }
}

/// Splits the number x into its normalized fraction f and exponent e, such that
/// x = f * 2^e and 0.5 <= f < 1. The function returns f and the exponent e. If
/// x is zero, both f and e are set to zero.
pub fn gslmath_frexp(x: f64) -> (f64, i32)
{
    unsafe {
        let mut _e: i32 = 0;
        let f = gsl_frexp(x, &mut _e as *mut c_int) as f64;
        (f, _e)
    }
}

/// Pattern-matched, efficient power function for positive powers
pub fn gslmath_pow(x: f64, n: u32) -> f64
{
    unsafe {
        match n {
            1u32 => x,
            2u32 => gsl_pow_2(x),
            3u32 => gsl_pow_3(x),
            4u32 => gsl_pow_4(x),
            5u32 => gsl_pow_5(x),
            6u32 => gsl_pow_6(x),
            7u32 => gsl_pow_7(x),
            8u32 => gsl_pow_8(x),
            9u32 => gsl_pow_9(x),
            _ => gsl_pow_uint(x, n),
        }
    }
}

/// Efficient power function for positive and negative powers
pub fn gslmath_spow(x: f64, n: i32) -> f64
{
    unsafe {
        gsl_pow_int(x, n)
    }
}

/// Determines whether x and y are approximately equal to a relative accuracy epsilon
pub fn gslmath_fcmp(x: f64, y: f64, epsilon: f64) -> bool
{
    unsafe {
        gsl_fcmp(x, y, epsilon) == 0
    }
}

/// Convenience macro for doing assert_eq using fcmp
#[macro_export]
macro_rules! assert_epeq {
    ($left:expr , $right:expr, $epsilon:expr) => ({
        match (&($left), &($right), &($epsilon)) {
            (left_val, right_val, epsilon) => {
                if !gsl_math::gslmath_fcmp(*left_val, *right_val, *epsilon) {
                    panic!("assertion failed: `(left == right)` \
                           (left: `{:?}`, right: `{:?}`)", left_val, right_val)
                }
            }
        }
    })
}
