/// Polynomial function utilities.
use array::Array;
use array::Order;
use std::f64;
use libc::{c_int, c_double, size_t};

#[link(name = "gsl")]
extern
{
    /// The functions described here evaluate the polynomial:
    ///   P(x) = c[0] + c[1]*x + c[2]*x^2 + ... + c[len-1]*x^{len-1}
    /// using Horner’s method for stability.
    
    /// Evaluate a polynomial with real coefficients for the real variable x
    fn gsl_poly_eval(c: *const c_double, len: c_int, x: c_double) -> c_double;
    
    /// Evaluate a polynomial and its derivatives storing the results in the
    /// array res of size lenres. The output array contains the values of:
    ///   d^k P/d x^k
    /// for the specified value of x starting with k = 0.
    fn gsl_poly_eval_derivs(
        c: *const c_double,
        lenc: size_t,
        x: c_double,
        res: *mut c_double,
        lenres: size_t) -> c_int;

    /// The functions described here manipulate polynomials stored in Newton’s
    /// divided-difference representation. The use of divided-differences is
    /// described in Abramowitz & Stegun sections 25.1.4 and 25.2.26, and Burden
    /// and Faires, chapter 3, and discussed briefly below.
    /// Given a function f(x), an nth degree interpolating polynomial P_{n}(x)
    /// can be constructed which agrees with f at n+1 distinct points
    /// x_0,x_1,...,x_{n}. This polynomial can be written in a form known as
    /// Newton’s divided-difference representation:
    ///   P_n(x) = f(x_0) + \sum_(k=1)^n [x_0,x_1,...,x_k] (x-x_0)(x-x_1)...(x-x_(k-1))
    /// where the divided differences [x_0,x_1,...,x_k] are defined in section
    /// 25.1.4 of Abramowitz and Stegun. Additionally, it is possible to
    /// construct an interpolating polynomial of degree 2n+1 which also matches
    /// the first derivatives of f at the points x_0,x_1,...,x_n. This is called
    /// the Hermite interpolating polynomial and is defined as
    ///   H_(2n+1)(x) = f(z_0) + sum(k=1)^(2n+1) [z_0,z_1,...,z_k] (x-z_0)(x-z_1)...(x-z_(k-1))
    /// where the elements of z = {x_0,x_0,x_1,x_1,...,x_n,x_n} are defined by
    /// z_{2k} = z_{2k+1} = x_k. The divided-differences [z_0,z_1,...,z_k] are
    /// discussed in Burden and Faires, section 3.4.

    /// This function computes a divided-difference representation of the
    /// interpolating polynomial for the points (x, y) stored in the arrays xa
    /// and ya of length size. On output the divided-differences of (xa,ya) are
    /// stored in the array dd, also of length size. Using the notation above,
    /// dd[k] = [x_0,x_1,...,x_k].
    fn gsl_poly_dd_init(
        dd: *mut c_double,
        xa: *const c_double,
        ya: *const c_double,
        size: size_t) -> c_int;

    /// Evaluate the polynomial stored in divided-difference form in the arrays
    /// dd and xa of length size at the point x.
    fn gsl_poly_dd_eval(
        dd: *const c_double,
        xa: *const c_double,
        size: size_t,
        x: c_double) -> c_double;

    /// This function converts the divided-difference representation of a
    /// polynomial to a Taylor expansion. The divided-difference representation
    /// is supplied in the arrays dd and xa of length size. On output the Taylor
    /// coefficients of the polynomial expanded about the point xp are stored in
    /// the array c also of length size. A workspace of length size must be
    /// provided in the array w.
    fn gsl_poly_dd_taylor(
        c: *mut c_double,
        xp: c_double,
        dd: *const c_double,
        xa: *const c_double,
        size: size_t,
        w: *const c_double) -> c_int;

    /// This function computes a divided-difference representation of the
    /// interpolating Hermite polynomial for the points (x, y) stored in the
    /// arrays xa and ya of length size. Hermite interpolation constructs
    /// polynomials which also match first derivatives dy/dx which are provided
    /// in the array dya also of length size. The first derivatives can be
    /// incorported into the usual divided-difference algorithm by forming a new
    /// dataset z = \{x_0,x_0,x_1,x_1,...\}, which is stored in the array za of
    /// length 2*size on output. On output the divided-differences of the
    /// Hermite representation are stored in the array dd, also of length
    /// 2*size. Using the notation above, dd[k] = [z_0,z_1,...,z_k]. The
    /// resulting Hermite polynomial can be evaluated by calling
    /// gsl_poly_dd_eval and using za for the input argument xa.
    fn gsl_poly_dd_hermite_init(
        dd: *mut c_double,
        za: *mut c_double,
        xa: *const c_double,
        ya: *const c_double,
        dya: *const c_double,
        size: size_t) -> c_int;

    /// This function finds the real roots of the quadratic equation,
    ///   a x^2 + b x + c = 0
    /// The number of real roots (either zero, one or two) is returned, and
    /// their locations are stored in x0 and x1. If no real roots are found then
    /// x0 and x1 are not modified. If one real root is found (i.e. if a=0) then
    /// it is stored in x0. When two real roots are found they are stored in x0
    /// and x1 in ascending order. The case of coincident roots is not
    /// considered special. For example (x-1)^2=0 will have two roots, which
    /// happen to have exactly equal values.  The number of roots found depends
    /// on the sign of the discriminant b^2 - 4 a c. This will be subject to
    /// rounding and cancellation errors when computed in double precision, and
    /// will also be subject to errors if the coefficients of the polynomial are
    /// inexact. These errors may cause a discrete change in the number of
    /// roots. However, for polynomials with small integer coefficients the
    /// discriminant can always be computed exactly.
    fn gsl_poly_solve_quadratic(
        a: c_double,
        b: c_double,
        c: c_double,
        x0: *mut c_double,
        x1: *mut c_double) -> c_int;

    /// Find the real roots of the cubic equation,
    ///   x^3 + a x^2 + b x + c = 0
    /// with a leading coefficient of unity. The number of real roots (either
    /// one or three) is returned, and their locations are stored in x0, x1 and
    /// x2. If one real root is found then only x0 is modified. When three real
    /// roots are found they are stored in x0, x1 and x2 in ascending order. The
    /// case of coincident roots is not considered special. For example, the
    /// equation (x-1)^3=0 will have three roots with exactly equal values. As
    /// in the quadratic case, finite precision may cause equal or
    /// closely-spaced real roots to move off the real axis into the complex
    /// plane, leading to a discrete change in the number of real roots.
    fn gsl_poly_solve_cubic(
        a: c_double,
        b: c_double,
        c: c_double,
        x0: *mut c_double,
        x1: *mut c_double,
        x2: *mut c_double) -> c_int;

}

/// Evaluate a polynomial with real coefficients for the real variable x
pub fn poly_eval(c: &Array<f64>, x: f64) -> f64
{
    unsafe {
        gsl_poly_eval(c.as_ptr(), c.len() as c_int, x as c_double)
    }
}

/// Evaluate a polynomial and its k derivatives, placing it in an array
pub fn poly_eval_derivs(c: &Array<f64>, x: f64, k: usize) -> Array<f64>
{
    if k == 0
    {
        panic!("Cannot take the 0th derivative");
    }
    let mut derivs = Array::<f64>::new_filled(0f64, k, Order::Row);

    unsafe {
        let result = gsl_poly_eval_derivs(
            c.as_ptr(),
            c.len() as size_t,
            x as c_double,
            derivs.as_mut_ptr(),
            k as size_t);
        if result != 0
        {
            panic!("Failed to calculate the derivative at x: {}", x);
        }
    }

    return derivs;
}

/// Computes a divided-difference representation of the interpolating polynomial
/// for the points (x, y) in the form dd[k] = [x_0,x_1,...,x_k]
pub fn poly_divdiff_init(xa: &Array<f64>, ya: &Array<f64>) -> Array<f64>
{
    if xa.len() != ya.len()
    {
        panic!("Inconsistent sizes for xa, ya: {} != {}", xa.len(), ya.len());
    }
    let mut dd = Array::<f64>::new_filled(0f64, xa.len(), Order::Row);
    unsafe {
        let result = gsl_poly_dd_init(
            dd.as_mut_ptr(),
            xa.as_ptr(),
            ya.as_ptr(),
            xa.len() as size_t);
        if result != 0
        {
            panic!("An error occured while ");
        }
    }
    return dd;
}

// Evaluates the divided-difference polynomial at x
pub fn poly_divdiff_eval(dd: &Array<f64>, xa: &Array<f64>, x: f64) -> f64
{
    if dd.len() != xa.len()
    {
        panic!("Inconsistent sizes for xa, ya, dd");
    }
    unsafe{
        gsl_poly_dd_eval(
            dd.as_ptr(),
            xa.as_ptr(),
            dd.len() as size_t,
            x as c_double) as f64
    }
}

/// Converts the divided-difference polynomial to a Taylor expansion around a point
pub fn poly_divdiff_to_taylor(xp: f64, dd: &Array<f64>, xa: &Array<f64>) -> Array<f64>
{
    let mut c = Array::<f64>::new_filled(0f64, xa.len(), Order::Row);
    unsafe {
        let result = gsl_poly_dd_taylor(
            c.as_mut_ptr(),
            xp,
            dd.as_ptr(),
            xa.as_ptr(),
            xa.len() as size_t,
            Array::<f64>::new_filled(0f64, xa.len(), Order::Row).as_mut_ptr());
        if result != 0
        {
            panic!("Talyor expansion conversion failed");
        }
    }
    return c;
}

/// Computes a divided-difference representation of the interpolating Hermite
/// polynomial for the points (x, y).
/// The resulting Hermite polynomial can be evaluated by calling
/// gsl_poly_dd_eval and using za for the input argument xa.
pub fn poly_divdiff_to_hermite(
    xa: &Array<f64>,
    ya: &Array<f64>,
    dya: &Array<f64>) -> (Array<f64>, Array<f64>)
{
    if xa.len() != ya.len() || dya.len() != xa.len()
    {
        panic!("Inconsistent sizes for xa, ya, dya");
    }
    let mut dd = Array::<f64>::new_filled(0f64, xa.len(), Order::Row);
    let mut za = Array::<f64>::new_filled(0f64, xa.len(), Order::Row);

    unsafe {
        let result = gsl_poly_dd_hermite_init(
            dd.as_mut_ptr(),
            za.as_mut_ptr(),
            xa.as_ptr(),
            ya.as_ptr(),
            dya.as_ptr(),
            xa.len());
        if result != 0
        {
            panic!("Computing representation for Hermite polynomial failed.");
        }
    }

    return (dd, za);
}

//////////////////////
// Quadratic Solver //
//////////////////////

// Finds the real roots of the quadratic equation a*x^2 + b*x + c = 0
pub fn poly_solve_quadratic(coeffs: [f64; 3]) -> Option<Vec<f64>>
{
    let mut x0 = f64::NAN;
    let mut x1 = f64::NAN;
    unsafe {
        let num_roots = gsl_poly_solve_quadratic(
            coeffs[0],
            coeffs[1],
            coeffs[2],
            &mut x0 as *mut c_double,
            &mut x1 as *mut c_double);
        if num_roots == 0
        {
            return None;
        } else if num_roots == 1 {
            return Some(vec![x0]);
        } else {
            return Some(vec![x0, x1]);
        }
    }
}

/// Find the real roots of the cubic equation, x^3 + a*x^2 + b*x + c = 0
pub fn poly_solve_cubic(coeffs: [f64; 3]) -> Vec<f64>
{
    let mut x0 = f64::NAN;
    let mut x1 = f64::NAN;
    let mut x2 = f64::NAN;
    unsafe {
        let num_roots = gsl_poly_solve_cubic(
            coeffs[0],
            coeffs[1],
            coeffs[2],
            &mut x0 as *mut c_double,
            &mut x1 as *mut c_double,
            &mut x2 as *mut c_double
        );
        if num_roots == 1 {
            return vec![x0];
        } else {
            return vec![x0, x1, x2];
        }
    }
}
