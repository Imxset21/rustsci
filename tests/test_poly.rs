/// Polynomial evaluation testing
#[macro_use]
extern crate rustsci;

use rustsci::array;
use rustsci::gsl_math;
use rustsci::gsl_poly;

const EPS: f64 = 2.2204460492503131e-14;

// Polynomial evaluation

#[test]
fn test_3_poly()
{
    let c = arr![1.0, 0.5, 0.3];
    let x: f64 = 0.5;
    let y: f64 = gsl_poly::poly_eval(&c, x);
    assert_epeq!(y, 1f64 + 0.5f64 * x + 0.3f64 * x * x, EPS);
}

#[test]
pub fn test_11_poly()
{
    let d = arr![1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1.];
    let x: f64 = 1.0;
    let y: f64 = gsl_poly::poly_eval(&d, x);
    assert_epeq!(y, 1.0f64, EPS);
}

// Quadratic
#[test]
pub fn test_quadratic_solve_1()
{
    match gsl_poly::poly_solve_quadratic([4.0, -20.0, 26.0]) {
        None => (),
        Some(_) => panic!("Quadratic test failed: invalid number of roots")
    };
}

#[test]
pub fn test_quadratic_solve_2()
{
    let (x0, x1) = match gsl_poly::poly_solve_quadratic([4.0, -20.0, 25.0]) {
        None => panic!("Quadratic test failed: invalid number of roots"),
        Some(x) => (x[0], x[1])
    };

    assert_epeq!(x0, 2.5, 1e-9);
    assert_epeq!(x1, 2.5, 1e-9);
    assert_eq!(x0, x1);
}

#[test]
pub fn test_quadratic_solve_3()
{
    let (x0, x1) = match gsl_poly::poly_solve_quadratic([4.0, -20.0, 21.0]) {
        None => panic!("Quadratic test failed: invalid number of roots"),
        Some(x) => (x[0], x[1])
    };

    assert_epeq!(x0, 1.5, 1e-9);
    assert_epeq!(x1, 3.5, 1e-9);
}

#[test]
pub fn test_quadratic_solve_4()
{
    let (x0, x1) = match gsl_poly::poly_solve_quadratic([4.0, 7.0, 0.0]) {
        None => panic!("Quadratic test failed: invalid number of roots"),
        Some(x) => (x[0], x[1])
    };

    assert_epeq!(x0, -1.75, 1e-9);
    assert_epeq!(x1, 0.0, 1e-9);
}

#[test]
pub fn test_quadratic_solve_5()
{
    let (x0, x1) = match gsl_poly::poly_solve_quadratic([5.0, 0.0, -20.0]) {
        None => panic!("Quadratic test failed: invalid number of roots"),
        Some(x) => (x[0], x[1])
    };
    
    assert_epeq!(x0, -2.0, 1e-9);
    assert_epeq!(x1, 2.0, 1e-9);
}

#[test]
pub fn test_quadratic_solve_6()
{
    let x0 = match gsl_poly::poly_solve_quadratic([0.0, 3.0, -21.0]) {
        None => panic!("Quadratic test failed: invalid number of roots"),
        Some(x) => x[0]
    };
    
    assert_epeq!(x0, 7.0, 1e-9);
}

#[test]
pub fn test_quadratic_solve_7()
{
    match gsl_poly::poly_solve_quadratic([0.0, 0.0, 1.0]) {
        None => (),
        Some(_) => panic!("Quadratic test failed: invalid number of roots")
    };
}


///////////
// Cubic //
///////////

#[test]
pub fn test_cubic_solve_1()
{
    let roots = gsl_poly::poly_solve_cubic([0.0, 0.0, -27.0]);

    assert_eq!(roots.len(), 1);
    assert_epeq!(roots[0], 3.0, 1e-9);
}

#[test]
pub fn test_cubic_solve_2()
{
    let roots = gsl_poly::poly_solve_cubic([-51.0, 867.0, -4913.0]);

    assert_eq!(roots.len(), 3);
    assert_epeq!(roots[0], 17.0, 1e-9);
    assert_epeq!(roots[1], 17.0, 1e-9);
    assert_epeq!(roots[2], 17.0, 1e-9);
}

#[test]
pub fn test_cubic_solve_3()
{
    let roots = gsl_poly::poly_solve_cubic([-57.0, 1071.0, -6647.0]);

    assert_eq!(roots.len(), 3);
    assert_epeq!(roots[0], 17.0, 1e-9);
    assert_epeq!(roots[1], 17.0, 1e-9);
    assert_epeq!(roots[2], 23.0, 1e-9);
}

#[test]
pub fn test_cubic_solve_4()
{
    let roots = gsl_poly::poly_solve_cubic([-11.0, -493.0, 6647.0]);

    assert_eq!(roots.len(), 3);
    assert_epeq!(roots[0], -23.0, 1e-9);
    assert_epeq!(roots[1], 17.0, 1e-9);
    assert_epeq!(roots[2], 17.0, 1e-9);
}

#[test]
pub fn test_cubic_solve_5()
{
    let roots = gsl_poly::poly_solve_cubic([-143.0, 5087.0, -50065.0]);

    assert_eq!(roots.len(), 3);
    assert_epeq!(roots[0], 17.0, 1e-9);
    assert_epeq!(roots[1], 31.0, 1e-9);
    assert_epeq!(roots[2], 95.0, 1e-9);
}

#[test]
pub fn test_cubic_solve_6()
{
    let roots = gsl_poly::poly_solve_cubic([-109.0, 803.0, 50065.0]);

    assert_eq!(roots.len(), 3);
    assert_epeq!(roots[0], -17.0, 1e-9);
    assert_epeq!(roots[1], 31.0, 1e-9);
    assert_epeq!(roots[2], 95.0, 1e-9);
}

// DD Polynomials

#[test]
pub fn test_dd_eval_1()
{
    let xa = arr![0.16, 0.97, 1.94, 2.74, 3.58, 3.73, 4.70];
    let ya = arr![0.73, 1.11, 1.49, 1.84, 2.30, 2.41, 3.07];

    let dd_expected: array::Array<f64> = arr![
        7.30000000000000e-01,
        4.69135802469136e-01,
        -4.34737219941284e-02,
        2.68681098870099e-02,
        -3.22937056934996e-03,
        6.12763259971375e-03,
        -6.45402453527083e-03];

    let dd = gsl_poly::poly_divdiff_init(&xa, &ya);

    for (dd_v, dd_e) in (&dd).into_iter().zip((&dd_expected).into_iter())
    {
        assert_epeq!(*dd_v, *dd_e, 1e-10);
    }
    
    for (xa_i, ya_i) in (&xa).into_iter().zip((&ya).into_iter())
    {
        let y = gsl_poly::poly_divdiff_eval(&dd, &xa, *xa_i);
        assert_epeq!(y, *ya_i, 1e-10);
    }

    let coeff = gsl_poly::poly_divdiff_to_taylor(1.5f64, &dd, &xa);
    
    for (xa_i, ya_i) in (&xa).into_iter().zip((&ya).into_iter())
    {
        let y = gsl_poly::poly_eval(&coeff, *xa_i - 1.5f64);
        assert_epeq!(y, *ya_i, 1e-10);
    }
}

#[test]
#[ignore]
pub fn test_dd_eval_2()
{
    // TODO: Find out why this isn't working
    let xa: array::Array<f64> = arr![1.3, 1.6, 1.9];
    let ya: array::Array<f64> = arr![0.6200860, 0.4554022, 0.2818186];
    let dya: array::Array<f64> = arr![-0.5220232, -0.5698959, -0.5811571];

    let dd_expected: array::Array<f64> =
        arr![6.200860000000e-01,
             -5.220232000000e-01,
             -8.974266666667e-02,
             6.636555555556e-02,
             2.666666666662e-03,
             -2.774691357989e-03];

    let (dd, za) = gsl_poly::poly_divdiff_to_hermite(&xa, &ya, &dya);
    println!("{:?}", za);
    for (dd_v, dd_e) in (&dd).into_iter().zip((&dd_expected).into_iter())
    {
        assert_epeq!(*dd_v, *dd_e, 1e-10);
    }

    for (xa_i, ya_i) in (&xa).into_iter().zip((&ya).into_iter())
    {
        let y = gsl_poly::poly_divdiff_eval(&dd, &za, *xa_i);
        assert_epeq!(y, *ya_i, 1e-10);
    }

    for (xa_i, dya_i) in (&xa).into_iter().zip((&dya).into_iter())
    {
        let coeff = gsl_poly::poly_divdiff_to_taylor(*xa_i, &dd, &za);
        assert_epeq!(coeff[1], *dya_i, 1e-10);
    }
}

#[test]
pub fn test_poly_eval_derivs()
{
    let c = arr![1.0, -2.0, 3.0, -4.0, 5.0, -6.0];
    let x = -0.5f64;
    let dc = gsl_poly::poly_eval_derivs(&c, x, 6);

    assert_epeq!(dc[0], c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x + c[4]*x*x*x*x + c[5]*x*x*x*x*x, EPS);
    assert_epeq!(dc[1], c[1] + 2.0*c[2]*x + 3.0*c[3]*x*x + 4.0*c[4]*x*x*x + 5.0*c[5]*x*x*x*x, EPS);
    assert_epeq!(dc[2], 2.0*c[2] + 3.0*2.0*c[3]*x + 4.0*3.0*c[4]*x*x + 5.0*4.0*c[5]*x*x*x , EPS);
    assert_epeq!(dc[3], 3.0*2.0*c[3] + 4.0*3.0*2.0*c[4]*x + 5.0*4.0*3.0*c[5]*x*x, EPS);
    assert_epeq!(dc[4], 4.0*3.0*2.0*c[4] + 5.0*4.0*3.0*2.0*c[5]*x, EPS);
    assert_epeq!(dc[5], 5.0*4.0*3.0*2.0*c[5], EPS);
}
