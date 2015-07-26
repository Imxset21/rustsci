// extern crate libc; // <- Uncomment for static analysis
// mod array;         // <- Uncomment for static analysis

use array::Array;
use libc::{c_int, c_double, c_float};

#[link(name = "blas")]
extern
{
    /// Computes a vector-vector dot product for f64 ("double") precision.
    fn cblas_ddot(N: c_int,
                  X: *const c_double, incX: c_int,
                  Y: *const c_double, incY: c_int) -> c_double;

    /// Computes a vector-vector dot product for f32 ("single") precision.
    fn cblas_sdot(N: c_int,
                  X: *const c_float, incX: c_int,
                  Y: *const c_float, incY: c_int) -> c_float;

    /// Computes the sum of magnitudes of the vector elements, single precision
    fn cblas_sasum(N: c_int, X: *const c_float, incX: c_int) -> c_float;
    
    /// Computes the sum of magnitudes of the vector elements, single precision
    fn cblas_dasum(N: c_int, X: *const c_double, incX: c_int) -> c_double;

    /// Computes a vector-scalar product and adds the result to a vector (f32)
    fn cblas_saxpy(n: c_int,
                   alpha: c_float, X: *const c_float, incX: c_int,
                   Y: *mut c_float, incY: c_int);

    /// Computes a vector-scalar product and adds the result to a vector (f64)
    fn cblas_daxpy(n: c_int,
                   alpha: c_double, X: *const c_double, incX: c_int,
                   Y: *mut c_double, incY: c_int);

    /// Computes the Euclidean norm of a vector (f32)
    fn cblas_snrm2(n: c_int, X: *const f32, incx: c_int) -> c_float;

    /// Computes the Euclidean norm of a vector (f64)
    fn cblas_dnrm2(n: c_int, X: *const f64, incx: c_int) -> c_double;
}

/// OpenBLAS computation of the dot product of two vectors, double precision
pub fn openblas_ddot(arr_a: &Array<f64>, arr_b: &Array<f64>) -> f64
{
    if arr_a.len() != arr_b.len()
    {
        panic!("Dot product of vectors of different lengths is undefined.")
    }

    let n : c_int = arr_a.len() as c_int;
    let inc_x : c_int = 1 as c_int;
    let inc_y : c_int = 1 as c_int;
    
    unsafe {
        cblas_ddot (n,
                    arr_a.as_ptr(), inc_x,
                    arr_b.as_ptr(), inc_y)
    }
}

/// OpenBLAS computation of the dot product of two vectors, single precision
pub fn openblas_sdot(arr_a: &Array<f32>, arr_b: &Array<f32>) -> f32
{
    if arr_a.len() != arr_b.len()
    {
        panic!("Dot product of vectors of different lengths is undefined.")
    }

    let n : c_int = arr_a.len() as c_int;
    let inc_x : c_int = 1 as c_int;
    let inc_y : c_int = 1 as c_int;
    
    unsafe {
        cblas_sdot(n,
                   arr_a.as_ptr(), inc_x,
                   arr_b.as_ptr(), inc_y)
    }
}

/// OpenBLAS computation of sum of magnitudes of the vector elements (f32)
pub fn openblas_sasum (arr_a: &Array<f32>) -> f32
{
    let n : c_int = arr_a.len() as c_int;
    let inc_x : c_int = 1 as c_int;
    
    unsafe {
        cblas_sasum(n, arr_a.as_ptr(), inc_x)
    }
}

/// OpenBLAS computation of sum of magnitudes of the vector elements (f64)
pub fn openblas_dasum (arr_a: &Array<f64>) -> f64
{
    let n : c_int = arr_a.len() as c_int;
    let inc_x : c_int = 1 as c_int;
    
    unsafe {
        cblas_dasum(n, arr_a.as_ptr(), inc_x)
    }
}

/// OpenBLAS computation of vector-scalar product, adds the result to a vector (f64).
pub fn openblas_daxpy(arr_x: &Array<f64>,
                      alpha: f64,
                      arr_y: &mut Array<f64>)
{
    if arr_x.len() != arr_y.len()
    {
        panic!("Vector addition for vectors of different lengths is undefined.")
    }

    let n : c_int = arr_x.len() as c_int;
    let inc_x : c_int = 1 as c_int;
    let inc_y : c_int = 1 as c_int;

    unsafe
    {
        cblas_daxpy(n,
                    alpha as c_double,
                    arr_x.as_ptr(), inc_x,
                    arr_y.as_mut_ptr(), inc_y);
    }
}

/// OpenBLAS computation of vector-scalar product, adds the result to a vector (f32)
pub fn openblas_saxpy(arr_x: &Array<f32>,
                      alpha: f32,
                      arr_y: &mut Array<f32>)
{
    if arr_x.len() != arr_y.len()
    {
        panic!("Vector addition for vectors of different lengths is undefined.")
    }

    let n : c_int = arr_x.len() as c_int;
    let inc_x : c_int = 1 as c_int;
    let inc_y : c_int = 1 as c_int;

    unsafe
    {
        cblas_saxpy(n,
                    alpha as c_float,
                    arr_x.as_ptr(), inc_x,
                    arr_y.as_mut_ptr(), inc_y);
    }
}

/// Computes the Euclidean norm of a vector (f32)
pub fn openblas_snrm2(arr_a: &Array<f32>) -> f32
{
    let n = arr_a.len() as c_int;
    let incx = 1 as c_int;
    
    unsafe {
        return cblas_snrm2(n, arr_a.as_ptr() as *const f32, incx) as f32;
    }
}

/// Computes the Euclidean norm of a vector (f64)
pub fn openblas_dnrm2(arr_a: &Array<f64>) -> f64
{
    let n = arr_a.len() as c_int;
    let incx = 1 as c_int;
    
    unsafe {
        return cblas_dnrm2(n, arr_a.as_ptr() as *const f64, incx) as f64;
    }
}
