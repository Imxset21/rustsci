// #![feature(libc)]   // <- Uncomment for static analysis
// extern crate libc; // <- Uncomment for static analysis
// mod array;  // <- Uncomment for static analysis

use array::Array;
use libc::{c_int, c_double};

#[link(name = "blas")]
extern
{
    /// Computes a vector-vector dot product for f64 ("double") precision.
    fn cblas_ddot (N: c_int,
                   X: *const c_double, incX: c_int,
                   Y: *const c_double, incY: c_int) -> c_double;
}

/// OpenBLAS implementation of the dot product of two vectors, double precision
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
