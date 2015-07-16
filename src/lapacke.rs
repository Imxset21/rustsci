// #![feature(libc)]   // <- Uncomment for static analysis
// extern crate libc; // <- Uncomment for static analysis
// mod matrix;  // <- Uncomment for static analysis

use matrix::Matrix;
use libc::{c_int, c_double, c_char};
// use std::mem;

#[link(name = "lapack")]
extern
{
    fn dpotrf_(uplo: *mut c_char,
               n: *mut c_int,
               a: *mut c_double,
               lda: *mut c_int,
               info: *mut c_int);
}

pub fn cholesky_decomposition(mat: &mut Matrix<f64>)
{
    if !mat.is_square()
    {
        panic!("Cannot perform Cholesky decomposition on non-square matrix.");
    }
    let mut uplo : c_char = 'L' as c_char;
    let (_n, _) = mat.get_dims();
    let mut info : c_int = 0;
    let mut n : c_int = _n as c_int;
    let mut lda : c_int = n;

    // Unsafe due to being an FFI, and because of .as_mut_ptr()
    unsafe {
        // let vec_cpy : *mut c_double =
        //     libc::calloc(mem::size_of::<c_double>() as size_t, n as size_t)
        //     as *mut c_double;
        // dpotrf_(&mut uplo, &mut n, vec_cpy, &mut lda, &mut info);
        dpotrf_(&mut uplo, &mut n, mat.as_mut_ptr(), &mut lda, &mut info);
    }
    mat.zero_trigonal_lower(0.0);
}
