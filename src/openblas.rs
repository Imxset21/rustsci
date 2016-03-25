use array::Array;
use matrix::Matrix;
use libc::{c_int, c_double, c_float};

enum CblasOrder {CblasRowMajor=101, CblasColMajor=102}
enum CblasTranspose {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113}
enum CblasUplo {CblasUpper=121, CblasLower=122}
enum CblasDiag {CblasNonUnit=131, CblasUnit=132}
enum CblasSide {CblasLeft=141, CblasRight=142}

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
    
    /// Computes the sum of magnitudes of the vector elements, double precision
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

    /// Matrix vector multiply (f32) - alpha * A * x + beta * y
    fn cblas_sgemv(
        order: c_int,
        TransA: c_int, M: c_int, N: c_int,
        alpha: c_float, A: *const f32, lda: c_int,
        X: *const f32, incX: c_int, beta: c_float,
        Y: *mut f32, incY: c_int);

    /// Matrix vector multiply (f64) - alpha * A * x + beta * y
    fn cblas_dgemv(
        order: c_int,
        TransA: c_int, M: c_int, N: c_int,
        alpha: c_double, A: *const f64, lda: c_int,
        X: *const f64, incX: c_int, beta: c_double,
        Y: *mut f64, incY: c_int);

    /// Matrix-Matrix multiply (f32) - alpha * A * B + beta * C
    fn cblas_sgemm(
        order: c_int, TransA: c_int, TransB: c_int,
        M: c_int, N: c_int, K: c_int,
        alpha: c_float, A: *const f32,
        lda: c_int, B: *const f32, ldb: c_int,
        beta: c_float, C: *mut f32, ldc: c_int
    );

    /// Matrix-Matrix multiply (f64) - alpha * A * B + beta * C
    fn cblas_dgemm(
        order: c_int, TransA: c_int, TransB: c_int,
        M: c_int, N: c_int, K: c_int,
        alpha: c_double, A: *const f64,
        lda: c_int, B: *const f64, ldb: c_int,
        beta: c_double, C: *mut f64, ldc: c_int
    );
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

/// Computes alpha * A * x + beta * y (f32)
/// A is an m-by-n matrix, x and y are vectors, alpha and beta are scalars.
pub fn openblas_sgemv(
    mat_a: &Matrix<f32>,
    arr_x: &Array<f32>,
    arr_y: &Array<f32>,
    alpha: f32,
    beta: f32) -> Array<f32>
{
    let (nrows, ncols) = mat_a.get_dims();
    // TODO: Fix this check to actually check the right dimensions
    if nrows != arr_x.len()
    {
        panic!("Mistmatched matrix and array dimensions: {} vs {}", nrows, arr_x.len());
    }
    let mut m1 = arr_y.clone();

    unsafe {
        cblas_sgemv(
            CblasOrder::CblasRowMajor as c_int,
            CblasTranspose::CblasNoTrans as c_int, nrows as c_int, ncols as c_int,
            alpha as c_float, mat_a.as_ptr() as *const f32, nrows as c_int,
            arr_x.as_ptr() as *const f32, 1 as c_int, beta as c_float,
            m1.as_mut_ptr() as *mut f32, 1 as c_int);
    }
    return m1;
}

/// Computes alpha * A * x + beta * y (f64)
/// A is an m-by-n matrix, x and y are vectors, alpha and beta are scalars.
pub fn openblas_dgemv(
    mat_a: &Matrix<f64>,
    arr_x: &Array<f64>,
    arr_y: &Array<f64>,
    alpha: f64,
    beta: f64) -> Array<f64>
{
    let (nrows, ncols) = mat_a.get_dims();
    // TODO: Fix this check to actually check the right dimensions
    if nrows != arr_x.len()
    {
        panic!("Mistmatched matrix and array dimensions: {} vs {}",
               nrows, arr_x.len());
    }
    let mut m1 = arr_y.clone();

    unsafe {
        cblas_dgemv(
            CblasOrder::CblasRowMajor as c_int,
            CblasTranspose::CblasNoTrans as c_int, nrows as c_int, ncols as c_int,
            alpha as c_double, mat_a.as_ptr() as *const f64, nrows as c_int,
            arr_x.as_ptr() as *const f64, 1 as c_int, beta as c_double,
            m1.as_mut_ptr() as *mut f64, 1 as c_int);
    }
    return m1;
}

/// Computes alpha * A * B + beta * C (f32)
/// Alpha and beta are scalars, and A, B and C are matrices, with A an m-by-k
/// matrix, B a k-by-n matrix and C an m-by-n matrix.
pub fn openblas_sgemm(
    mat_a: &Matrix<f32>,
    mat_b: &Matrix<f32>,
    mat_c: &Matrix<f32>,
    alpha: f32,
    beta: f32) -> Matrix<f32>
{
    let (a_nrows, a_ncols) = mat_a.get_dims();
    let (b_nrows, b_ncols) = mat_b.get_dims();
    let (c_nrows, c_ncols) = mat_c.get_dims();

    if a_nrows != c_nrows
    {
        panic!("A's number of rows doesn't match C's number of rows: {} != {}",
               a_nrows, c_nrows);
    }
    if a_ncols != b_nrows
    {
        panic!("B's number of rows doesn't match A's number of cols: {} != {}",
               b_nrows, a_ncols);
    }
    if b_ncols != c_ncols
    {
        panic!("C's number of rows doesn't match B's number of cols: {} != {}",
               c_nrows, b_ncols);
    }

    let mut mat_c_out = mat_c.clone();

    unsafe {
        cblas_sgemm(
            CblasOrder::CblasRowMajor as c_int,
            CblasTranspose::CblasNoTrans as c_int,
            CblasTranspose::CblasNoTrans as c_int,
            a_nrows as c_int, b_ncols as c_int, b_nrows as c_int,
            alpha as c_float,
            mat_a.as_ptr() as *const f32, a_nrows as c_int,
            mat_b.as_ptr() as *const f32, b_nrows as c_int,
            beta as c_float,
            mat_c_out.as_mut_ptr() as *mut f32, c_nrows as c_int
        );
    }
    return mat_c_out;
}

/// Computes alpha * A * B + beta * C (f64)
/// Alpha and beta are scalars, and A, B and C are matrices, with A an m-by-k
/// matrix, B a k-by-n matrix and C an m-by-n matrix.
pub fn openblas_dgemm(
    mat_a: &Matrix<f64>,
    mat_b: &Matrix<f64>,
    mat_c: &Matrix<f64>,
    alpha: f64,
    beta: f64) -> Matrix<f64>
{
    let (a_nrows, a_ncols) = mat_a.get_dims();
    let (b_nrows, b_ncols) = mat_b.get_dims();
    let (c_nrows, c_ncols) = mat_c.get_dims();

    if a_nrows != c_nrows
    {
        panic!("A's number of rows doesn't match C's number of rows: {} != {}",
               a_nrows, c_nrows);
    }
    if a_ncols != b_nrows
    {
        panic!("B's number of rows doesn't match A's number of cols: {} != {}",
               b_nrows, a_ncols);
    }
    if b_ncols != c_ncols
    {
        panic!("C's number of rows doesn't match B's number of cols: {} != {}",
               c_nrows, b_ncols);
    }

    let mut mat_c_out = mat_c.clone();

    unsafe {
        cblas_dgemm(
            CblasOrder::CblasRowMajor as c_int,
            CblasTranspose::CblasNoTrans as c_int,
            CblasTranspose::CblasNoTrans as c_int,
            a_nrows as c_int, b_ncols as c_int, b_nrows as c_int,
            alpha as c_double,
            mat_a.as_ptr() as *const f64, a_nrows as c_int,
            mat_b.as_ptr() as *const f64, b_nrows as c_int,
            beta as c_double,
            mat_c_out.as_mut_ptr() as *mut f64, c_nrows as c_int
        );
    }
    return mat_c_out;
}
