use matrix::Matrix;
use array::Array;
use array::Order;
use std::cmp::min;
use libc::{c_int, c_double, c_float, c_char};

#[link(name = "lapack")]
extern
{
    fn dpotrf_(uplo: *mut c_char,
               n: *mut c_int,
               a: *mut c_double,
               lda: *mut c_int,
               info: *mut c_int);

    /// SGEBRD reduces a general real M-by-N matrix A to upper or lower
    /// bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
    /// If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    fn sgebrd_(m: *mut c_int, n: *mut c_int,
               a: *mut c_float, lda: *mut c_int, d: *mut c_float, e: *mut c_float,
               tauq: *mut c_float, taup: *mut c_float, work: *mut c_float,
               lwork: *mut c_int, info: *mut c_int);

    // SBDSQR computes the singular values and, optionally, the right and/or
    // left singular vectors from the singular value decomposition (SVD) of a
    // real N-by-N (upper or lower) bidiagonal matrix B using the implicit
    // zero-shift QR algorithm.  The SVD of B has the form
    //    B = Q * S * P**T
    // where S is the diagonal matrix of singular values, Q is an orthogonal
    // matrix of left singular vectors, and P is an orthogonal matrix of right
    // singular vectors.  If left singular vectors are requested, this
    // subroutine actually returns U*Q instead of Q, and, if right singular
    // vectors are requested, this subroutine returns P**T*VT instead of P**T,
    // for given real input matrices U and VT.  When U and VT are the orthogonal
    // matrices that reduce a general matrix A to bidiagonal form: A = U*B*VT,
    // as computed by SGEBRD, then
    //    A = (U*Q) * S * (P**T*VT)
    // is the SVD of A.  Optionally, the subroutine may also compute Q**T*C for
    // a given real input matrix C.
    fn sbdsqr_(uplo: *mut c_char, n: *mut c_int, ncvt: *mut c_int,
               nru: *mut c_int, ncc: *mut c_int, d: *mut c_float, e: *mut c_float,
               vt: *mut c_float, ldvt: *mut c_int, u: *mut c_float, ldu: *mut c_int,
               c: *mut c_float, ldc: *mut c_int, work: *mut c_float, info: *mut c_int);

}

/// Performs a Cholesky decomposition on a matrix in-place
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
        dpotrf_(&mut uplo, &mut n, mat.as_mut_ptr(), &mut lda, &mut info);
    }
    mat.zero_trigonal_lower(0.0);
}

/// Reduces a M-by-N matrix A to a bidiagonal matrix
pub fn bidiagonal_reduction(a: &mut Matrix<f32>)
{
    let (_m, _n) = a.get_dims();
    let mut m = _m as c_int;
    let mut n = _n as c_int;
    let mut lda: c_int = n as c_int;
    // D is REAL array, dimension (min(M,N))
    // The diagonal elements of the bidiagonal matrix B:
    // D(i) = A(i,i).
    let mut d = Array::<f32>::new_filled(0f32, min(_m, _n), Order::Row);
    // E is REAL array, dimension (min(M,N)-1)
    // The off-diagonal elements of the bidiagonal matrix B:
    // if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
    // if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
    let mut e = Array::<f32>::new_filled(0f32, min(_m, _n) - 1, Order::Row);
    // TAUQ is REAL array dimension (min(M,N))
    // The scalar factors of the elementary reflectors which
    // represent the orthogonal matrix Q.
    let mut tauq = Array::<f32>::new_filled(0f32, min(_m, _n), Order::Row);
    // TAUP is REAL array, dimension (min(M,N))
    // The scalar factors of the elementary reflectors which
    // represent the orthogonal matrix P.
    let mut taup = Array::<f32>::new_filled(0f32, min(_m, _n), Order::Row);
    // WORK is REAL array, dimension (M + N)
    let mut lwork: c_int = m + n;
    let mut work = Array::<f32>::new_filled(0f32, _m + _n, Order::Row);
    // INFO is INTEGER
    // = 0:  successful exit 
    // < 0:  if INFO = -i, the i-th argument had an illegal value.
    let mut info: c_int = 0;

    unsafe {
        sgebrd_(&mut m, &mut n,
                a.as_mut_ptr(), &mut lda, d.as_mut_ptr(), e.as_mut_ptr(),
                tauq.as_mut_ptr(), taup.as_mut_ptr(), work.as_mut_ptr(),
                &mut lwork, &mut info);
    }
}

// Computes the singular values and, optionally, the right and/or left singular
// vectors from the singular value decomposition (SVD) of a real N-by-N (upper
// or lower) bidiagonal matrix B using the implicit zero-shift QR algorithm.
pub fn q_svd_f32(b: &mut Matrix<f32>)
{
    if !b.is_square()
    {
        panic!("Cannot perform SVD on non-square matrix.");
    }
    unimplemented!();
}
