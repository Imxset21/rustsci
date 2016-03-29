/// OpenBLAS testing
#[macro_use]
extern crate rustsci;

use rustsci::array;
use rustsci::matrix;
use rustsci::openblas;

////////////////////////
// BLAS Level 1 Tests //
////////////////////////

#[test]
fn test_openblas_ddot()
{
    let a1 = arr![1f64, 2f64, 3f64];
    let a2 = arr![4f64, 5f64, 6f64];
    let blas_result = openblas::openblas_ddot(&a1, &a2);
    let gen_result = a1 * a2;  // Generic dot product
    assert_eq!(gen_result, blas_result);
}

#[test]
fn test_openblas_sdot()
{
    let a1 = arr![1f32, 2f32, 3f32];
    let a2 = arr![4f32, 5f32, 6f32];
    let blas_result = openblas::openblas_sdot(&a1, &a2);
    let gen_result = a1 * a2;  // Generic dot product
    assert_eq!(gen_result, blas_result);
}

#[test]
fn test_openblas_sasum()
{
    let a1 = arr![1f32, 2f32, 3f32];
    let elem_sum = 1f32 + 2f32 + 3f32;
    let blas_result = openblas::openblas_sasum(&a1);
    assert_eq!(elem_sum, blas_result);
}

#[test]
fn test_openblas_dasum()
{
    let a1 = arr![1f64, 2f64, 3f64];
    let elem_sum = 1f64 + 2f64 + 3f64;
    let blas_result = openblas::openblas_dasum(&a1);
    assert_eq!(elem_sum, blas_result);
}

#[test]
fn test_openblas_daxpy()
{
    let a1 = arr![1f64, 2f64, 3f64];
    let mut a2 = arr![0f64, 0f64, 0f64];
    let alpha = 5f64;
    
    openblas::openblas_daxpy(&a1, alpha, &mut a2);
    
    let result_arr = arr![5f64, 10f64, 15f64];
    assert_eq!(result_arr, a2);
}

#[test]
fn test_openblas_saxpy()
{
    let a1 = arr![1f32, 2f32, 3f32];
    let mut a2 = arr![1f32, 2f32, 3f32];
    let alpha = 5f32;
    
    openblas::openblas_saxpy(&a1, alpha, &mut a2);
    
    let result_arr = arr![6f32, 12f32, 18f32];
    assert_eq!(result_arr, a2);
}

#[test]
fn test_openblas_snrm2()
{
    let arr = arr![8f32, 4f32, 1f32, 0f32];
    let norm : f32 = openblas::openblas_snrm2(&arr);
    assert_eq!(norm, 9f32);
}

#[test]
fn test_openblas_dnrm2()
{
    let arr = arr![8f64, 4f64, 1f64, 0f64];
    let norm : f64 = openblas::openblas_dnrm2(&arr);
    assert_eq!(norm, 9f64);
}

////////////////////////
// BLAS Level 2 Tests //
////////////////////////

#[test]
fn test_filled_new()
{
    let m1 = matrix::Matrix::<f32>::new_filled(0f32, 3, 3);
    let m2 = mat![[0f32, 0f32, 0f32],
                  [0f32, 0f32, 0f32],
                  [0f32, 0f32, 0f32]];
    assert_eq!(m1, m2);
}

#[test]
fn test_openblas_sgemv()
{
    let mat_a = mat![[25f32, 15f32, -5f32],
                     [15f32, 18f32, 0f32],
                     [-5f32, 0f32, 11f32]];
    let arr_x = arr![8f32, 4f32, 1f32];
    let arr_y = array::Array::<f32>::new_filled(0f32, 3, array::Order::Row);
    let alpha = 1f32;
    let beta = 1f32;

    let ymat = openblas::openblas_sgemv(&mat_a, &arr_x, &arr_y, alpha, beta);
    let resultm = arr![255f32, 192f32, -29f32];
    assert_eq!(ymat, resultm);
}

////////////////////////
// BLAS Level 3 Tests //
////////////////////////

#[test]
fn test_openblas_sgemm()
{
    let mat_a = mat![[25f32, 15f32, -5f32],
                     [15f32, 18f32,  0f32],
                     [-5f32,  0f32, 11f32]];
    let mat_b = mat![[7f32, 15f32, -12f32],
                     [1f32, -14f32,  1f32],
                     [-5f32,  1f32, -11f32]];
    let mat_c = mat![[1f32, 0f32, 0f32],
                     [0f32, 1f32, 0f32],
                     [0f32, 0f32, 1f32]];
    let alpha = 1f32;
    let beta = 1f32;

    let c_mat = openblas::openblas_sgemm(&mat_a, &mat_b, &mat_c, alpha, beta);

    let result_mat = mat![[216f32, 160f32, -230f32],
                          [123f32, -26f32, -162f32],
                          [-90f32, -64f32, -60f32]];

    assert_eq!(c_mat, result_mat);
}
