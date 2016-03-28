#![feature(step_by)]
#![feature(libc)]
#![allow(dead_code)]

extern crate libc;

#[macro_use] pub mod array;
#[macro_use] pub mod matrix;
pub mod lapacke;
pub mod openblas;
pub mod gsl_poly;
#[macro_use] pub mod gsl_math;
pub mod gsl_sf;
pub mod gsl_airy;
pub mod gsl_bessel;

#[cfg(test)]
mod test
{
    use array;
    use matrix;
    use lapacke;
    use openblas;
    use gsl_poly;
    use gsl_math;
    use gsl_airy;
    use std::mem;

    /////////////////
    // Array Tests //
    /////////////////
    
    #[test]
    fn test_arr_eq()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        assert_eq!(a, b);
    }
    
    #[test]
    fn test_arr_sub()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        let c = array::Array::<i32>::new(vec![3, 3, 3], array::Order::Row);
        assert_eq!(b - a, c);
    }

    #[test]
    fn test_arr_add()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        let c = array::Array::<i32>::new(vec![5, 7, 9], array::Order::Row);
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_arr_memuse()
    {
        let a = array::Array::<i32>::new(vec![0; 100], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        assert_eq!(mem::size_of_val(&a), mem::size_of_val(&b));
    }

    #[test]
    fn test_arr_indx()
    {
        let mut a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        a[0] = 5;
        assert_eq!(a[0], 5);
    }

    #[test]
    fn test_arr_macro()
    {
        let a = arr![1, 2, 3];
        let b = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        assert_eq!(a, b);
    }

    //////////////////
    // Matrix Tests //
    //////////////////
    
    #[test]
    fn test_mat_add()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![5, 6], vec![7, 8]]);
        let c = matrix::Matrix::<i32>::new(vec![vec![6, 8], vec![10, 12]]);
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_mat_sub()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![5, 6], vec![7, 8]]);
        let c = matrix::Matrix::<i32>::new(vec![vec![-4, -4], vec![-4, -4]]);
        assert_eq!(a - b, c);
    }

    #[test]
    fn test_mat_eq()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![3, 4, 5]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![3, 4, 5]]);
        assert_eq!(a, b);
    }

    #[test]
    #[should_panic]
    fn test_mat_invalid_new()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4, 6]]);
        assert!(a == a); // Unreachable
    }

    #[test]
    #[should_panic]
    fn test_mat_invalid_add()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![3, 4, 5]]);
        let result = a + b;
        assert!(result == result);  // Unreachable
    }

    #[test]
    fn test_mat_t_2x2()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        let a_t = matrix::Matrix::<i32>::new(vec![vec![1, 3], vec![2, 4]]);
        assert_eq!(a_t, a.transpose())
    }

    #[test]
    fn test_mat_t_2x3()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![3, 4, 5]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![1,3], vec![2, 4], vec![3, 5]]);
        assert_eq!(b, a.transpose());
    }

    #[test]
    fn test_mat_t_3x2()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4], vec![5, 6]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![4, 5, 6]]);
        assert_eq!(b, a.transpose());
    }

    #[test]
    fn test_mat_zero_tri_lower()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3],
                                                vec![0, 5, 6],
                                                vec![0, 0, 9]]);
        let mut b = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3],
                                                    vec![4, 5, 6],
                                                    vec![7, 8, 9]]);
        matrix::Matrix::<i32>::zero_trigonal_lower(&mut b, 0);
        assert_eq!(a, b);
    }

    #[test]
    fn test_square_mat_macro()
    {
        let n = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![4, 5, 6]]);
        let m = mat!([1, 2, 3], [4, 5, 6]);
        assert_eq!(m, n);
    }

    #[test]
    fn test_mat_indx_get()
    {
        let b = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        assert_eq!(b[(1, 1)], 4);
    }

    #[test]
    fn test_mat_indx_set()
    {
        let mut b = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        b[(1, 0)] = 5;
        assert_eq!(b[(1, 0)], 5);
    }

    #[test]
    fn test_mat_set()
    {
        let mut b = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3],
                                                    vec![4, 5, 6],
                                                    vec![7, 8, 9]]);
        b.set(1, 1, 0);
        assert_eq!(*b.get(1, 1), 0);
    }

    #[test]
    fn test_mat_get()
    {
        let b = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3],
                                                vec![4, 5, 6],
                                                vec![7, 8, 9]]);
        assert_eq!(*b.get(1, 1), 5);
    }

    #[test]
    fn test_mat_debug()
    {
        let m1 = mat![[1, 2, 3],
                      [4, 5, 6],
                      [7, 8, 9]];
        println!("{:?}", m1);
    }

    #[test]
    fn test_arr_dot_product_generic()
    {
        let a = arr![1, 2, 3];
        let b = arr![4, 5, 6];
        assert_eq!(32, a * b);
    }

    #[test]
    fn test_matrix_mul_generic_square()
    {
        let m1 = mat![[1, 2],
                      [3, 4]];
        let m2 = mat![[5, 6],
                      [7, 8]];
        let m3 = mat![[19, 22],
                      [43, 50]];
        assert_eq!(m3, m1 * m2);
    }

    #[test]
    fn test_matrix_mul_generic_nonsquare()
    {
        let m1 = mat![[1, 2, 3],
                      [4, 5, 6]];
        let m2 = mat![[7, 8],
                      [9, 10],
                      [11, 12]];
        let m3 = mat![[58, 64],
                      [139, 154]];
        assert_eq!(m3, m1 * m2);
    }

    ////////////////////////////
    // Symmetric Matrix Tests //
    ////////////////////////////

    #[test]
    fn test_symmat_new_diag_with_fill()
    {        
        let m2 = matrix::SymMat::<f32>::new_diag_with_fill(1f32, 0f32, 1);
        assert_eq!(m2[(0,0)], 1f32);
        let m1 = matrix::SymMat::<f32>::new_diag_with_fill(1f32, 0f32, 3);
        assert_eq!(m1[(0,0)], 1f32);
        assert_eq!(m1[(0,1)], 0f32);
        assert_eq!(m1[(0,2)], 0f32);
        assert_eq!(m1[(1,0)], 0f32);
        assert_eq!(m1[(1,1)], 1f32);
        assert_eq!(m1[(1,2)], 0f32);
        assert_eq!(m1[(2,0)], 0f32);
        assert_eq!(m1[(2,1)], 0f32);
        assert_eq!(m1[(2,2)], 1f32);
    }

    #[test]
    fn test_symmat_new_filled()
    {
        let m1 = matrix::SymMat::<f32>::new_filled(1f32, 3);
        let m2 = mat![[1f32, 1f32, 1f32],
                      [0f32, 1f32, 1f32],
                      [0f32, 0f32, 1f32]];
        let expected = matrix::SymMat::<f32>::new_from_upper_trig(&m2);
        assert!(m1 == expected);
    }

    #[test]
    fn test_symmat_convert()
    {
        let m1 = matrix::SymMat::<f32>::new_diag_with_fill(1f32, 0f32, 3);
        let expected = mat![[1f32, 0f32, 0f32],
                            [0f32, 1f32, 0f32],
                            [0f32, 0f32, 1f32]];
        let result = m1.to_mat();
        assert_eq!(expected, result);
    }

    #[test]
    fn test_symmat_upper_trig()
    {
        let m1 = mat![[1f32, 2f32, 3f32],
                      [0f32, 1f32, 4f32],
                      [0f32, 0f32, 1f32]];
        let sym = matrix::SymMat::<f32>::new_from_upper_trig(&m1);
        assert_eq!(sym[(0,1)], m1[(0,1)]);
    }

    ///////////////////
    // LAPACKE Tests //
    ///////////////////

    #[test]
    fn test_mat_cholesky()
    {
        let mut m = mat![[25., 15., -5.],
                         [15., 18., 0.],
                         [-5., 0., 11.]];
        let result = mat![[5., 3., -1.],
                          [0., 3., 1.],
                          [0., 0., 3.]];
        lapacke::cholesky_decomposition(&mut m);
        assert_eq!(m, result);
    }

    #[test]
    fn test_bidiag_reduc()
    {
        let mut m = mat![[12f32, -51f32, 4f32],
                         [6f32, 167f32, -68f32],
                         [-4f32, 24f32, -41f32]];
        lapacke::bidiagonal_reduction(&mut m);
        // TODO: Finish testing bidiag_reduc
    }

    ////////////////////
    // OpenBLAS Tests //
    ////////////////////
    
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

    //////////////////////////
    // GSL Polynomial Tests //
    //////////////////////////

    #[test]
    fn test_poly_eval()
    {
        // f(x) = 1 + x + x^2
        let c = arr![1f64, 1f64, 1f64];
        let result = gsl_poly::poly_eval(&c, 1.0);
        assert_eq!(result, 3f64);
    }

    #[test]
    fn test_poly_eval_derivs()
    {
        // f(x) = 1 + x + x^2
        let c = arr![1f64, 1f64, 1f64];
        let result = gsl_poly::poly_eval_derivs(&c, 1.0, 1);
        assert_eq!(result[0], 3f64);
    }

    #[test]
    fn test_poly_divdiff_init()
    {
        let xa = arr![1f64, 2f64, 3f64];
        let ya = arr![1f64, 2f64, 3f64];
        let poly = gsl_poly::poly_divdiff_init(&xa, &ya);
        assert_eq!(poly, arr![1f64, 1f64, 0f64]);
    }

    #[test]
    fn test_poly_divdiff_eval()
    {
        let xa = arr![1f64, 2f64, 3f64];
        let ya = arr![1f64, 2f64, 3f64];
        let poly = gsl_poly::poly_divdiff_init(&xa, &ya);
        let result = gsl_poly::poly_divdiff_eval(&poly, &xa, 1f64);
        assert_eq!(1f64, result);
    }

    #[test]
    fn test_poly_divdiff_to_taylor()
    {
        let xa = arr![1f64, 2f64, 3f64];
        let ya = arr![1f64, 2f64, 3f64];
        let poly = gsl_poly::poly_divdiff_init(&xa, &ya);
        let taylor = gsl_poly::poly_divdiff_to_taylor(0f64, &poly, &xa);
        assert_eq!(arr![0f64, 1f64, 0f64], taylor);
    }

    #[test]
    fn test_poly_divdiff_to_hermite()
    {
        // TODO: Write test for Hermite polynomial
    }

    #[test]
    fn test_poly_solve_quadratic()
    {
        let coeffs: [f64; 3] = [-4f64, -3f64, 1f64];
        let roots = gsl_poly::poly_solve_quadratic(coeffs);
        let y = match roots {
            Some(x) => x,
            None    => panic!("Test failed: no roots found"),
        };
        assert_eq!(y[0], -1f64);
        assert_eq!(y[1], 0.25f64);
    }

    #[test]
    fn test_poly_solve_cubic()
    {
        let coeffs: [f64; 3] = [-4f64, 1f64, 6f64];
        let y = gsl_poly::poly_solve_cubic(coeffs);
        assert!(gsl_math::gslmath_fcmp(y[0], -1f64, 0.001f64));
        assert!(gsl_math::gslmath_fcmp(y[1], 2f64, 0.001f64));
        assert!(gsl_math::gslmath_fcmp(y[2], 3f64, 0.001f64));
    }

    ///////////////////////////////////////
    // Special Functions: Airy Functions //
    ///////////////////////////////////////

    #[test]
    fn test_gslairy_ai()
    {
        let (val, _) = gsl_airy::gslairy_ai(1f64);
        assert!(gsl_math::gslmath_fcmp(val, 0.13529241631288141f64, 0.001f64));
    }
}


