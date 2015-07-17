#![feature(step_by)]
#![feature(libc)]
#![allow(dead_code)]

extern crate libc;

#[macro_use] mod array;
#[macro_use] mod matrix;
mod common;
mod lapacke;

#[cfg(test)]
mod test
{
    use array;
    use matrix;
    use lapacke;
    use common::Transposable;
    use std::mem;

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
    fn test_arr_dot_product_generic()
    {
        let a = arr![1, 2, 3];
        let b = arr![4, 5, 6];
        assert_eq!(32, a * b);
    }
}

