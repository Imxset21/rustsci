mod array;
mod matrix;

#[cfg(test)]
mod test
{
    use array;
    use matrix;
    use std::mem;

    #[test]
    fn test_arr_eq()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        assert!(a == b);
    }
    
    #[test]
    fn test_arr_sub()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        let c = array::Array::<i32>::new(vec![3, 3, 3], array::Order::Row);
        assert!(b - a == c);
    }

    #[test]
    fn test_arr_add()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        let c = array::Array::<i32>::new(vec![5, 7, 9], array::Order::Row);
        assert!(a + b == c);
    }

    #[test]
    fn test_arr_memuse()
    {
        let a = array::Array::<i32>::new(vec![0; 100], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        assert!(mem::size_of_val(&a) == mem::size_of_val(&b));
    }

    #[test]
    fn test_mat_add()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![5, 6], vec![7, 8]]);
        let c = matrix::Matrix::<i32>::new(vec![vec![6, 8], vec![10, 12]]);
        assert!(a + b == c);
    }

    #[test]
    fn test_mat_sub()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2], vec![3, 4]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![5, 6], vec![7, 8]]);
        let c = matrix::Matrix::<i32>::new(vec![vec![-4, -4], vec![-4, -4]]);
        assert!(a - b == c);
    }

    #[test]
    fn test_mat_eq()
    {
        let a = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![3, 4, 5]]);
        let b = matrix::Matrix::<i32>::new(vec![vec![1, 2, 3], vec![3, 4, 5]]);
        assert!(a == b);
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
}
