mod array;

#[cfg(test)]
mod test
{
    use array;
    use std::mem;

    #[test]
    fn test_eq()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        assert!(a == b);
    }
    
    #[test]
    fn test_sub()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        let c = array::Array::<i32>::new(vec![3, 3, 3], array::Order::Row);
        assert!(b - a == c);
    }

    #[test]
    fn test_add()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        let c = array::Array::<i32>::new(vec![5, 7, 9], array::Order::Row);
        assert!(a + b == c);
    }

    #[test]
    fn test_memuse()
    {
        let a = array::Array::<i32>::new(vec![0; 100], array::Order::Row);
        let b = array::Array::<i32>::new(vec![4, 5, 6], array::Order::Row);
        assert!(mem::size_of_val(&a) == mem::size_of_val(&b));
    }
}
