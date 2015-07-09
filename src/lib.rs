mod array;

#[cfg(test)]
mod test
{
    use array;

    #[test]
    fn test_eq()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3]);
        let b = array::Array::<i32>::new(vec![1, 2, 3]);
        assert!(a == b);
    }
    
    #[test]
    fn test_sub()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3]);
        let b = array::Array::<i32>::new(vec![4, 5, 6]);
        let c = array::Array::<i32>::new(vec![3, 3, 3]);
        assert!(b - a == c);
    }

    #[test]
    fn test_add()
    {
        let a = array::Array::<i32>::new(vec![1, 2, 3]);
        let b = array::Array::<i32>::new(vec![4, 5, 6]);
        let c = array::Array::<i32>::new(vec![5, 7, 9]);
        assert!(a + b == c);
    }

}
