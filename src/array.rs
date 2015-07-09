/// Array type

use std::ops::Add;
use std::ops::Sub;
use std::cmp::PartialEq;

#[derive(Debug)]
/// Wraps an array in a box
pub struct Array<T: Add + Sub + Copy + PartialEq>
{
    my_vec: Box<Vec<T>>,
}

impl<T: Add + Copy + Sub + PartialEq> Array<T>
{
    /// A public constructor
    pub fn new(contents: Vec<T>) -> Array<T>
    {
        Array
        {
            // I guess this is a copy constructor?
            my_vec: Box::new(contents),
        }
    }
}

impl<T: Add<Output=T> + Sub + Copy + PartialEq> Add for Array<T>
{
    type Output = Array<T>;

    fn add(self, other: Array<T>) -> Array<T>
    {
        if self.my_vec.len() != other.my_vec.len()
        {
            panic!("Cannot add two Arrays of different sizes.");
        } else if self.my_vec.len() == 0 {
            return Array { my_vec: Box::<Vec<T>>::new(vec![]) }
        }
        
        let mut b = Vec::with_capacity(self.my_vec.len());
        for (num1, num2) in self.my_vec.iter().zip(other.my_vec.iter())
        {
            b.push(*num1 + *num2);
        }
        
        return Array { my_vec : Box::<Vec<T>>::new(b)}
    }
}

impl<T: Add + Sub<Output=T> + Copy + PartialEq> Sub for Array<T>
{
    type Output = Array<T>;

    fn sub(self, other: Array<T>) -> Array<T>
    {
        if self.my_vec.len() != other.my_vec.len()
        {
            panic!("Cannot add two Arrays of different sizes.");
        } else if self.my_vec.len() == 0 {
            return Array { my_vec: Box::<Vec<T>>::new(vec![]) }
        }
        
        let mut b = Vec::with_capacity(self.my_vec.len());
        for (num1, num2) in self.my_vec.iter().zip(other.my_vec.iter())
        {
            b.push(*num1 - *num2);
        }
        
        return Array { my_vec : Box::<Vec<T>>::new(b)}
    }
}

impl<T: Add + Sub + Copy + PartialEq> PartialEq for Array<T>
{
    fn eq(&self, other: &Array<T>) -> bool
    {
        for (num1, num2) in self.my_vec.iter().zip(other.my_vec.iter())
        {
            if *num1 != *num2
            {
                return false
            }
        }

        return true
    }
}
