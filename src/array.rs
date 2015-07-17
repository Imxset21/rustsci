// Array type

// mod common;  // <- Uncomment for static analysis tools

use common::Transposable;

use std::ops::{Add, Sub, Index, IndexMut, Mul};
use std::cmp::PartialEq;

/// Enumerator for whether the vector is horizontal or vertical
#[derive(Debug, PartialEq)]
pub enum Order { Column, Row}

#[derive(Debug)]
/// An array is a boxed vector with additional information (e.g. stride).
/// Arrays support generics as long as they meet the following criteria:
/// 1. 
pub struct Array<T> where T: Add + Sub + Copy + PartialEq
{
    my_vec: Vec<T>,
    order: Order,
}

impl<T> Array<T> where T: Add + Sub + Copy + PartialEq
{
    /// A public constructor
    pub fn new(contents: Vec<T>, order: Order) -> Array<T>
    {
        Array
        {
            // I guess this is a copy constructor?
            my_vec: contents,
            order: order
        }
    }
}

impl<T> Transposable for Array<T>
    where T: Add<Output=T> + Sub + Copy + PartialEq
{
    fn transpose(self) -> Array<T>
    {
        let mut order = self.order;
        if order == Order::Column
        {
            order = Order::Row;
        } else {
            order = Order::Column;
        }

        Array {my_vec: self.my_vec, order: order}
    }
}

/// Implementation of the Add trait
impl<T> Add for Array<T> where T: Add<Output=T> + Sub + Copy + PartialEq
{
    type Output = Array<T>;

    fn add(self, other: Array<T>) -> Array<T>
    {
        if self.order != other.order
        {
            panic!("Cannot add two Arrays of different order.");
        } else if self.my_vec.len() != other.my_vec.len() {
            panic!("Cannot add two Arrays of different sizes.");
        } else if self.my_vec.len() == 0 {
            return Array { my_vec: vec![], order: self.order }
        }

        fn add_two<'a, T>((a, b) : (&'a T, &'a T),) -> T
            where T: Add<Output=T> + Sub + Copy + PartialEq
        {
            *a + *b
        }

        // Make a zipped iterator containing (vec1[i], vec2[i]), then produce
        // a map iterator which yields their sum, and collect into a new vec.
        return Array { my_vec : self.my_vec.iter().zip(other.my_vec.iter())
                       .map(add_two::<T>).collect::<Vec<T>>(),
                       order: self.order
        }
    }
}

/// Implementation of the Subtract trait
impl<T> Sub for Array<T> where T: Add + Sub<Output=T> + Copy + PartialEq
{
    type Output = Array<T>;

    fn sub(self, other: Array<T>) -> Array<T>
    {
        if self.order != other.order
        {
            panic!("Cannot add two Arrays of different order.");
        } else if self.my_vec.len() != other.my_vec.len() {
            panic!("Cannot add two Arrays of different sizes.");
        } else if self.my_vec.len() == 0 {
            return Array { my_vec: vec![], order: self.order }
        }

        fn sub_two<'a, T>((a, b) : (&'a T, &'a T),) -> T
            where T: Add + Sub<Output=T> + Copy + PartialEq
        {
            *a - *b
        }
        
        // Make a zipped iterator containing (vec1[i], vec2[i]), then produce
        // a map iterator which yields their sub, and collect into a new vec.
        return Array { my_vec : self.my_vec.iter().zip(other.my_vec.iter())
                       .map(sub_two::<T>).collect::<Vec<T>>(),
                       order: self.order
        }
    }
}

/// Implementation of the PartialEq trait, allowing for a = b, a != b
impl<T> PartialEq for Array<T> where T: Add + Sub + Copy + PartialEq
{
    fn eq(&self, other: &Array<T>) -> bool
    {
        if self.order != other.order
        {
            panic!("Cannot add two Arrays of different order.");
        } else if self.my_vec.len() != other.my_vec.len() {
            panic!("Cannot compare arrays of different dimensions.")
        }
        
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

/// Enables getting an array value via index, e.g. a[0]
impl <T> Index<usize> for Array<T>
    where T: Add + Sub + Copy + PartialEq
{
    type Output = T;

    fn index<'a>(&'a self, index: usize) -> &'a T
    {
        if index > self.my_vec.len()
        {
            panic!("Index out of bounds.");
        }
        return &self.my_vec[index];
    }
}

/// Enables things like a[0] = 3, where a is an Array
impl <T> IndexMut<usize> for Array<T>
    where T: Add + Sub + Copy + PartialEq
{
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut T
    {
        if index > self.my_vec.len()
        {
            panic!("Index out of bounds.");
        }
        return &mut self.my_vec[index];
    }
}

/// Implements the multiplication trait as a dot product
impl <T> Mul for Array<T>
    where T: Add<Output=T> + Sub + Copy + PartialEq + Mul<Output=T>
{
    type Output = T;

    fn mul(self, _rhs: Array<T>) -> T
    {
        if self.my_vec.len() != _rhs.my_vec.len()
        {
            panic!("Cannot dot product different-length arrays.")
        }
        
        // let sum : Option<T> = None;
        // for (a_i, b_i) in self.my_vec.into_iter().zip(_rhs.my_vec.into_iter())
        // {
        //     sum = match sum {
        //         None => Some(a_i * b_i),
        //         Some(sum_val) => Some(sum_val + (a_i * b_i))
        //     };
        // }

        let sum: Option<T> =
            self.my_vec.into_iter().zip(_rhs.my_vec.into_iter()).fold(
                None,
                |sum_val, (a_i, b_i)| match sum_val {
                    None => Some(a_i * b_i),
                    Some(curr_val) => Some(curr_val + (a_i * b_i))
                }
            );
        
        // Return value
        match sum {
            Some(sum_val) => sum_val,
            None => panic!("Unable to calculate dot product")
        }
    }
}

/// Convenience macro for creating arrays in a vec!-like way
macro_rules! arr
{
    ( $( $x:expr ),* ) =>
    {
        {
            let mut temp_vec = Vec::new();
            $(
                temp_vec.push($x);
                )*
                array::Array::new(temp_vec, array::Order::Row)
        }
    };
}
