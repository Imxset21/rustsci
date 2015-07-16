// Array type

// mod common;  // <- Uncomment for static analysis tools

use common::Transposable;

use std::ops::{Add, Sub, Index, IndexMut};
use std::cmp::PartialEq;

/// Enumerator for whether the vector is horizontal or vertical
#[allow(dead_code)]
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
    #[allow(dead_code)]
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

/// Implementation of the PartialEq trait
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
