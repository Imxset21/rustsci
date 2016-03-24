// Array type

use std::ops::{Add, Sub, Index, IndexMut, Mul};
use std::cmp::PartialEq;
use std::slice;

/// Enumerator for whether the vector is horizontal or vertical
#[derive(Debug, PartialEq, Copy, Clone)]
pub enum Order { Column, Row}

#[derive(Debug, Clone)]
/// An array is a size-static vector with type restrictions.
pub struct Array<T> where T: Add + Sub + Copy + PartialEq
{
    /// Internal vector used for backing store
    my_vec: Vec<T>,
    /// Order of the array
    order: Order,
}

/// Generic implementation for any T that satisfies copying, add, sub, and eq
impl<T> Array<T> where T: Add + Sub + Copy + PartialEq
{
    /// Default public constructor, moves the contents from the vec by default
    pub fn new(contents: Vec<T>, order: Order) -> Array<T>
    {
        Array
        {
            my_vec: contents,
            order: order
        }
    }

    /// Constructs an array with a given value, length, and order
    pub fn new_filled(val: T, len: usize, order: Order) -> Array<T>
    {
        Array
        {
            my_vec: vec![val; len],
            order: order
        }
    }

    /// Transposes an array by flipping its orientation.
    pub fn transpose(self) -> Array<T>
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

    /// Returns a slice (in-memory view) of the entire array
    fn as_slice(&self) -> &[T]
    {
        return self.my_vec.as_slice();
    }

    /// Returns an iterator over size elements of the slice at a time. The
    /// chunks do not overlap. If size does not divide the length of the slice,
    /// then the last chunk will not have length size.
    fn chunks(&self, size: usize) -> slice::Chunks<T>
    {
        return self.my_vec.chunks(size);
    }
}

/// Implementation of the Add trait.
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

        #[inline]
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

/// Implementation of the Subtract trait.
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

/// Enables setting an array value via index, e.g. a[0] = 3
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

/// Implements the multiplication trait as a generic dot product
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

impl <T> Array<T> where T: Add + Sub + PartialEq + Copy
{
    /// Gets array contents as raw pointer
    pub unsafe fn as_ptr(&self) -> *const T
    {
        self.my_vec.as_ptr()
    }

    /// Gets array contents as a mutable pointer
    pub unsafe fn as_mut_ptr(&mut self) -> *mut T
    {
        self.my_vec.as_mut_ptr()
    }

    /// Gets length of underlying vector
    pub fn len(&self) -> usize
    {
        self.my_vec.len()
    }
}

/// Convenience macro for creating arrays in a vec!-like way
#[macro_export]
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
