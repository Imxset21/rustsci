// Matrix type

use std::ops::Add;
use std::ops::Sub;
use std::cmp::PartialEq;

/// Regular ol' matrix. No sparseness or any of that garbage, that's for later.
#[derive(Debug)]
pub struct Matrix<T> where T: Add + Sub + Copy + PartialEq
{
    my_dat: Vec<T>,
    col_ordering: bool,
    num_rows: usize,
    num_cols: usize,
}

/// Private (?) helper function
#[allow(dead_code)]
fn get_content_dims<T>(contents: &Vec<Vec<T>>) -> (usize, usize)
    where T: Add + Sub + Copy + PartialEq
{
    let num_rows = contents.len();
    if num_rows == 0
    {
        return (0, 0)
    } else {
        return (num_rows, contents[0].len());
    }
}

/// Validates whether a vector of vectors is consistent
#[allow(dead_code)]
fn validate_contents<T>(contents: &Vec<Vec<T>>) -> bool
    where T: Add + Sub + Copy + PartialEq
{
    let num_rows = contents.len();
    if num_rows == 0
    {
        return false;
    } else {
        let num_cols = contents[0].len();
        for vec in contents.iter()
        {
            if vec.len() != num_cols
            {
                return false;
            }
        }
        return true;
    }
}

impl<T> Matrix<T> where T: Add + Sub + Copy + PartialEq
{
    #[allow(dead_code)]
    pub fn new(contents: Vec<Vec<T>>) -> Matrix<T>
    {
        if !validate_contents(&contents)
        {
            panic!("Matrix contents are invalid: dimensions inconsistent.")
        }
        let (num_rows, num_cols) = get_content_dims(&contents);

        // Flattened iterator that goes through the entire matrix at once
        let flattened = contents.iter().flat_map(|y| y.iter());

        let mut new_vec = Vec::<T>::with_capacity(num_rows * num_cols);
        for i in flattened
        {
            new_vec.push(*i);
        }

        // Return new matrix
        Matrix
        {
            my_dat: new_vec,
            col_ordering: true,
            num_rows: num_rows,
            num_cols: num_cols,
        }
    }

    #[allow(dead_code)]    
    pub fn new_from_vec(contents: Vec<T>,
                        num_rows: usize,
                        num_cols: usize) -> Matrix<T>
    {
        if contents.len() % (num_rows * num_cols) != 0
        {
            panic!("Mismatch between matrix contents and dimensions.");
        }

        Matrix
        {
            my_dat: contents,
            col_ordering: true,
            num_rows: num_rows,
            num_cols: num_cols,
        }
    }
}

/// Implementation of the Add trait
impl<T> Add for Matrix<T> where T: Add<Output=T> + Sub + Copy + PartialEq
{
    type Output = Matrix<T>;

    fn add(self, other: Matrix<T>) -> Matrix<T>
    {
        if self.num_rows != other.num_rows ||
            self.num_cols != other.num_cols
        {
            panic!("Cannot add two Matrixs of different dimensions.");
        }

        fn add_two<'a, T>((a, b) : (&'a T, &'a T),) -> T
            where T: Add<Output=T> + Sub + Copy + PartialEq
        {
            *a + *b
        }
        
        // Make a zipped iterator containing (vec1[i], vec2[i]), then produce
        // a map iterator which yields their sum, and collect into a new vec.
        let my_vec = self.my_dat.iter().zip(other.my_dat.iter())
                       .map(add_two::<T>).collect::<Vec<T>>();

        Matrix { my_dat: my_vec, .. self }
    }
}

/// Implementation of the Add trait
impl<T> Sub for Matrix<T> where T: Add + Sub<Output=T> + Copy + PartialEq
{
    type Output = Matrix<T>;

    fn sub(self, other: Matrix<T>) -> Matrix<T>
    {
        if self.num_rows != other.num_rows ||
            self.num_cols != other.num_cols
        {
            panic!("Cannot add two Matrixs of different dimensions.");
        }

        fn sub_two<'a, T>((a, b) : (&'a T, &'a T),) -> T
            where T: Add + Sub<Output=T> + Copy + PartialEq
        {
            *a - *b
        }
        
        // Make a zipped iterator containing (vec1[i], vec2[i]), then produce
        // a map iterator which yields their sum, and collect into a new vec.
        let my_vec = self.my_dat.iter().zip(other.my_dat.iter())
                       .map(sub_two::<T>).collect::<Vec<T>>();

        Matrix { my_dat: my_vec, .. self }
    }
}


impl<T> PartialEq for Matrix<T> where T: Add<Output=T> + Sub + Copy + PartialEq
{
    fn eq(&self, other: &Matrix<T>) -> bool
    {
        // Using fold, though idiomatic, has average case of O(n)
        // Here we choose to use the standard for loop since avg. is O(log n)
        for (val1, val2) in self.my_dat.iter().zip(other.my_dat.iter())
        {
            if *val1 != *val2
            {
                return false;
            }
        }

        return true;
    }
}
