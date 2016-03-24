// Matrix type

// #![feature(step_by)]  // <- Uncomment for static analysis tools

use std::ops::{Add, Sub, Index, IndexMut, Mul};
use std::cmp::PartialEq;

/// Regular ol' matrix. No sparseness or any of that garbage, that's for later.
#[derive(Debug, Clone)]
pub struct Matrix<T> where T: Add + Sub + Copy + PartialEq
{
    my_dat: Vec<T>,
    col_ordering: bool,
    num_rows: usize,
    num_cols: usize,
}

/// Returns the dimensions of a vector of vectors as a two-tuple.
fn get_content_dims<T>(contents: &Vec<Vec<T>>) -> (usize, usize)
    where T: Add + Sub + Copy + PartialEq
{
    let num_rows = contents.len();
    if num_rows == 0
    {
        return (0, 0);
    } else {
        return (num_rows, contents[0].len());
    }
}

/// Validates whether all vectors in a vector of vectors are size-consistent.
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

/// Implementation of the Sub trait
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

/// Implementation of partial equality testing
impl<T> PartialEq for Matrix<T> where T: Add<Output=T> + Sub + Copy + PartialEq
{
    fn eq(&self, other: &Matrix<T>) -> bool
    {
        // Using fold, though idiomatic, has average case of O(n).
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

/// Primary implementation of public methods for generic matrices
impl <T> Matrix<T> where T: Add + Sub + Copy + PartialEq
{
    /// Initialization method to create a Matrix from a vector of vectors.
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

    /// Creates a new Matrix filled with a given value, with given dimensions
    pub fn new_filled(value: T, num_rows: usize, num_cols: usize) -> Matrix<T>
    {
        Matrix
        {
            my_dat: vec![value; num_rows * num_cols],
            col_ordering: true,
            num_rows: num_rows,
            num_cols: num_cols,
        }
    }

    /// Creates a new Matrix from a Vector, with given dimensions.
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
    
    /// Transposes a matrix by 'rotating' it by its diagonal.
    pub fn transpose(self) -> Matrix<T>
    {
        let mut new_vec = Vec::<T>::with_capacity(self.my_dat.len());
        let new_rows = self.num_cols;
        let new_cols = self.num_rows;

        for i in 0..new_rows
        {
            for new_indx in (i..new_vec.len()).step_by(new_rows)
            {
                new_vec.push(self.my_dat[new_indx])
            }
        }

        Matrix
        {
            my_dat: new_vec,
            col_ordering: true,
            num_rows: new_rows,
            num_cols: new_cols,
        }
    }
    
    /// Returns true if-and-only-if the matrix is square.
    pub fn is_square(&self) -> bool
    {
        self.num_rows == self.num_cols
    }

    /// Zero-out upper trigonal "corner" of the matrix with the given "zero" value.
    pub fn zero_trigonal_upper(&mut self, zero_val: T)
    {
        if !self.is_square()
        {
            panic!("Matrix must be square.");
        }
        let n = self.num_rows;
        for i in 0..n
        {
            for j in 0..n
            {
                if i > j
                {
                    self.my_dat[j * n + i] = zero_val;
                }
            }
        }
    }

    /// Zero-out lower trigonal "corner" of the matrix with the given "zero" value.
    pub fn zero_trigonal_lower(&mut self, zero_val: T)
    {
        if !self.is_square()
        {
            panic!("Matrix must be square.");
        }
        let n = self.num_rows;
        for i in 0..n
        {
            for j in 0..n
            {
                if i < j
                {
                    self.my_dat[j * n + i] = zero_val;
                }
            }
        }
    }

    /// Gets a reference to the value at the matrix's index coordinates.
    pub fn get(&self, i: usize, j: usize) -> &T
    {
        if i > self.num_rows || j > self.num_cols
        {
            panic!("Index out of bounds.");
        }

        return &self.my_dat[(i * self.num_cols) + j];
    }

    /// Sets the value at the matrix's index coordinates.
    pub fn set(&mut self, i: usize, j: usize, val: T)
    {
        if i > self.num_rows || j > self.num_cols
        {
            panic!("Indices out of bounds.");
        }

        self.my_dat[(i * self.num_cols) + j] = val;
    }

    /// Gets the dimension of the matrix as a two-tuple (num_rows, num_cols).
    pub fn get_dims(&self) -> (usize, usize)
    {
        (self.num_rows, self.num_cols)
    }

    /// Gets the underlying matrix data as a raw pointer.
    pub unsafe fn as_ptr(&self) -> *const T
    {
        self.my_dat.as_ptr()
    }

    /// Gets the underlying matrix data as a mutable pointer.
    pub unsafe fn as_mut_ptr(&mut self) -> *mut T
    {
        self.my_dat.as_mut_ptr()
    }
}

/// Allows for using a tuple as indexing, e.g. m[(1, 2)]
impl <T> Index<(usize, usize)> for Matrix<T>
    where T: Add + Sub + Copy + PartialEq
{
    type Output = T;

    fn index<'a>(&'a self, indices: (usize, usize)) -> &'a T
    {
        let (i, j) = indices;
        return self.get(i, j);
    }
}

/// Allows for using a tuple as mutable indexing, e.g. m[(1, 2)] = 4
impl <T> IndexMut<(usize, usize)> for Matrix<T>
    where T: Add + Sub + Copy + PartialEq
{
    fn index_mut<'a>(&'a mut self, indices: (usize, usize)) -> &'a mut T
    {
        let (i, j) = indices;
        if i > self.num_rows || j > self.num_cols
        {
            panic!("Indices out of bounds.");
        }
        return &mut self.my_dat[i * self.num_rows + j];
    }
}

/// Generic matrix product implementation (for non-BLAS types).
impl <T> Mul for Matrix<T>
    where T: Add<Output=T> + Sub + Copy + PartialEq + Mul<Output=T>
{
    type Output = Matrix<T>;

    fn mul(self, _rhs: Matrix<T>) -> Matrix<T>
    {
        // A is n x m, B is m x p.
        if self.num_cols != _rhs.num_rows  // m != m
        {
            panic!("Invalid dimensions for matrix multiplication.")
        }

        let n = self.num_rows;
        let m = self.num_cols;
        let p = _rhs.num_cols;

        // TODO: Rewrite as a series of folds, colletcs, etc.
        let mut toplevel: Vec<Vec<T>> = Vec::with_capacity(n);
        for i in 0..n
        {
            let mut curr_col: Vec<T> = Vec::with_capacity(p);
            for j in 0..p
            {
                // AB[(i, j)] = ...
                let sum : Option<T> = (0..m).fold(
                    None,
                    |sum_val, k| match sum_val {
                        None => Some(self[(i, k)] * _rhs[(k, j)]),
                        Some(curr_val) => Some(curr_val + (self[(i, k)] * _rhs[(k, j)]))
                    }
                    );
                
                match sum {
                    Some(sum_val) => curr_col.push(sum_val),
                    None => panic!("Error computing matrix product")
                };
            }
            toplevel.push(curr_col);
        }
        
        // Return new matrix
        return Matrix::<T>::new(toplevel);
    }
}

/// Convenience macro for creating matrices, like vec!
#[macro_export]
macro_rules! mat
{
    (
        $(
            [ $( $x:expr ),* ]
        ),*
    ) => {
        matrix::Matrix::new(vec![ $( vec![$( $x ),*]),*])
    }
}
