// Common traits

/// Transposable trait for matrices and arrays
pub trait Transposable
{
    /// Transposes a matrix or array
    fn transpose(self) -> Self;
}
