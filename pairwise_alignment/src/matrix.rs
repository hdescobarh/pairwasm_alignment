//! Module for matrix data structure

use std::ops::{Index, IndexMut};

pub struct Matrix<T> {
    rows: usize,
    cols: usize,
    container: Vec<T>,
}

impl<T: std::clone::Clone> Matrix<T> {
    pub fn empty(rows: usize, cols: usize) -> Self {
        let container: Vec<T> = Vec::with_capacity(rows * cols);
        Self {
            rows,
            cols,
            container,
        }
    }

    pub fn full(value: T, rows: usize, cols: usize) -> Self {
        let container: Vec<T> =
            Vec::from_iter(std::iter::repeat(value).take(rows * cols));
        Self {
            rows,
            cols,
            container,
        }
    }

    pub fn get(&self, row: usize, col: usize) -> Result<&T, MatError> {
        let index = self.map_2dim_to_1dim_index(row, col)?;
        if index >= self.container.len() {
            return Err(MatError::new(ErrorKind::EmptyAtIndex((
                [row, col],
                [self.rows, self.cols],
            ))));
        }
        Ok(&self[[row, col]])
    }

    pub fn get_mut(&mut self, row: usize, col: usize) -> Result<&mut T, MatError> {
        let index = self.map_2dim_to_1dim_index(row, col)?;
        if index >= self.container.len() {
            return Err(MatError::new(ErrorKind::EmptyAtIndex((
                [row, col],
                [self.rows, self.cols],
            ))));
        }
        Ok(&mut self[[row, col]])
    }

    fn map_2dim_to_1dim_index(&self, row: usize, col: usize) -> Result<usize, MatError> {
        if (row >= self.rows) || (col >= self.cols) {
            return Err(MatError::new(ErrorKind::OutOfDimension((
                [row, col],
                [self.rows, self.cols],
            ))));
        }
        Ok(row * self.cols + col)
    }
}

impl<T: std::clone::Clone> Index<[usize; 2]> for Matrix<T> {
    type Output = T;
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let i = match self.map_2dim_to_1dim_index(index[0], index[1]) {
            Ok(i) => i,
            Err(e) => panic!("{:?}", e),
        };
        &self.container[i]
    }
}

impl<T: std::clone::Clone> IndexMut<[usize; 2]> for Matrix<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        let i = match self.map_2dim_to_1dim_index(index[0], index[1]) {
            Ok(i) => i,
            Err(e) => panic!("{:?}", e),
        };
        &mut self.container[i]
    }
}

#[derive(Debug)]
pub struct MatError {
    kind: ErrorKind,
    message: String,
}

#[non_exhaustive]
#[derive(Debug, PartialEq)]
pub enum ErrorKind {
    /// (AttemptedIndex, ActualDimension)
    OutOfDimension(([usize; 2], [usize; 2])),
    EmptyAtIndex(([usize; 2], [usize; 2])),
}

impl MatError {
    fn new(kind: ErrorKind) -> Self {
        let message: String = match kind {
            ErrorKind::OutOfDimension((attempted, actual)) => format!(
                "The index '{attempted:?}' is out of the Matrix bounds.\
                 The Matrix dimension is '{actual:?}'."
            ),
            ErrorKind::EmptyAtIndex((attempted, actual)) => format!(
                "The index '{attempted:?}' belongs to Matrix bounds ('{actual:?}').\
                 but the location is empty."
            ),
        };

        Self { kind, message }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn creates_empty_matrix() {
        let empty_array: [u8; 0] = [];
        for (rows, cols) in [(0, 0), (3, 2), (5, 5)] {
            let matrix: Matrix<u8> = Matrix::empty(rows, cols);
            assert_eq!(rows * cols, matrix.container.capacity());
            assert_eq!(&empty_array, matrix.container.as_slice());
        }
    }

    #[test]
    fn creates_full_matrix() {
        // Zero dimension
        let matrix: Matrix<&str> = Matrix::full("🦀", 2, 0);
        assert_eq!(0, matrix.container.len());
        let empty_array: [&str; 0] = [];
        assert_eq!(&empty_array, matrix.container.as_slice());

        // non-zero dim

        let matrix: Matrix<&str> = Matrix::full("🦀", 3, 3);
        assert_eq!(9, matrix.container.len());
        assert_eq!(
            &["🦀", "🦀", "🦀", "🦀", "🦀", "🦀", "🦀", "🦀", "🦀"],
            matrix.container.as_slice()
        );
    }

    #[test]
    fn failed_get() {
        let mut matrix: Matrix<&str> = Matrix::full("🐢", 3, 3);
        // out of bounds by row
        assert!(matrix
            .get(3, 0)
            .is_err_and(|e| e.kind == ErrorKind::OutOfDimension(([3, 0], [3, 3]))));

        // out of bounds by col
        assert!(matrix
            .get(0, 3)
            .is_err_and(|e| e.kind == ErrorKind::OutOfDimension(([0, 3], [3, 3]))));

        // out of bounds because empty
        matrix.container.pop();
        assert!(matrix
            .get(2, 2)
            .is_err_and(|e| e.kind == ErrorKind::EmptyAtIndex(([2, 2], [3, 3]))));
    }

    #[test]
    fn successful_get() {
        let mut matrix: Matrix<&str> = Matrix::full("🐢", 3, 4);
        *matrix.get_mut(2, 0).unwrap() = "🦀";
        *matrix.get_mut(1, 2).unwrap() = "🦕";
        assert_eq!(
            &[
                "🐢", "🐢", "🐢", "🐢", "🐢", "🐢", "🦕", "🐢", "🦀", "🐢", "🐢", "🐢",
            ],
            matrix.container.as_slice()
        );
        assert!(matrix.get(0, 1).is_ok_and(|v| *v == "🐢"));
        assert!(matrix.get(2, 0).is_ok_and(|v| *v == "🦀"));
        assert!(matrix.get(1, 2).is_ok_and(|v| *v == "🦕"));
    }
}
