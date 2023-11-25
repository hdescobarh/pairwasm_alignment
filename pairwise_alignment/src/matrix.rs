//! Module for matrix data structure

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
        match self.container.get(index) {
            Some(v) => Ok(v),
            None => Err(MatError::new(ErrorKind::EmptyAtIndex((
                [row, col],
                [self.rows, self.cols],
            )))),
        }
    }

    pub fn get_mut(
        &mut self,
        row: usize,
        col: usize,
    ) -> Result<&mut T, MatError> {
        let index = self.map_2dim_to_1dim_index(row, col)?;
        match self.container.get_mut(index) {
            Some(v) => Ok(v),
            None => Err(MatError::new(ErrorKind::EmptyAtIndex((
                [row, col],
                [self.rows, self.cols],
            )))),
        }
    }

    fn map_2dim_to_1dim_index(
        &self,
        row: usize,
        col: usize,
    ) -> Result<usize, MatError> {
        if (row >= self.rows) || (col >= self.cols) {
            return Err(MatError::new(ErrorKind::OutOfDimension((
                [row, col],
                [self.rows, self.cols],
            ))));
        }
        Ok(row * self.cols + col)
    }
}

#[derive(Debug)]
pub struct MatError {
    kind: ErrorKind,
    message: String,
}

#[derive(Debug)]
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
        // out of bounds by row
        // out of bounds by col
        // out of bounds because empty
    }
}
