//! common data structures and functions used for multiple align algorithms

/// Represent values for backtracking
#[derive(Clone)]
pub enum BacktrackChoice {
    /// Without assigned value. Used for initialize data.
    Empty = 0,
    /// Up-left diagonal: Match or Mismatch
    D,
    /// Left: gap in side sequence (rows)
    L,
    /// Up: gap in top sequence (cols)
    U,
}
