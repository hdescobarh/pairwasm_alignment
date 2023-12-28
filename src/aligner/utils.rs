//! common data structures and functions used for multiple align algorithms

/// Represent values for backtracking
#[derive(Clone)]
pub enum BacktrackChoice {
    /// Without assigned value. Used for initialize data.
    Empty = 0,
    /// Up
    U,
    /// Up-left diagonal
    D,
    /// Left
    L,
}
