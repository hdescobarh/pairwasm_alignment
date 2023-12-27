//! Gap penalty schemas.
// A gap is a spaces chain, the length is the number of spaces. So,
// the minimum possible length of a gap is 1, then the length >= 1

use super::CostType;
use super::GapPenalty;

pub const MIN_OPEN_COST: CostType = 1.0;
pub const MAX_OPEN_COST: CostType = 100.0;
pub const MIN_EXTEND_COST: CostType = 0.0;
pub const MAX_EXTEND_COST: CostType = 10.0;

/// Implements affine gap model.
/// f(length) = open_cost + extend_cost * length, length \in Z+
pub struct Affine {
    open_cost: CostType,
    extend_cost: CostType,
}

impl Affine {
    pub fn new(open_cost: CostType, extend_cost: CostType) -> Self {
        check_open_cost(&open_cost);
        check_extend_cost(&extend_cost);
        Self {
            open_cost: open_cost,
            extend_cost: extend_cost,
        }
    }
}

impl GapPenalty for Affine {
    fn function(&self, length: usize) -> CostType {
        check_length(length);
        self.open_cost + (self.extend_cost * length as CostType)
    }
    fn open(&self) -> CostType {
        self.open_cost
    }
    fn extend(&self) -> CostType {
        self.extend_cost
    }
}

/// Implements linear gap model.
/// f(length) = extend_cost * length, length \in Z+.
/// open_cost is a constant 0.
pub struct Linear {
    extend_cost: CostType,
}

impl Linear {
    pub fn new(extend_cost: CostType) -> Self {
        check_extend_cost(&extend_cost);
        Self {
            extend_cost: extend_cost,
        }
    }
}

impl GapPenalty for Linear {
    fn function(&self, length: usize) -> CostType {
        check_length(length);
        self.extend_cost * length as CostType
    }

    fn open(&self) -> CostType {
        0.0
    }

    fn extend(&self) -> CostType {
        self.extend_cost
    }
}

fn check_length(length: usize) {
    if length == 0 {
        panic!("A length of 0 is not allowed.")
    }
}

fn check_open_cost(open_cost: &CostType) {
    if !(MIN_OPEN_COST..=MAX_OPEN_COST).contains(open_cost) {
        panic!(
            "Invalid gap open cost ({}). It must be in the closed interval [{}, {}]",
            open_cost, MIN_OPEN_COST, MAX_OPEN_COST
        )
    }
}

fn check_extend_cost(extend_cost: &CostType) {
    if !(MIN_EXTEND_COST..=MAX_EXTEND_COST).contains(extend_cost) {
        panic!(
            "Invalid gap extend cost ({}). It must be in the closed interval [{}, {}]",
            extend_cost, MIN_EXTEND_COST, MAX_EXTEND_COST
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;
}
