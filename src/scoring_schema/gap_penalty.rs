//! Gap penalty schemas.
// A gap is a spaces chain, the lenght is the number of spaces. So,
// the minimum possible lenght of a gap is 1, then the lenght >= 1

use super::GapPenalty;

// This must be a primitive float.
type CostType = f64;

/// Implements affine gap model.
/// f(lenght) = open_cost + extend_cost * lenght, lenght \in Z+
pub struct Affine {
    open_cost: CostType,
    extend_cost: CostType,
}

impl Affine {
    pub fn new(open_cost: f32, extend_cost: f32) -> Self {
        Self {
            open_cost: open_cost as CostType,
            extend_cost: extend_cost as CostType,
        }
    }
}

impl GapPenalty for Affine {
    fn function(&self, lenght: usize) -> CostType {
        self.open_cost + (self.extend_cost * lenght as CostType)
    }
    fn open(&self) -> CostType {
        self.open_cost
    }
    fn extend(&self) -> CostType {
        self.extend_cost
    }
}

/// Implements linear gap model.
/// f(lenght) = extend_cost * lenght, lenght \in Z+.
/// open_cost is a constant 0.
pub struct Linear {
    extend_cost: CostType,
}

impl Linear {
    pub fn new(extend_cost: f32) -> Self {
        Self {
            extend_cost: extend_cost as CostType,
        }
    }
}

impl GapPenalty for Linear {
    fn function(&self, lenght: usize) -> CostType {
        self.extend_cost * lenght as CostType
    }

    fn open(&self) -> CostType {
        0.0
    }

    fn extend(&self) -> CostType {
        self.extend_cost
    }
}

#[cfg(test)]
mod test {
    use super::*;
}
