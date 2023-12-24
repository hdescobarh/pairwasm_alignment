//! Gap penalty schemas
// A gap is a spaces chain, the lenght is the number of spaces. So,
// the minimum possible lenght of a gap is 1, then the lenght >= 1

// Flag for gap penalty schemas
pub trait Penalty {
    // The gap penalty is a map $(lenght) \mapto \mathbb(R)$.
    fn gap_penalty(&self, lenght: usize) -> f64;
}

// Implements affine gap model.
pub struct Affine {
    gap_open: f64,
    gap_extend: f64,
}

impl Affine {
    pub fn new(gap_open: f32, gap_extend: f32) -> Self {
        Self {
            gap_open: gap_open as f64,
            gap_extend: gap_extend as f64,
        }
    }

    pub fn get_open(&self) -> f64 {
        self.gap_open
    }

    pub fn get_extend(&self) -> f64 {
        self.gap_extend
    }
}

impl Penalty for Affine {
    fn gap_penalty(&self, lenght: usize) -> f64 {
        self.gap_open + (self.gap_extend * lenght as f64)
    }
}

pub struct Constant {
    value: f64,
}

impl Constant {
    pub fn new(value: f32) -> Self {
        Self {
            value: value as f64,
        }
    }
}

impl Penalty for Constant {
    fn gap_penalty(&self, _lenght: usize) -> f64 {
        self.value
    }
}

#[cfg(test)]
mod test {
    use super::*;
}
