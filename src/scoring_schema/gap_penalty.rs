//! Gap penalty schemas.
// A gap is a spaces chain, the length is the number of spaces. So,
// the minimum possible length of a gap is 1, then the length >= 1

use super::CostType;
use super::GapPenalty;

pub const MIN_OPEN_COST: CostType = 1.0;
pub const MAX_OPEN_COST: CostType = 100.0;
pub const MIN_EXTEND_COST: CostType = 0.0;
pub const MAX_EXTEND_COST: CostType = 10.0;

pub enum PenaltyKind {
    // open_cost: CostType, extend_cost: CostType
    Affine(CostType, CostType),
    // extend_cost: CostType
    Linear(CostType),
}

/// Penalty constructor
pub fn penalty_builder(kind: PenaltyKind) -> Box<dyn GapPenalty> {
    match kind {
        PenaltyKind::Affine(open_cost, extend_cost) => {
            Box::new(Affine::new(open_cost, extend_cost))
        }
        PenaltyKind::Linear(extend_cost) => Box::new(Linear::new(extend_cost)),
    }
}

/// Implements affine gap model.
/// f(length) = open_cost + extend_cost * length, length \in Z+
pub struct Affine {
    open_cost: CostType,
    extend_cost: CostType,
}

impl Affine {
    fn new(open_cost: CostType, extend_cost: CostType) -> Self {
        check_open_cost(&open_cost);
        check_extend_cost(&extend_cost);
        Self {
            open_cost,
            extend_cost,
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
    fn new(extend_cost: CostType) -> Self {
        check_extend_cost(&extend_cost);
        Self { extend_cost }
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
    // To guard in case of a future implementation changes length type
    #[allow(clippy::absurd_extreme_comparisons)]
    if length <= 0 {
        panic!("Length must be a positive value.")
    }
}

fn check_open_cost(open_cost: &CostType) {
    if !(MIN_OPEN_COST..=MAX_OPEN_COST).contains(open_cost) {
        panic!(
            "Invalid gap open cost ({}). It must be in the closed interval [{}, {}].",
            open_cost, MIN_OPEN_COST, MAX_OPEN_COST
        )
    }
}

fn check_extend_cost(extend_cost: &CostType) {
    if !(MIN_EXTEND_COST..=MAX_EXTEND_COST).contains(extend_cost) {
        panic!(
            "Invalid gap extend cost ({}). It must be in the closed interval [{}, {}].",
            extend_cost, MIN_EXTEND_COST, MAX_EXTEND_COST
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn valid_affine() {
        let gap_model = penalty_builder(PenaltyKind::Affine(1.0, 0.5));
        assert_eq!(6.0, gap_model.function(10));
        assert_eq!(6.0, gap_model.open() + gap_model.extend() * 10.0);

        let gap_model = penalty_builder(PenaltyKind::Affine(15.0, 2.0));
        assert_eq!(21.0, gap_model.function(3));
        assert_eq!(21.0, gap_model.open() + gap_model.extend() * 3.0)
    }

    #[test]
    #[should_panic(
        expected = "Invalid gap open cost (-1). It must be in the closed interval [1, 100]."
    )]
    fn invalid_affine_param_open_negative() {
        penalty_builder(PenaltyKind::Affine(-1.0, 0.5));
    }

    #[test]
    #[should_panic(
        expected = "Invalid gap open cost (0). It must be in the closed interval [1, 100]."
    )]
    fn invalid_affine_param_open_zero() {
        penalty_builder(PenaltyKind::Affine(0.0, 0.5));
    }

    #[test]
    #[should_panic(
        expected = "Invalid gap extend cost (-0.5). It must be in the closed interval [0, 10]."
    )]
    fn invalid_affine_param_extend_negative() {
        penalty_builder(PenaltyKind::Affine(1.0, -0.5));
    }

    #[test]
    #[should_panic(expected = "Length must be a positive value.")]
    fn invalid_affine_length() {
        let gap_model = penalty_builder(PenaltyKind::Affine(1.0, 0.5));
        gap_model.function(0);
    }

    #[test]
    fn valid_linear() {
        let gap_model = penalty_builder(PenaltyKind::Linear(0.5));
        assert_eq!(4.5, gap_model.function(9));
        assert_eq!(4.5, gap_model.open() + gap_model.extend() * 9.0);

        let gap_model = penalty_builder(PenaltyKind::Linear(9.0));
        assert_eq!(27.0, gap_model.function(3));
        assert_eq!(27.0, gap_model.open() + gap_model.extend() * 3.0)
    }

    #[test]
    #[should_panic(
        expected = "Invalid gap extend cost (-0.5). It must be in the closed interval [0, 10]."
    )]
    fn invalid_linear_param_extend_negative() {
        penalty_builder(PenaltyKind::Linear(-0.5));
    }

    #[test]
    #[should_panic(expected = "Length must be a positive value.")]
    fn invalid_linear_length() {
        let gap_model = penalty_builder(PenaltyKind::Linear(7.0));
        gap_model.function(0);
    }
}
