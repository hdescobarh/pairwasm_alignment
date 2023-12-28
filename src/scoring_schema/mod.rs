//! Scoring systems for sequence alignment.
//!
//! Each scoring schema is composed of a measure of similarity between the
//! sequence units (e.g. substitution and transition/transversion matrices), and
//! penalties for the presence of gaps.

mod aminoacid_data;
pub mod aminoacid_schema;
pub mod gap_penalty;

use crate::bioseq::Aac;
use crate::utils::AlignmentUnit;

use self::aminoacid_schema::AaScoringKind;
use self::gap_penalty::PenaltyKind;

type CostType = f32;
type SimilarityType = i8;

/// Scoring schema's similarity measure component
pub trait Similarity<A>
where
    A: AlignmentUnit,
{
    fn read_score(&self, code_1: A, code_2: A) -> SimilarityType;
}

/// Scoring schema's gap penalty component
pub trait GapPenalty {
    /// The gap penalty is a map $(length) \mapto \mathbb(R)$.
    fn function(&self, length: usize) -> CostType;

    /// Get the open gap parameter. Be aware that under some gap penalty models
    /// this value can be different from calling f(1).
    /// For example, under affine model, f(1) = open_cost + extend_cost * 1.
    fn open(&self) -> CostType;
    /// Get the extend gap parameter.
    fn extend(&self) -> CostType;
}

/// Scoring schema getters
pub trait ScoringSchema<A>
where
    A: AlignmentUnit,
{
    fn get_score(&self, code_1: A, code_2: A) -> SimilarityType;

    fn get_function(&self, length: usize) -> CostType;

    fn get_open(&self) -> CostType;

    fn get_extend(&self) -> CostType;
}

/// Amino acid sequence scoring schema

pub struct AaScoringSchema {
    substitution: Box<dyn Similarity<Aac>>,
    penalty: Box<dyn GapPenalty>,
}

impl AaScoringSchema {
    pub fn new(score_kind: AaScoringKind, penalty_kind: PenaltyKind) -> Self {
        let substitution = aminoacid_schema::similarity_builder(score_kind);
        let penalty = gap_penalty::penalty_builder(penalty_kind);
        Self {
            substitution,
            penalty,
        }
    }
}

impl ScoringSchema<Aac> for AaScoringSchema {
    fn get_score(&self, code_1: Aac, code_2: Aac) -> SimilarityType {
        self.substitution.read_score(code_1, code_2)
    }

    fn get_function(&self, length: usize) -> CostType {
        self.penalty.function(length)
    }

    fn get_open(&self) -> CostType {
        self.penalty.open()
    }

    fn get_extend(&self) -> CostType {
        self.penalty.extend()
    }
}
