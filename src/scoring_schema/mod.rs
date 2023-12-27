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
use aminoacid_schema::*;
use gap_penalty::*;

/// Scoring schema's similarity measure component
pub trait Similarity<A>
where
    A: AlignmentUnit,
{
    fn read_score(&self, code_1: A, code_2: A) -> i8;
}

/// Scoring schema's gap penalty component
pub trait GapPenalty {
    /// The gap penalty is a map $(length) \mapto \mathbb(R)$.
    fn function(&self, length: usize) -> f64;

    /// Get the open gap parameter. Be aware that under some gap penalty models
    /// this value can be different from calling f(1).
    /// For example, under affine model, f(1) = open_cost + extend_cost * 1.
    fn open(&self) -> f64;
    /// Get the extend gap parameter.
    fn extend(&self) -> f64;
}

/// Scoring schema getters
pub trait ScoringSchema<A>
where
    A: AlignmentUnit,
{
    fn get_score(&self, code_1: A, code_2: A) -> i8;

    fn get_function(&self, length: usize) -> f64;

    fn get_open(&self) -> f64;

    fn get_extend(&self) -> f64;
}

/// Amino acid sequence scoring schema
pub struct AaScoringSchema<S, P>
where
    P: GapPenalty,
    S: Similarity<Aac>,
{
    substitution: S,
    penalty: P,
}

impl<S, P> ScoringSchema<Aac> for AaScoringSchema<S, P>
where
    P: GapPenalty,
    S: Similarity<Aac>,
{
    fn get_score(&self, code_1: Aac, code_2: Aac) -> i8 {
        self.substitution.read_score(code_1, code_2)
    }

    fn get_function(&self, length: usize) -> f64 {
        self.penalty.function(length)
    }

    fn get_open(&self) -> f64 {
        self.penalty.open()
    }

    fn get_extend(&self) -> f64 {
        self.penalty.open()
    }
}

impl AaScoringSchema<Blosum45, Affine> {
    pub fn new(open_cost: f32, extend_cost: f32) -> Self {
        Self {
            substitution: Blosum45 {},
            penalty: Affine::new(open_cost, extend_cost),
        }
    }
}

impl AaScoringSchema<Blosum45, Linear> {
    pub fn new(extend_cost: f32) -> Self {
        Self {
            substitution: Blosum45 {},
            penalty: Linear::new(extend_cost),
        }
    }
}

impl AaScoringSchema<Blosum62, Affine> {
    pub fn new(open_cost: f32, extend_cost: f32) -> Self {
        Self {
            substitution: Blosum62 {},
            penalty: Affine::new(open_cost, extend_cost),
        }
    }
}

impl AaScoringSchema<Blosum62, Linear> {
    pub fn new(extend_cost: f32) -> Self {
        Self {
            substitution: Blosum62 {},
            penalty: Linear::new(extend_cost),
        }
    }
}

impl AaScoringSchema<Pam160, Affine> {
    pub fn new(open_cost: f32, extend_cost: f32) -> Self {
        Self {
            substitution: Pam160 {},
            penalty: Affine::new(open_cost, extend_cost),
        }
    }
}

impl AaScoringSchema<Pam160, Linear> {
    pub fn new(extend_cost: f32) -> Self {
        Self {
            substitution: Pam160 {},
            penalty: Linear::new(extend_cost),
        }
    }
}
