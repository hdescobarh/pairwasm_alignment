//! Scoring schemas used for creating cost and distance functions
mod aminoacid_data;
pub mod aminoacid_schema;
pub mod gap_penalty;

use crate::bioseq::Aac;
use crate::utils::AlignmentUnit;
use aminoacid_schema::*;
use gap_penalty::*;

/// Flag for substitution matrices
pub trait SubstitutionSchema<A>
where
    A: AlignmentUnit,
{
    fn get_score(&self, code_1: A, code_2: A) -> i8;
}

pub trait ScoringSchema<A>
where
    A: AlignmentUnit,
{
    fn get_score(&self, code_1: A, code_2: A) -> i8;

    fn get_function(&self, lenght: usize) -> f64;

    fn get_open(&self) -> f64;

    fn get_extend(&self) -> f64;
}

/// Score schema to be used with amino acid related alignments
pub struct AaScoringSchema<S, P>
where
    P: Penalty,
    S: SubstitutionSchema<Aac>,
{
    substitution: S,
    penalty: P,
}

impl<S, P> ScoringSchema<Aac> for AaScoringSchema<S, P>
where
    P: Penalty,
    S: SubstitutionSchema<Aac>,
{
    fn get_score(&self, code_1: Aac, code_2: Aac) -> i8 {
        self.substitution.get_score(code_1, code_2)
    }

    fn get_function(&self, lenght: usize) -> f64 {
        self.penalty.function(lenght)
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
