//! Scoring schemas used for creating cost and distance functions
mod aminacid_data;
pub mod aminoacid_schema;
pub mod gap_penalty;
pub mod substitution; // delete this

use crate::utils::AlignmentUnit;
use gap_penalty::*;

pub trait SubstitutionSchema<A>
where
    A: AlignmentUnit,
{
    fn get_score(&self, code_1: A, code_2: A) -> i8;
}

pub struct SchoringSchema<S, P>
where
    P: Penalty,
    S: Sized,
{
    score_table: S,
    gap: P,
}
