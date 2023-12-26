//! Scoring schemas used for creating cost and distance functions
mod aminoacid_data;
mod aminoacid_schema;
mod gap_penalty;

use crate::bioseq::Aac;
use crate::utils::AlignmentUnit;
use aminoacid_schema::*;
use gap_penalty::*;

pub trait SubstitutionSchema<A>
where
    A: AlignmentUnit,
{
    fn get_score(&self, code_1: A, code_2: A) -> i8;
}

pub struct ScoringSchema<S, P>
where
    P: Penalty,
    S: Sized,
{
    score_table: S,
    gap: P,
}

impl<S, P> ScoringSchema<S, P>
where
    S: SubstitutionSchema<Aac>,
    P: Penalty,
{
    pub fn new(score_table: S, gap: P) -> Self {
        Self { score_table, gap }
    }

    pub fn prueba() {
        let new = ScoringSchema::new(Blosum62 {}, Affine::new(0.5, 0.7));
        let x = new.score_table.get_score(Aac::A, Aac::D);
    }
}
