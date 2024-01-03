//! Algorithms for local alignment

use super::utils::{AffineTransversalOrder, AlignmentSequence, BackTrack};
use crate::bioseq::{Aac, HasSequence};
use crate::matrix::Matrix;
use crate::scoring_schema::aminoacid_schema::AaScoringKind;
use crate::scoring_schema::gap_penalty::PenaltyKind;
use crate::scoring_schema::AaScoringSchema;
use crate::{scoring_schema::ScoringSchema, utils::AlignmentUnit};

pub struct SmithWaterman<'a, A>
where
    A: AlignmentUnit,
{
    sequence_left: &'a dyn HasSequence<A>,
    sequence_top: &'a dyn HasSequence<A>,
    scoring_schema: Box<dyn ScoringSchema<A>>,
    matrix: Matrix<BackTrack>,
    maximum_score: f32,
    starting_indices: Vec<[usize; 2]>,
}

impl<'a, A> SmithWaterman<'a, A>
where
    A: AlignmentUnit,
{
    pub fn run(&mut self) -> Vec<AlignmentSequence<A>> {
        todo!()
    }

    fn initialize(&mut self) {
        todo!()
    }

    fn solve_subproblems(&mut self) {
        todo!()
    }
}

impl<'a, A> AffineTransversalOrder<A> for SmithWaterman<'a, A> where A: AlignmentUnit {}

impl<'a> SmithWaterman<'a, Aac> {
    // S -> row sequence, T -> col sequence
    pub fn new(
        sequence_left: &'a dyn HasSequence<Aac>,
        sequence_top: &'a dyn HasSequence<Aac>,
        score_kind: AaScoringKind,
        penalty_kind: PenaltyKind,
    ) -> Self {
        // Implementation only valid for linear and affine
        #[allow(unreachable_patterns)]
        match penalty_kind {
            PenaltyKind::Affine(_, _) => (),
            PenaltyKind::Linear(_) => (),
            _ => panic!("Only allowed for Affine and Linear gap models."),
        }
        let scoring_schema = Box::new(AaScoringSchema::new(score_kind, penalty_kind));
        let rows = 1 + sequence_left.seq().len();
        let cols = 1 + sequence_top.seq().len();

        Self {
            sequence_left,
            sequence_top,
            scoring_schema: scoring_schema as Box<dyn ScoringSchema<Aac>>,
            matrix: Matrix::full(BackTrack::Empty, rows, cols),
            maximum_score: f32::NEG_INFINITY,
            starting_indices: Vec::new(),
        }
    }
}
