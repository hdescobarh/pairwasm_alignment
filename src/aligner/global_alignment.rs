// The original Needleman-Wunsch uses a linear gap penalty
use crate::bioseq::{Aac, HasSequence};
use crate::matrix::Matrix;
use crate::scoring_schema::aminoacid_schema::{self, AaScoringKind};
use crate::scoring_schema::gap_penalty::PenaltyKind;
use crate::scoring_schema::AaScoringSchema;
use crate::{scoring_schema::ScoringSchema, utils::AlignmentUnit};

use super::utils::BacktrackChoice;

pub struct NeedlemanWunsch<A>
where
    A: AlignmentUnit,
{
    sequence_s: Box<dyn HasSequence<A>>,
    sequence_t: Box<dyn HasSequence<A>>,
    scoring_schema: Box<dyn ScoringSchema<A>>,
    scores_matrix: Matrix<f32>,
    backtracking_matrix: Matrix<super::utils::BacktrackChoice>,
}

impl NeedlemanWunsch<Aac> {
    // S -> row sequence, T -> col sequence
    pub fn new(
        sequence_s: Box<dyn HasSequence<Aac>>,
        sequence_t: Box<dyn HasSequence<Aac>>,
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
        let rows = 1 + sequence_s.seq().len();
        let cols = 1 + sequence_t.seq().len();

        Self {
            sequence_s,
            sequence_t,
            scoring_schema: scoring_schema as Box<dyn ScoringSchema<Aac>>,
            scores_matrix: Matrix::full(0.0, rows, cols),
            backtracking_matrix: Matrix::full(BacktrackChoice::Empty, rows, cols),
        }
    }
}
