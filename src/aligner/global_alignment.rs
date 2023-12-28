// The original Needleman-Wunsch uses a linear gap penalty
use crate::bioseq::HasSequence;
use crate::matrix::Matrix;
use crate::{scoring_schema::ScoringSchema, utils::AlignmentUnit};

use super::utils::BacktrackChoice;

pub struct NeedlemanWunsch<A>
where
    A: AlignmentUnit,
{
    sequence: Box<dyn HasSequence<A>>,
    scoring_schema: Box<dyn ScoringSchema<A>>,
    scores_matrix: Matrix<f32>,
    backtracking_matrix: Matrix<super::utils::BacktrackChoice>,
}

impl<A> NeedlemanWunsch<A>
where
    A: AlignmentUnit,
{
    pub fn new(
        sequence: Box<dyn HasSequence<A>>,
        scoring_schema: Box<dyn ScoringSchema<A>>,
    ) -> Self {
        let size = sequence.seq().len();
        Self {
            sequence,
            scoring_schema,
            scores_matrix: Matrix::full(0.0, size, size),
            backtracking_matrix: Matrix::full(BacktrackChoice::Empty, size, size),
        }
    }
}
