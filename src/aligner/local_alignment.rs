//! Algorithms for local alignment

use super::utils::{AffineTransversalOrder, AlignmentSequence, BackTrack};
use crate::bioseq::{Aac, HasSequence};
use crate::matrix::Matrix;
use crate::scoring_schema::aminoacid_schema::AaScoringKind;
use crate::scoring_schema::gap_penalty::PenaltyKind;
use crate::scoring_schema::AaScoringSchema;
use crate::{scoring_schema::ScoringSchema, utils::AlignmentUnit};

/// Smith Waterman original algorithm. Returns the longest and best local alignment.
pub struct SmithWaterman<'a, A>
where
    A: AlignmentUnit,
{
    sequence_left: &'a dyn HasSequence<A>,
    sequence_top: &'a dyn HasSequence<A>,
    scoring_schema: Box<dyn ScoringSchema<A>>,
    matrix: Matrix<BackTrack>,
    /// The highest found score
    global_maximum: f32,
    /// Indices whose score is the global maximum
    maximum_indices: Vec<[usize; 2]>,
}

impl<'a, A> SmithWaterman<'a, A>
where
    A: AlignmentUnit,
{
    pub fn run(&mut self) -> Vec<AlignmentSequence<A>> {
        self.initialize();
        self.solve_subproblems();

        let mut all_paths: Vec<Vec<[usize; 2]>> = Vec::new();

        for [init_row, init_col] in &self.maximum_indices {
            let mut path =
                BackTrack::backtracking(&self.matrix, *init_row, *init_col, 0.0);
            all_paths.append(&mut path);
        }

        let longest_path = all_paths
            .into_iter()
            .reduce(|acc, e| if acc.len() > e.len() { acc } else { e })
            .unwrap();

        let alignments: Vec<AlignmentSequence<A>> = vec![AlignmentSequence::new(
            longest_path,
            self.sequence_left,
            self.sequence_top,
        )];
        alignments
    }

    fn initialize(&mut self) {
        self.matrix[[0, 0]] = BackTrack::D(0.0);
        let [rows, cols] = self.matrix.dim();

        for i in 1..rows {
            self.matrix[[i, 0]] = BackTrack::T(0.0);
        }

        for j in 1..cols {
            self.matrix[[0, j]] = BackTrack::L(0.0);
        }
    }

    fn solve_subproblems(&mut self) {
        let [rows, cols] = self.matrix.dim();
        for i in 1..rows {
            for j in 1..cols {
                let diagonal = Self::diagonal_score(
                    self.sequence_left,
                    self.sequence_top,
                    &self.scoring_schema,
                    &self.matrix,
                    i,
                    j,
                )
                .max(0.0);
                let top =
                    Self::top_score(&self.scoring_schema, &self.matrix, i, j).max(0.0);
                let left =
                    Self::left_score(&self.scoring_schema, &self.matrix, i, j).max(0.0);

                let (backtrack, current_maximum) =
                    BackTrack::make_backtrack(top, diagonal, left);
                self.update_maximum_entries(current_maximum, i, j);
                self.matrix[[i, j]] = backtrack;
            }
        }
    }

    fn update_maximum_entries(&mut self, current_maximum: f32, i: usize, j: usize) {
        // This comparisons may need to be improved because similarity is calculated with subtractions
        // and they are ill-conditioned.
        if current_maximum > self.global_maximum {
            self.global_maximum = current_maximum;
            self.maximum_indices = vec![[i, j]]
        } else if current_maximum == self.global_maximum {
            self.maximum_indices.push([i, j])
        }
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
            global_maximum: f32::NEG_INFINITY,
            maximum_indices: Vec::new(),
        }
    }
}
