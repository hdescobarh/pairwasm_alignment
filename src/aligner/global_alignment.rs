// The original Needleman-Wunsch uses a linear gap penalty
use super::utils::{AlignmentSequence, BackTrack};
use crate::bioseq::{Aac, HasSequence};
use crate::matrix::Matrix;
use crate::scoring_schema::aminoacid_schema::AaScoringKind;
use crate::scoring_schema::gap_penalty::PenaltyKind;
use crate::scoring_schema::AaScoringSchema;
use crate::{scoring_schema::ScoringSchema, utils::AlignmentUnit};

pub struct NeedlemanWunsch<'a, A>
where
    A: AlignmentUnit,
{
    sequence_left: &'a dyn HasSequence<A>,
    sequence_top: &'a dyn HasSequence<A>,
    scoring_schema: Box<dyn ScoringSchema<A>>,
    matrix: Matrix<BackTrack>,
}

impl<'a> NeedlemanWunsch<'a, Aac> {
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
        }
    }
}

impl<'a, A> NeedlemanWunsch<'a, A>
where
    A: AlignmentUnit,
{
    pub fn run(&mut self) -> Vec<AlignmentSequence<A>> {
        self.initialize();
        self.solve_subproblems();
        let [row_dim, col_dim] = self.matrix.dim();
        let [init_row, init_col] = [row_dim - 1, col_dim - 1];
        let all_paths = BackTrack::backtracking(&self.matrix, init_row, init_col);
        let mut alignments: Vec<AlignmentSequence<A>> =
            Vec::with_capacity(all_paths.len());
        for backtrack_path in all_paths {
            let new_alignment = AlignmentSequence::new(
                backtrack_path,
                self.sequence_left,
                self.sequence_top,
            );
            alignments.push(new_alignment);
        }
        alignments
    }
    fn initialize(&mut self) {
        self.matrix[[0, 0]] = BackTrack::D(0.0);
        let [rows, cols] = self.matrix.dim();

        for i in 1..rows {
            self.matrix[[i, 0]] = BackTrack::T(self.scoring_schema.get_function(i));
        }

        for j in 1..cols {
            self.matrix[[0, j]] = BackTrack::L(self.scoring_schema.get_function(j));
        }
    }

    fn solve_subproblems(&mut self) {
        let [rows, cols] = self.matrix.dim();
        for i in 1..rows {
            for j in 1..cols {
                let diagonal = self.diagonal_score(i, j);
                let top = self.top_score(i, j);
                let left = self.left_score(i, j);
                self.matrix[[i, j]] = BackTrack::make_backtrack(top, diagonal, left);
            }
        }
    }

    fn read_sequences(&self, i: usize, j: usize) -> [A; 2] {
        let left_alignable: A = self.sequence_left.seq()[i - 1];
        let top_alignable: A = self.sequence_top.seq()[j - 1];
        [left_alignable, top_alignable]
    }

    fn diagonal_score(&self, i: usize, j: usize) -> f32 {
        let [left_alignable, top_alignable] = self.read_sequences(i, j);
        let score_ij = self.scoring_schema.get_score(left_alignable, top_alignable);
        let value = match self.matrix[[i - 1, j - 1]] {
            BackTrack::T(v) => v,
            BackTrack::D(v) => v,
            BackTrack::L(v) => v,
            BackTrack::DT(v) => v,
            BackTrack::DL(v) => v,
            BackTrack::TL(v) => v,
            BackTrack::All(v) => v,
            BackTrack::Empty => {
                panic!("This must be unreachable.Check the transversal order.")
            }
        };
        value + score_ij as f32
    }

    fn top_score(&self, i: usize, j: usize) -> f32 {
        // i-1, j
        match self.matrix[[i - 1, j]] {
            // top_gap + top_gap is an extension
            BackTrack::T(v) => v - self.scoring_schema.get_extend(),
            // not(top_gap) + top_gap is and opening
            BackTrack::D(v) => v - self.scoring_schema.get_function(1),
            BackTrack::L(v) => v - self.scoring_schema.get_function(1),
            BackTrack::DL(v) => v - self.scoring_schema.get_function(1),
            // Max(v - extend_existing_gap, v - add_new_gap) = v - extend_gap
            // because extend_existing_gap <= add_new_gap
            BackTrack::DT(v) => v - self.scoring_schema.get_extend(),
            BackTrack::TL(v) => v - self.scoring_schema.get_extend(),
            BackTrack::All(v) => v - self.scoring_schema.get_extend(),
            BackTrack::Empty => {
                panic!("This must be unreachable.Check the transversal order.")
            }
        }
    }

    fn left_score(&self, i: usize, j: usize) -> f32 {
        // i, j-1
        match self.matrix[[i, j - 1]] {
            // left_gap + left_gap is and extension
            BackTrack::L(v) => v - self.scoring_schema.get_extend(),
            // not(left_gap) + left_gap is a new gap
            BackTrack::T(v) => v - self.scoring_schema.get_function(1),
            BackTrack::D(v) => v - self.scoring_schema.get_function(1),
            BackTrack::DT(v) => v - self.scoring_schema.get_function(1),
            // Max(v - extend_existing_gap, v - add_new_gap) = v - extend_gap
            // because extend_existing_gap <= add_new_gap
            BackTrack::DL(v) => v - self.scoring_schema.get_extend(),
            BackTrack::TL(v) => v - self.scoring_schema.get_extend(),
            BackTrack::All(v) => v - self.scoring_schema.get_extend(),
            BackTrack::Empty => {
                panic!("This must be unreachable.Check the transversal order.")
            }
        }
    }
}
