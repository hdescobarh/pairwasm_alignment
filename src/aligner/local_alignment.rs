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

#[cfg(test)]
mod test {
    use crate::bioseq::Protein;

    use super::*;

    #[test]
    fn sw_blosum62_affine() {
        let left_string: &str =
            "MSGLRVYSTSVTGSREIKSQQSEVTRILDGKRIQYQLVDISQDNALRDEMRALAGNPKAT\
            PPQIVNGDQYCGDYELFVEAVEQNTLQEFLKLA";

        let top_string: &str =
            "MVIRVYIASSSGSTAIKKKQQDVLCFLEANKIGFEEKDIAANEENRKWMRENVPEDSRPS\
            TGYPLPPQIFNECQYRGDYDAFFEARENNAVYAFLGLTAPPGSKEAEAQANQQA";

        let sequence_left = Protein::new(left_string).unwrap();
        let sequence_top = Protein::new(top_string).unwrap();
        let mut sw = SmithWaterman::new(
            &sequence_left,
            &sequence_top,
            AaScoringKind::Blosum62,
            PenaltyKind::Affine(10.0, 1.0),
        );

        let alignments = sw.run();

        let expected_alignment = [
            [Some(Aac::L), Some(Aac::I)],
            [Some(Aac::R), Some(Aac::R)],
            [Some(Aac::V), Some(Aac::V)],
            [Some(Aac::Y), Some(Aac::Y)],
            [Some(Aac::S), Some(Aac::I)],
            [Some(Aac::T), Some(Aac::A)],
            [Some(Aac::S), Some(Aac::S)],
            [Some(Aac::V), Some(Aac::S)],
            [Some(Aac::T), Some(Aac::S)],
            [Some(Aac::G), Some(Aac::G)],
            [Some(Aac::S), Some(Aac::S)],
            [Some(Aac::R), Some(Aac::T)],
            [Some(Aac::E), Some(Aac::A)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::K), Some(Aac::K)],
            [Some(Aac::S), Some(Aac::K)],
            [Some(Aac::Q), Some(Aac::K)],
            [Some(Aac::Q), Some(Aac::Q)],
            [Some(Aac::S), Some(Aac::Q)],
            [Some(Aac::E), Some(Aac::D)],
            [Some(Aac::V), Some(Aac::V)],
            [Some(Aac::T), Some(Aac::L)],
            [Some(Aac::R), Some(Aac::C)],
            [Some(Aac::I), Some(Aac::F)],
            [Some(Aac::L), Some(Aac::L)],
            [Some(Aac::D), Some(Aac::E)],
            [Some(Aac::G), Some(Aac::A)],
            [Some(Aac::K), Some(Aac::N)],
            [Some(Aac::R), Some(Aac::K)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::Q), Some(Aac::G)],
            [Some(Aac::Y), Some(Aac::F)],
            [Some(Aac::Q), Some(Aac::E)],
            [Some(Aac::L), Some(Aac::E)],
            [Some(Aac::V), Some(Aac::K)],
            [Some(Aac::D), Some(Aac::D)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::S), Some(Aac::A)],
            [Some(Aac::Q), Some(Aac::A)],
            [Some(Aac::D), Some(Aac::N)],
            [Some(Aac::N), Some(Aac::E)],
            [Some(Aac::A), Some(Aac::E)],
            [Some(Aac::L), Some(Aac::N)],
            [Some(Aac::R), Some(Aac::R)],
            [None, Some(Aac::K)],
            [None, Some(Aac::W)],
            [None, Some(Aac::M)],
            [None, Some(Aac::R)],
            [None, Some(Aac::E)],
            [None, Some(Aac::N)],
            [None, Some(Aac::V)],
            [None, Some(Aac::P)],
            [Some(Aac::D), Some(Aac::E)],
            [Some(Aac::E), Some(Aac::D)],
            [Some(Aac::M), Some(Aac::S)],
            [Some(Aac::R), Some(Aac::R)],
            [Some(Aac::A), Some(Aac::P)],
            [Some(Aac::L), Some(Aac::S)],
            [Some(Aac::A), Some(Aac::T)],
            [Some(Aac::G), Some(Aac::G)],
            [Some(Aac::N), Some(Aac::Y)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::K), None],
            [Some(Aac::A), None],
            [Some(Aac::T), Some(Aac::L)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::Q), Some(Aac::Q)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::V), Some(Aac::F)],
            [Some(Aac::N), Some(Aac::N)],
            [Some(Aac::G), Some(Aac::E)],
            [Some(Aac::D), Some(Aac::C)],
            [Some(Aac::Q), Some(Aac::Q)],
            [Some(Aac::Y), Some(Aac::Y)],
            [Some(Aac::C), Some(Aac::R)],
            [Some(Aac::G), Some(Aac::G)],
            [Some(Aac::D), Some(Aac::D)],
            [Some(Aac::Y), Some(Aac::Y)],
            [Some(Aac::E), Some(Aac::D)],
            [Some(Aac::L), Some(Aac::A)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::V), Some(Aac::F)],
            [Some(Aac::E), Some(Aac::E)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::V), Some(Aac::R)],
            [Some(Aac::E), Some(Aac::E)],
            [Some(Aac::Q), Some(Aac::N)],
            [Some(Aac::N), Some(Aac::N)],
            [Some(Aac::T), Some(Aac::A)],
            [Some(Aac::L), Some(Aac::V)],
            [Some(Aac::Q), Some(Aac::Y)],
            [Some(Aac::E), Some(Aac::A)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::L), Some(Aac::L)],
            [Some(Aac::K), Some(Aac::G)],
            [Some(Aac::L), Some(Aac::L)],
        ];

        let actual_alignment = alignments[0].read();

        for p in 0..expected_alignment.len() {
            assert!(
                expected_alignment[p][0] == actual_alignment[p][0]
                    && expected_alignment[p][1] == actual_alignment[p][1],
                "Error at position {}. Expected [{:?}. Got [{:?}]]",
                p,
                expected_alignment[p],
                actual_alignment[p]
            )
        }
    }
}
