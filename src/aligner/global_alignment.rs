use super::Aligner;
// The original Needleman-Wunsch uses a linear gap penalty
use super::utils::{AffineTransversalOrder, AlignmentSequence, BackTrack};
use crate::bioseq::{Aac, HasSequence};
use crate::matrix::Matrix;
use crate::scoring_schema::aminoacid_schema::AaScoringKind;
use crate::scoring_schema::gap_penalty::PenaltyKind;
use crate::scoring_schema::AaScoringSchema;
use crate::{scoring_schema::ScoringSchema, utils::AlignmentUnit};

pub struct NeedlemanWunsch<A>
where
    A: AlignmentUnit,
{
    sequence_left: Box<dyn HasSequence<A>>,
    sequence_top: Box<dyn HasSequence<A>>,
    scoring_schema: Box<dyn ScoringSchema<A>>,
    matrix: Matrix<BackTrack>,
}

impl NeedlemanWunsch<Aac> {
    // S -> row sequence, T -> col sequence
    pub fn new(
        sequence_left: impl HasSequence<Aac> + 'static,
        sequence_top: impl HasSequence<Aac> + 'static,
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
            sequence_left: Box::new(sequence_left),
            sequence_top: Box::new(sequence_top),
            scoring_schema: scoring_schema as Box<dyn ScoringSchema<Aac>>,
            matrix: Matrix::full(BackTrack::Empty, rows, cols),
        }
    }
}

impl<A> NeedlemanWunsch<A>
where
    A: AlignmentUnit,
{
    fn run(&mut self) -> Vec<AlignmentSequence<A>> {
        self.initialize();
        self.solve_subproblems();
        let [row_dim, col_dim] = self.matrix.dim();
        let [init_row, init_col] = [row_dim - 1, col_dim - 1];
        let all_paths =
            BackTrack::backtracking(&self.matrix, init_row, init_col, f32::NEG_INFINITY);
        //let mut alignments: Vec<AlignmentSequence<A>> =
        //    Vec::with_capacity(all_paths.len());

        let longest_path = all_paths
            .into_iter()
            .reduce(|acc, e| if acc.len() > e.len() { acc } else { e })
            .unwrap();

        //for backtrack_path in all_paths {
        //    let new_alignment = AlignmentSequence::new(
        //        backtrack_path,
        //        self.sequence_left.as_ref(),
        //        self.sequence_top.as_ref(),
        //    );
        //    alignments.push(new_alignment);

        let alignments: Vec<AlignmentSequence<A>> = vec![AlignmentSequence::new(
            longest_path,
            self.sequence_left.as_ref(),
            self.sequence_top.as_ref(),
        )];
        alignments
    }
    fn initialize(&mut self) {
        self.matrix[[0, 0]] = BackTrack::D(0.0);
        let [rows, cols] = self.matrix.dim();

        for i in 1..rows {
            self.matrix[[i, 0]] = BackTrack::T(-self.scoring_schema.get_function(i));
        }

        for j in 1..cols {
            self.matrix[[0, j]] = BackTrack::L(-self.scoring_schema.get_function(j));
        }
    }

    fn solve_subproblems(&mut self) {
        let [rows, cols] = self.matrix.dim();
        for i in 1..rows {
            for j in 1..cols {
                let diagonal = Self::diagonal_score(
                    self.sequence_left.as_ref(),
                    self.sequence_top.as_ref(),
                    &self.scoring_schema,
                    &self.matrix,
                    i,
                    j,
                );
                let top = Self::top_score(&self.scoring_schema, &self.matrix, i, j);
                let left = Self::left_score(&self.scoring_schema, &self.matrix, i, j);
                self.matrix[[i, j]] = BackTrack::make_backtrack(top, diagonal, left).0;
            }
        }
    }
}

impl<A> AffineTransversalOrder<A> for NeedlemanWunsch<A> where A: AlignmentUnit {}

impl<A> Aligner<A> for NeedlemanWunsch<A>
where
    A: AlignmentUnit,
{
    fn run(&mut self) -> Vec<AlignmentSequence<A>> {
        self.run()
    }
}

#[cfg(test)]
mod test {
    use crate::{
        bioseq::{Aac, Protein},
        scoring_schema::{aminoacid_schema::AaScoringKind, gap_penalty::PenaltyKind},
    };

    use super::NeedlemanWunsch;

    #[test]
    fn nw_blossum62_affine() {
        let left_string: &str =
            "ATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAI\
        YNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA";

        let top_string: &str =
            "ETTQRAEREVTRMVIIMFLICWLPYAGVAWYIFTHQGSEFGPVFMTLPAFFAKTSAV\
        YNPCIYICMNKQFRHCMITTLCCGKNPFEEEEGASTTASKTEASSVSSSSVSPA";

        let sequence_left = Protein::new(left_string).unwrap();
        let sequence_top = Protein::new(top_string).unwrap();
        let mut nw = NeedlemanWunsch::new(
            sequence_left,
            sequence_top,
            AaScoringKind::Blosum62,
            PenaltyKind::Affine(10.0, 1.0),
        );

        let alignments = nw.run();

        assert_eq!(1, alignments.len());

        let expected_alignment: [[Option<Aac>; 2]; 114] = [
            [Some(Aac::A), Some(Aac::E)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::Q), Some(Aac::Q)],
            [Some(Aac::K), Some(Aac::R)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::E), Some(Aac::E)],
            [Some(Aac::K), Some(Aac::R)],
            [Some(Aac::E), Some(Aac::E)],
            [Some(Aac::V), Some(Aac::V)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::R), Some(Aac::R)],
            [Some(Aac::M), Some(Aac::M)],
            [Some(Aac::V), Some(Aac::V)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::M), Some(Aac::M)],
            [Some(Aac::V), None],
            [Some(Aac::I), None],
            [Some(Aac::A), None],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::L), Some(Aac::L)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::C), Some(Aac::C)],
            [Some(Aac::W), Some(Aac::W)],
            [Some(Aac::V), Some(Aac::L)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::Y), Some(Aac::Y)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::S), Some(Aac::G)],
            [Some(Aac::V), Some(Aac::V)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::F), Some(Aac::W)],
            [Some(Aac::Y), Some(Aac::Y)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::H), Some(Aac::H)],
            [Some(Aac::Q), Some(Aac::Q)],
            [Some(Aac::G), Some(Aac::G)],
            [Some(Aac::S), Some(Aac::S)],
            [Some(Aac::N), Some(Aac::E)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::G), Some(Aac::G)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::I), Some(Aac::V)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::M), Some(Aac::M)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::I), Some(Aac::L)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::K), Some(Aac::K)],
            [Some(Aac::S), Some(Aac::T)],
            [Some(Aac::A), Some(Aac::S)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::I), Some(Aac::V)],
            [Some(Aac::Y), Some(Aac::Y)],
            [Some(Aac::N), Some(Aac::N)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::V), Some(Aac::C)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::Y), Some(Aac::Y)],
            [Some(Aac::I), Some(Aac::I)],
            [Some(Aac::M), Some(Aac::C)],
            [Some(Aac::M), Some(Aac::M)],
            [Some(Aac::N), Some(Aac::N)],
            [Some(Aac::K), Some(Aac::K)],
            [Some(Aac::Q), Some(Aac::Q)],
            [Some(Aac::F), Some(Aac::F)],
            [Some(Aac::R), Some(Aac::R)],
            [Some(Aac::N), Some(Aac::H)],
            [Some(Aac::C), Some(Aac::C)],
            [Some(Aac::M), Some(Aac::M)],
            [Some(Aac::L), Some(Aac::I)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::I), Some(Aac::L)],
            [Some(Aac::C), Some(Aac::C)],
            [Some(Aac::C), Some(Aac::C)],
            [Some(Aac::G), Some(Aac::G)],
            [Some(Aac::K), Some(Aac::K)],
            [Some(Aac::N), Some(Aac::N)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::L), Some(Aac::F)],
            [Some(Aac::G), Some(Aac::E)],
            [Some(Aac::D), Some(Aac::E)],
            [Some(Aac::D), Some(Aac::E)],
            [Some(Aac::E), Some(Aac::E)],
            [None, Some(Aac::G)],
            [Some(Aac::A), Some(Aac::A)],
            [Some(Aac::S), Some(Aac::S)],
            [Some(Aac::A), Some(Aac::T)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::V), Some(Aac::A)],
            [Some(Aac::S), Some(Aac::S)],
            [Some(Aac::K), Some(Aac::K)],
            [Some(Aac::T), Some(Aac::T)],
            [Some(Aac::E), Some(Aac::E)],
            [None, Some(Aac::A)],
            [None, Some(Aac::S)],
            [None, Some(Aac::S)],
            [None, Some(Aac::V)],
            [None, Some(Aac::S)],
            [Some(Aac::T), Some(Aac::S)],
            [Some(Aac::S), Some(Aac::S)],
            [Some(Aac::Q), Some(Aac::S)],
            [Some(Aac::V), Some(Aac::V)],
            [Some(Aac::A), Some(Aac::S)],
            [Some(Aac::P), Some(Aac::P)],
            [Some(Aac::A), Some(Aac::A)],
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
