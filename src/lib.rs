pub mod aligner;
pub mod bioseq;
pub mod formatter;
pub mod matrix;
pub mod scoring_schema;
mod utils;

#[cfg(test)]
pub mod tests;

use aligner::{utils::AlignmentSequence, Aligner, AlignerKind};
use bioseq::{Aac, HasSequence};
use scoring_schema::{aminoacid_schema::AaScoringKind, gap_penalty::PenaltyKind};
use utils::AlignmentUnit;
use wasm_bindgen::prelude::*;

use crate::bioseq::Protein;

#[wasm_bindgen]
pub fn do_protein_alignment(
    string_1: &str,
    string_2: &str,
    open_cost: f32,
    extend_cost: f32,
    substitution_matrix: u8,
    algorithm: u8,
) -> String {
    let sequence_1 = Protein::new(string_1).unwrap();
    let sequence_2 = Protein::new(string_2).unwrap();

    let penalty_kind = PenaltyKind::Affine(open_cost, extend_cost);
    let score_kind = match substitution_matrix {
        b'\x01' => AaScoringKind::Blosum45,
        _ => panic!("Invalid option"),
    };

    let aligner_kind = match algorithm {
        b'\x01' => AlignerKind::SmithWaterman,
        b'\x02' => AlignerKind::NeedlemanWunsch,
        _ => panic!("Invalid option"),
    };

    let session: AlignerSession<Aac> = AlignerSession::new(
        sequence_1,
        sequence_2,
        aligner_kind,
        score_kind,
        penalty_kind,
    );

    session
        .run()
        .into_iter()
        .fold(String::new(), |acc, e| acc + &format!("{}", e))
}

pub struct AlignerSession<A>
where
    A: AlignmentUnit,
{
    aligner: Option<Box<dyn Aligner<A>>>,
}

impl<A: AlignmentUnit> AlignerSession<A> {
    fn run(self) -> Vec<AlignmentSequence<A>> {
        match self.aligner {
            Some(mut aligner) => aligner.run(),
            None => panic!("Missing aligner."),
        }
    }
}

impl AlignerSession<Aac> {
    pub fn new(
        sequence_1: impl HasSequence<Aac> + 'static,
        sequence_2: impl HasSequence<Aac> + 'static,
        aligner_kind: AlignerKind,
        score_kind: AaScoringKind,
        penalty_kind: PenaltyKind,
    ) -> Self {
        let aligner = aligner::aminoacid_align_builder(
            aligner_kind,
            sequence_1,
            sequence_2,
            score_kind,
            penalty_kind,
        );

        Self {
            aligner: Some(aligner),
        }
    }
}
