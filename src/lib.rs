pub mod aligner;
pub mod bioseq;
pub mod formatter;
pub mod matrix;
pub mod scoring_schema;
mod utils;

#[cfg(test)]
pub mod tests;

use aligner::AlignerKind;
use scoring_schema::{aminoacid_schema::AaScoringKind, gap_penalty::PenaltyKind};
use utils::set_panic_hook;
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
) -> Result<String, JsError> {
    // set panic_hook
    set_panic_hook();

    // set and run alignment
    let sequence_1 = Protein::new(string_1)?;
    let sequence_2 = Protein::new(string_2)?;

    let penalty_kind = PenaltyKind::Affine(open_cost, extend_cost);
    let score_kind = match substitution_matrix {
        b'\x01' => AaScoringKind::Blosum45,
        b'\x02' => AaScoringKind::Blosum62,
        b'\x03' => AaScoringKind::Pam160,
        _ => panic!("Invalid option"),
    };

    let aligner_kind = match algorithm {
        b'\x01' => AlignerKind::SmithWaterman,
        b'\x02' => AlignerKind::NeedlemanWunsch,
        _ => panic!("Invalid option"),
    };

    let mut aligner_instance = aligner::aminoacid_align_builder(
        aligner_kind,
        sequence_1,
        sequence_2,
        score_kind,
        penalty_kind,
    );

    Ok(aligner_instance
        .run()
        .into_iter()
        .fold(String::new(), |acc, e| acc + &format!("{}", e)))
}
