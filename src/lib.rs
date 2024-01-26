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
use std::{error, fmt};
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
        _ => Err(InputError::new(InputErrorKind::ScoringMatrixNotExist))?,
    };

    let aligner_kind = match algorithm {
        b'\x01' => AlignerKind::SmithWaterman,
        b'\x02' => AlignerKind::NeedlemanWunsch,
        _ => Err(InputError::new(InputErrorKind::AlignerNotExist))?,
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

#[derive(Debug)]
/// Error type for input operations.
pub struct InputError {
    kind: InputErrorKind,
    message: String,
}

#[non_exhaustive]
#[derive(Debug, PartialEq)]
/// A list specifying general error categories of InputError.
pub enum InputErrorKind {
    AlignerNotExist,
    GapModelNotExist,
    ScoringMatrixNotExist,
}

impl InputError {
    fn new(kind: InputErrorKind) -> Self {
        let mut message: String = match kind {
            InputErrorKind::AlignerNotExist => {
                "The chosen aligner algorithm does not exist.".to_string()
            }
            InputErrorKind::GapModelNotExist => {
                "The chosen gap model does not exist.".to_string()
            }
            InputErrorKind::ScoringMatrixNotExist => {
                "The chosen scoring matrix does not exist.".to_string()
            }
        };

        message.push_str(" Please check the documentation for more information.");

        Self { kind, message }
    }
}

impl fmt::Display for InputError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({:?}) {}", self.kind, self.message)
    }
}

impl error::Error for InputError {}
