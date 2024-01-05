//! Alignment algorithms

use crate::{
    bioseq::{Aac, HasSequence},
    scoring_schema::{aminoacid_schema::AaScoringKind, gap_penalty::PenaltyKind},
    utils::AlignmentUnit,
};

use self::{
    global_alignment::NeedlemanWunsch, local_alignment::SmithWaterman,
    utils::AlignmentSequence,
};

mod global_alignment;
mod local_alignment;
pub mod utils;

/// Flag for alignment algorithm implementations
pub trait Aligner<A: AlignmentUnit> {
    fn run(&mut self) -> Vec<AlignmentSequence<A>>;
}

pub enum AlignerKind {
    NeedlemanWunsch,
    SmithWaterman,
}

/// Aligner constructor
pub fn aminoacid_align_builder(
    kind: AlignerKind,
    sequence_1: impl HasSequence<Aac> + 'static,
    sequence_2: impl HasSequence<Aac> + 'static,
    score_kind: AaScoringKind,
    penalty_kind: PenaltyKind,
) -> Box<dyn Aligner<Aac>> {
    match kind {
        AlignerKind::NeedlemanWunsch => Box::new(NeedlemanWunsch::new(
            sequence_1,
            sequence_2,
            score_kind,
            penalty_kind,
        )),
        AlignerKind::SmithWaterman => Box::new(SmithWaterman::new(
            sequence_1,
            sequence_2,
            score_kind,
            penalty_kind,
        )),
    }
}
