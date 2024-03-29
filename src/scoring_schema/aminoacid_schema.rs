//! Amino acid scoring schemas

use super::aminoacid_data;
use super::Similarity;
use super::SimilarityType;
use crate::bioseq::Aac;

/// Represents available Amino acid scoring matrices
pub enum AaScoringKind {
    Blosum45,
    Blosum62,
    Pam160,
}

/// Similarity schema constructor
pub fn similarity_builder(kind: AaScoringKind) -> Box<dyn Similarity<Aac>> {
    match kind {
        AaScoringKind::Blosum45 => Box::new(Blosum45 {}),
        AaScoringKind::Blosum62 => Box::new(Blosum62 {}),
        AaScoringKind::Pam160 => Box::new(Pam160 {}),
    }
}

pub struct Blosum45 {}

impl Similarity<Aac> for Blosum45 {
    fn read_score(&self, code_1: Aac, code_2: Aac) -> SimilarityType {
        aminoacid_data::read_blosum45(code_1, code_2)
    }
}

pub struct Blosum62 {}

impl Similarity<Aac> for Blosum62 {
    fn read_score(&self, code_1: Aac, code_2: Aac) -> SimilarityType {
        aminoacid_data::read_blosum62(code_1, code_2)
    }
}

pub struct Pam160 {}

impl Similarity<Aac> for Pam160 {
    fn read_score(&self, code_1: Aac, code_2: Aac) -> SimilarityType {
        aminoacid_data::read_pam160(code_1, code_2)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const ALL_AAC: [Aac; 20] = [
        Aac::A,
        Aac::C,
        Aac::D,
        Aac::E,
        Aac::F,
        Aac::G,
        Aac::H,
        Aac::I,
        Aac::K,
        Aac::L,
        Aac::M,
        Aac::N,
        Aac::P,
        Aac::Q,
        Aac::R,
        Aac::S,
        Aac::T,
        Aac::V,
        Aac::W,
        Aac::Y,
    ];

    #[test]
    fn check_some_blosum45() {
        let blosum = similarity_builder(AaScoringKind::Blosum45);

        let score_cases = [
            // Diagonal extremes and mid
            (5, Aac::A, Aac::A),
            (8, Aac::Y, Aac::Y),
            (6, Aac::N, Aac::N),
            // Check symmetria
            (-2, Aac::P, Aac::I),
            (-2, Aac::I, Aac::P),
            (-2, Aac::Q, Aac::G),
            (-2, Aac::G, Aac::Q),
        ];
        for (expected, code_1, code_2) in score_cases {
            assert_eq!(expected, blosum.read_score(code_1, code_2))
        }
    }

    #[test]
    fn blosum45_is_complete() {
        let blosum = similarity_builder(AaScoringKind::Blosum45);
        for code_1 in ALL_AAC {
            for code_2 in ALL_AAC {
                blosum.read_score(code_1, code_2);
            }
        }
    }

    #[test]
    fn check_some_blosum62() {
        let blosum = similarity_builder(AaScoringKind::Blosum62);
        let score_cases = [
            // Diagonal extremes and mid
            (4, Aac::A, Aac::A),
            (7, Aac::Y, Aac::Y),
            (6, Aac::N, Aac::N),
            // Check symmetria
            (-3, Aac::P, Aac::I),
            (-3, Aac::I, Aac::P),
            (-2, Aac::Q, Aac::G),
            (-2, Aac::G, Aac::Q),
        ];

        for (expected, code_1, code_2) in score_cases {
            assert_eq!(expected, blosum.read_score(code_1, code_2))
        }
    }

    #[test]
    fn blosum62_is_complete() {
        let blosum = similarity_builder(AaScoringKind::Blosum62);
        for code_1 in ALL_AAC {
            for code_2 in ALL_AAC {
                blosum.read_score(code_1, code_2);
            }
        }
    }

    #[test]
    fn check_some_pam160() {
        let pam = similarity_builder(AaScoringKind::Pam160);
        let score_cases = [
            // Diagonal extremes and mid
            (2, Aac::A, Aac::A),
            (8, Aac::Y, Aac::Y),
            (3, Aac::N, Aac::N),
            // Check symmetria
            (-2, Aac::P, Aac::I),
            (-2, Aac::I, Aac::P),
            (-2, Aac::Q, Aac::G),
            (-2, Aac::G, Aac::Q),
        ];

        for (expected, code_1, code_2) in score_cases {
            assert_eq!(expected, pam.read_score(code_1, code_2))
        }
    }

    #[test]
    fn read_pam160_is_complete() {
        let pam = similarity_builder(AaScoringKind::Pam160);
        for code_1 in ALL_AAC {
            for code_2 in ALL_AAC {
                pam.read_score(code_1, code_2);
            }
        }
    }
}
