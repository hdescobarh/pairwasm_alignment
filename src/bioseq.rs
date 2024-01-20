//! Data structures representing biological sequences and their building blocks.

use crate::utils::AlignmentUnit;
use std::fmt::Debug;

/// IUPAC Amino acid codes. Represents the basic 20 amino acids.
#[derive(Clone, Copy, PartialEq, Eq)]
#[cfg_attr(test, derive(Debug, PartialOrd, Ord))]
#[repr(u8)]
pub enum Aac {
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    P,
    Q,
    R,
    S,
    T,
    V,
    W,
    Y,
}

impl Aac {
    /// Creates a Aac (amino acid code) from a single character IUPAC code.
    /// The function is case-insensitive. Returns SeqError if the character is not a valid code.
    ///
    /// # Arguments
    /// + `char_code`: - A char representing a valid IUPAC code
    ///
    /// # Examples
    ///
    /// ```
    /// use pairwasm_alignment::bioseq::*;
    /// use std::mem;
    ///
    /// let lys = Aac::from_char('l').unwrap();
    /// assert_eq!(mem::discriminant(&Aac::L), mem::discriminant(&lys));
    /// assert_eq!(
    ///     mem::discriminant(&lys),
    ///     mem::discriminant(&Aac::from_char('L').unwrap())
    /// );
    /// assert!(Aac::from_char('Z').is_err())
    /// ```
    pub fn from_char(char_code: char) -> Result<Self, SeqError> {
        if !char_code.is_ascii() {
            return Err(SeqError::new(ErrorKind::NonAscii));
        }
        let char_code = char_code.to_ascii_uppercase();

        Self::char_mapping(char_code)
    }

    // Contains the map between the valid char values amino acid code and their enum representation
    fn char_mapping(char_code: char) -> Result<Self, SeqError> {
        match char_code {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'D' => Ok(Self::D),
            'E' => Ok(Self::E),
            'F' => Ok(Self::F),
            'G' => Ok(Self::G),
            'H' => Ok(Self::H),
            'I' => Ok(Self::I),
            'K' => Ok(Self::K),
            'L' => Ok(Self::L),
            'M' => Ok(Self::M),
            'N' => Ok(Self::N),
            'P' => Ok(Self::P),
            'Q' => Ok(Self::Q),
            'R' => Ok(Self::R),
            'S' => Ok(Self::S),
            'T' => Ok(Self::T),
            'V' => Ok(Self::V),
            'W' => Ok(Self::W),
            'Y' => Ok(Self::Y),
            _ => Err(SeqError::new(ErrorKind::InvalidCode)),
        }
    }
}

impl AlignmentUnit for Aac {}

/// Trait that allows to biological sequences to expose their content.
pub trait HasSequence<T>
where
    T: Copy + AlignmentUnit,
{
    /// returns the protein or nucleic acid sequence.
    fn seq(&self) -> &Vec<T>;
}

/// Representation of a protein.
pub struct Protein {
    /// Encodes the protein primary structure.
    sequence: Vec<Aac>,
}

impl Protein {
    /// Creates a Protein from a string. The function is case-insensitive.
    /// Returns SeqError if the string contains non-valid IUPAC codes.
    ///
    /// # Arguments
    ///
    /// * `string` - a text containing valid IUPAC amino acid code points. Only accepts ASCII characters.
    ///
    /// # Examples
    ///
    /// ```
    /// use pairwasm_alignment::bioseq::*;
    /// use std::mem;
    /// let protein = Protein::new("pVaGH").unwrap();
    /// let expected_sequence: Vec<Aac> = [Aac::P, Aac::V, Aac::A, Aac::G, Aac::H].to_vec();
    /// for (i, aminoacid) in protein.seq().iter().enumerate() {
    ///     assert_eq!(
    ///         mem::discriminant(aminoacid),
    ///         mem::discriminant(&expected_sequence[i])
    ///     );
    /// }
    /// assert!(Protein::new("pBaGH").is_err())
    /// ```
    pub fn new(string: &str) -> Result<Self, SeqError> {
        if string.is_empty() {
            return Err(SeqError::new(ErrorKind::EmptyString));
        }
        let mut sequence: Vec<Aac> = Vec::new();
        for c in string.chars() {
            sequence.push(Aac::from_char(c)?)
        }
        Ok(Self { sequence })
    }
}

impl HasSequence<Aac> for Protein {
    fn seq(&self) -> &Vec<Aac> {
        &self.sequence
    }
}

#[non_exhaustive]
#[derive(Debug, PartialEq)]
/// A list specifying general error categories of SeqError.
pub enum ErrorKind {
    EmptyString,
    InvalidCode,
    NonAscii,
}

#[derive(Debug)]
/// Error type for operations related to biological sequences manipulation.
pub struct SeqError {
    kind: ErrorKind,
    message: String,
}

impl SeqError {
    fn new(kind: ErrorKind) -> Self {
        let message: String = match kind {
            ErrorKind::InvalidCode => {
                "The string contains a non valid IUPAC code.".to_string()
            }
            ErrorKind::EmptyString => {
                "The string must containt at least one IUPAC code".to_string()
            }
            ErrorKind::NonAscii => {
                "All the IUPAC codes must be ASCII characters.".to_string()
            }
        };

        Self { kind, message }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn creates_protein_from_string() {
        assert_eq!(
            Vec::from([
                Aac::M,
                Aac::N,
                Aac::G,
                Aac::T,
                Aac::E,
                Aac::G,
                Aac::P,
                Aac::N,
                Aac::F,
                Aac::Y,
                Aac::V,
                Aac::P,
            ]),
            Protein::new("MnGTEgPNFyVp").unwrap().sequence
        );
    }

    #[test]
    fn empty_string_to_protein() {
        assert!(Protein::new("").is_err_and(|e| e.kind == ErrorKind::EmptyString))
    }

    #[test]
    fn bad_string_to_protein() {
        // Non-ASCII
        assert!(Protein::new("VTVQＨKKLRT").is_err_and(|e| e.kind == ErrorKind::NonAscii));
        // contains non IUPAC code characters
        assert!(
            Protein::new("VTVQBKKLRT").is_err_and(|e| e.kind == ErrorKind::InvalidCode)
        )
    }

    #[test]
    fn read_sequence_from_external() {
        let protein: Protein = Protein::new("pVaGH").unwrap();
        assert_eq!(
            [Aac::P, Aac::V, Aac::A, Aac::G, Aac::H].to_vec(),
            *protein.seq()
        )
    }
}
