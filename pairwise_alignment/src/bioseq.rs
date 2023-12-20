//! Module representing biological sequences

use std::fmt::Debug;

/// IUPAC Amino acid codes
/// Represents the basic 20 amino acids
#[derive(PartialEq, Debug, Clone, Copy, Eq, PartialOrd, Ord)]
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
    /// The function is case-insensitive. Returns SeqError if the character is not a valid code
    ///
    /// # Arguments
    /// + `char_code`: - A char representing a valid IUPAC code
    ///
    /// # Examples
    ///
    /// ```
    /// use pairwise_alignment::bioseq::*;
    ///
    /// let lys = Aac::from_char('l').unwrap();
    /// assert_eq!(Aac::L, lys);
    /// assert_eq!(lys, Aac::from_char('L').unwrap());
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

    /// It is a map (amino acid code, amino acid code) ↦ integer
    /// Represents each Aac duple as a unique integer identifier.
    /// The current implementation is the cantor pairing function,
    /// which is bijective and strictly monotonic.
    pub fn duple_pairing(code_1: &Self, code_2: &Self) -> u16 {
        let k1 = *code_1 as u16;
        let k2 = *code_2 as u16;
        let k12 = k1 + k2;
        (k12 * (k12 + 1)).div_euclid(2) + k2
    }
}

pub trait HasSequence<T> {
    /// returns the protein or nucleic acid sequence
    fn seq(&self) -> &Vec<T>
    where
        T: Eq + Ord + Debug + Copy + Clone;
}

/// Representation of a protein
pub struct Protein {
    /// Encodes the protein primary structure
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
    /// use pairwise_alignment::bioseq::*;
    /// let protein = Protein::new("pVaGH").unwrap();
    /// assert_eq!([Aac::P, Aac::V, Aac::A, Aac::G, Aac::H].to_vec(),*protein.seq());
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
    fn seq(&self) -> &Vec<Aac>
    where
        Aac: Eq + Ord + Debug + Copy + Clone,
    {
        &self.sequence
    }
}

#[non_exhaustive]
#[derive(Debug, PartialEq)]
pub enum ErrorKind {
    EmptyString,
    InvalidCode,
    NonAscii,
}

#[derive(Debug)]
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
    fn aminoacid_code_pairing() {
        // Cantor pairing function
        // (0,0) ↦ 0
        assert_eq!(0, Aac::duple_pairing(&Aac::A, &Aac::A));
        // (4,3) ↦ 31
        assert_eq!(31, Aac::duple_pairing(&Aac::F, &Aac::E));
        // (5,14) ↦ 204
        assert_eq!(204, Aac::duple_pairing(&Aac::G, &Aac::R));
        // (16,15) ↦ 511
        assert_eq!(511, Aac::duple_pairing(&Aac::T, &Aac::S));
        // (19,19) ↦ 760
        assert_eq!(760, Aac::duple_pairing(&Aac::Y, &Aac::Y));
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
