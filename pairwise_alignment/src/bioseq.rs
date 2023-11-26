//! Module representing biological sequences

/// IUPAC Amino acid codes
#[derive(PartialEq, Debug, Clone, Copy)]
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
    pub fn from_char(char_code: char) -> Result<Self, SeqError> {
        if !char_code.is_ascii() {
            return Err(SeqError::new(ErrorKind::NonAscii));
        }
        let char_code = char_code.to_ascii_uppercase();

        Self::char_mapping(char_code)
    }

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

    /// Cantor pairing function
    /// The function is bijective and strictly monotonic
    pub fn duple_pairing(code_1: &Self, code_2: &Self) -> u16 {
        let k1 = *code_1 as u16;
        let k2 = *code_2 as u16;
        let k12 = k1 + k2;
        (k12 * (k12 + 1)).div_euclid(2) + k2
    }
}

pub struct Protein {
    sequence: Vec<Aac>,
}

impl Protein {
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
}
