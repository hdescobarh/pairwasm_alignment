//! Module representing biological sequences

/// IUPAC Amino acid codes
#[derive(PartialEq, Debug)]
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
            ErrorKind::InvalidCode => "The string contains a non valid IUPAC code.".to_string(),
            ErrorKind::EmptyString => {
                "The string must containt at least one IUPAC code".to_string()
            }
            ErrorKind::NonAscii => "All the IUPAC codes must be ASCII characters.".to_string(),
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
        assert!(Protein::new("VTVQBKKLRT").is_err_and(|e| e.kind == ErrorKind::InvalidCode))
    }
}
