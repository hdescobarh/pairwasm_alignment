//! Module representing biological sequences

const IUPAC_AMINOACIDS: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

pub struct Protein {
    sequence: Vec<char>,
}

impl Protein {
    pub fn new(string: &str) -> Result<Self, SeqError> {
        if string.is_empty() {
            return Err(SeqError::new(ErrorKind::EmptyString));
        }

        let mut sequence: Vec<char> = Vec::new();
        for c in string.chars() {
            if !c.is_ascii() {
                return Err(SeqError::new(ErrorKind::NonAscii));
            }

            let upper_case = c.to_ascii_uppercase();

            if !IUPAC_AMINOACIDS.contains(&upper_case) {
                return Err(SeqError::new(ErrorKind::InvalidCode));
            }
            sequence.push(upper_case)
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
            Vec::from(['M', 'N', 'G', 'T', 'E', 'G', 'P', 'N', 'F', 'Y', 'V', 'P',]),
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
