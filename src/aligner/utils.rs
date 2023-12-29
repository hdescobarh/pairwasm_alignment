//! common data structures and functions used for multiple align algorithms

/// Represent values for backtracking
#[derive(Clone)]
pub enum BackTrack {
    /// Used for initialize collections
    Empty,
    /// Top: Gap at top sequence
    T(f32),
    /// Top-left: Match/Mismatch
    D(f32),
    /// Left: Gap at left sequence
    L(f32),
    /// Top-left and top
    DT(f32),
    /// Top-left and left
    DL(f32),
    /// Top and left
    TL(f32),
    /// All: top, left and top-left
    All(f32),
}

impl BackTrack {
    /// generate the back track direction from the scores
    pub fn make_backtrack(diagonal: f32, top: f32, left: f32) -> BackTrack {
        let max_score = [top, diagonal, left].into_iter().reduce(f32::min).unwrap();
        let mut which_back: u8 = 0b000;
        for (value, indicator) in [(top, 0b001), (diagonal, 0b010), (left, 0b100)] {
            if value == max_score {
                which_back &= indicator
            }
        }

        // [Empty, Top, Diag., Left, Diag-top, Diag-left, Top-left, Any]
        // [0b000, 0b001, 0b010, 0b100, 0b011, 0b110, 0b101, 0b111]
        match which_back {
            b'\x00' => BackTrack::Empty,
            b'\x01' => BackTrack::T(max_score),
            b'\x02' => BackTrack::D(max_score),
            b'\x04' => BackTrack::L(max_score),
            b'\x03' => BackTrack::DT(max_score),
            b'\x06' => BackTrack::DL(max_score),
            b'\x05' => BackTrack::TL(max_score),
            b'\x07' => BackTrack::All(max_score),
            _ => panic!("This must be unreachable. Check match variants."),
        }
    }
}
