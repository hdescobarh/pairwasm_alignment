//! common data structures and functions used for multiple align algorithms

use crate::{
    bioseq::HasSequence, matrix::Matrix, scoring_schema::ScoringSchema,
    utils::AlignmentUnit,
};
use std::mem::replace;

/// Represent values for backtracking
#[derive(Clone, Copy)]
#[cfg_attr(test, derive(PartialEq, Debug))]
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
    pub fn make_backtrack(top: f32, diagonal: f32, left: f32) -> (BackTrack, f32) {
        let max_score = [top, diagonal, left].into_iter().reduce(f32::max).unwrap();
        let mut which_back: u8 = 0b000;
        for (value, indicator) in [(top, 0b001), (diagonal, 0b010), (left, 0b100)] {
            if value == max_score {
                which_back |= indicator
            }
        }

        // [Empty, Top, Diag., Left, Diag-top, Diag-left, Top-left, Any]
        // [0b000, 0b001, 0b010, 0b100, 0b011, 0b110, 0b101, 0b111]
        let backtrack = match which_back {
            b'\x01' => BackTrack::T(max_score),
            b'\x02' => BackTrack::D(max_score),
            b'\x04' => BackTrack::L(max_score),
            b'\x03' => BackTrack::DT(max_score),
            b'\x06' => BackTrack::DL(max_score),
            b'\x05' => BackTrack::TL(max_score),
            b'\x07' => BackTrack::All(max_score),
            _ => panic!(
                "This must be unreachable. Check match variants.\
             Max score: {}. Back Value {}",
                max_score, which_back
            ),
        };
        (backtrack, max_score)
    }

    /// from an entry matrix cell tracks all the paths
    ///
    /// * `cutoff_score`: a lower bound for the score of a single matrix element.
    /// If the cell contains a score equal or lower than cutoff, then the backtrack
    /// in that branch stops. You can use f32::NEG_INFINITY if do not want to set any cutoff
    pub fn backtracking(
        matrix: &Matrix<BackTrack>,
        init_row: usize,
        init_col: usize,
        cutoff_score: f32,
    ) -> Vec<Vec<[usize; 2]>> {
        let mut paths: Vec<Vec<[usize; 2]>> = Vec::new();
        let mut pending_stack: Vec<Vec<[usize; 2]>> = Vec::new();
        let current_path: Vec<[usize; 2]> = vec![[init_row, init_col]];
        Self::find_paths(
            matrix,
            current_path,
            &mut pending_stack,
            &mut paths,
            cutoff_score,
        );
        paths
    }

    fn find_paths(
        matrix: &Matrix<BackTrack>,
        mut current_path: Vec<[usize; 2]>,
        pending_stack: &mut Vec<Vec<[usize; 2]>>,
        paths: &mut Vec<Vec<[usize; 2]>>,
        // The minimum allowed score of a single node. The path ends prematurely if
        // a score lower or equal is found.
        cutoff_score: f32,
    ) {
        let [row, col] = *current_path.last().unwrap();
        let (indicator, score) = Self::decompose(matrix[[row, col]]);
        if ((row == 0) && (col == 0)) || (score <= cutoff_score) {
            match pending_stack.pop() {
                Some(next_path) => {
                    let old_path = replace(&mut current_path, next_path);
                    paths.push(old_path);
                }
                None => {
                    paths.push(current_path);
                    return;
                }
            }
        } else {
            match indicator {
                // T
                b'\x01' => current_path.push([row - 1, col]),
                //D
                b'\x02' => current_path.push([row - 1, col - 1]),
                //L
                b'\x04' => current_path.push([row, col - 1]),
                //DT
                b'\x03' => {
                    let mut branch = current_path.clone();
                    branch.push([row - 1, col]);
                    pending_stack.push(branch);
                    current_path.push([row - 1, col - 1]);
                    },
                //DL
                b'\x06' => {
                    let mut branch = current_path.clone();
                    branch.push([row, col - 1]);
                    pending_stack.push(branch);
                    current_path.push([row - 1, col - 1]);
                }
                //TL
                b'\x05' => {
                    let mut branch = current_path.clone();
                    branch.push([row, col - 1]);
                    pending_stack.push(branch);
                    current_path.push([row - 1, col]);
                }
                //All
                b'\x07' => {
                    let mut branch = current_path.clone();
                    branch.push([row, col - 1]);
                    pending_stack.push(branch);

                    let mut branch = current_path.clone();
                    branch.push([row - 1, col]);
                    pending_stack.push(branch);

                    current_path.push([row - 1, col - 1])}

                _ => panic!(
                    "Empty at [{row}, {col}]. Any implementation must remove all Empty from the matrix."
                ),
            };
        };
        Self::find_paths(matrix, current_path, pending_stack, paths, cutoff_score);
    }

    /// Separates the BackTrack from its associated value. If BackTrack::Empty, returns NAN.
    fn decompose(backtrack: BackTrack) -> (u8, f32) {
        match backtrack {
            BackTrack::Empty => (0b000, f32::NAN),
            BackTrack::T(v) => (0b001, v),
            BackTrack::D(v) => (0b010, v),
            BackTrack::L(v) => (0b100, v),
            BackTrack::DT(v) => (0b011, v),
            BackTrack::DL(v) => (0b110, v),
            BackTrack::TL(v) => (0b101, v),
            BackTrack::All(v) => (0b111, v),
        }
    }
}

/// Represents a single alignment
#[cfg_attr(test, derive(Debug))]
pub struct AlignmentSequence<A>
where
    A: AlignmentUnit,
{
    pairs: Vec<[Option<A>; 2]>,
}

impl<A> AlignmentSequence<A>
where
    A: AlignmentUnit,
{
    pub fn new(
        // remember this is shifted: [i, j] means left.seq[i-1] and top.seq[j-1]
        backtrack_path: Vec<[usize; 2]>,
        sequence_left: &(impl HasSequence<A> + ?Sized),
        sequence_top: &(impl HasSequence<A> + ?Sized),
    ) -> Self
    where
        A: AlignmentUnit,
    {
        let mut pairs: Vec<[Option<A>; 2]> = Vec::with_capacity(backtrack_path.len());

        for index in (0..backtrack_path.len() - 1).rev() {
            let [row, col] = backtrack_path[index];
            let [last_row, last_col] = backtrack_path[index + 1];

            let next = {
                if row != last_row && col != last_col {
                    [
                        Some(sequence_left.seq()[row - 1]),
                        Some(sequence_top.seq()[col - 1]),
                    ]
                } else if row != last_row && col == last_col {
                    [Some(sequence_left.seq()[row - 1]), None]
                } else if row == last_row && col != last_col {
                    [None, Some(sequence_top.seq()[col - 1])]
                } else {
                    panic!("This must be unreachable. Does not exist a path such that repeats indices.")
                }
            };
            pairs.push(next)
        }

        Self { pairs }
    }

    pub fn read(&self) -> &Vec<[Option<A>; 2]> {
        &self.pairs
    }
}

// Be aware this implimentation is intended to be used with Affine gap models and
// their subsets; i.e., Linear and constant models.
pub trait AffineTransversalOrder<A>
where
    A: AlignmentUnit,
{
    fn diagonal_score(
        sequence_left: &(impl HasSequence<A> + ?Sized),
        sequence_top: &(impl HasSequence<A> + ?Sized),
        scoring_schema: &Box<dyn ScoringSchema<A>>,
        matrix: &Matrix<BackTrack>,
        i: usize,
        j: usize,
    ) -> f32 {
        // Read the sequences i,j element. Remember the Matrix has (n+1)(m+1) elements, with the
        // extra row and colum at the start.
        let left_alignable: A = sequence_left.seq()[i - 1];
        let top_alignable: A = sequence_top.seq()[j - 1];
        let score_ij = scoring_schema.get_score(left_alignable, top_alignable);
        let value = match matrix[[i - 1, j - 1]] {
            BackTrack::T(v) => v,
            BackTrack::D(v) => v,
            BackTrack::L(v) => v,
            BackTrack::DT(v) => v,
            BackTrack::DL(v) => v,
            BackTrack::TL(v) => v,
            BackTrack::All(v) => v,
            BackTrack::Empty => {
                panic!("This must be unreachable.Check the transversal order.")
            }
        };
        value + score_ij as f32
    }

    fn top_score(
        scoring_schema: &Box<dyn ScoringSchema<A>>,
        matrix: &Matrix<BackTrack>,
        i: usize,
        j: usize,
    ) -> f32 {
        // i-1, j
        match matrix[[i - 1, j]] {
            // top_gap + top_gap is an extension
            BackTrack::T(v) => v - scoring_schema.get_extend(),
            // not(top_gap) + top_gap is and opening
            BackTrack::D(v) => v - scoring_schema.get_function(1),
            BackTrack::L(v) => v - scoring_schema.get_function(1),
            BackTrack::DL(v) => v - scoring_schema.get_function(1),
            // Max(v - extend_existing_gap, v - add_new_gap) = v - extend_gap
            // because extend_existing_gap <= add_new_gap
            BackTrack::DT(v) => v - scoring_schema.get_extend(),
            BackTrack::TL(v) => v - scoring_schema.get_extend(),
            BackTrack::All(v) => v - scoring_schema.get_extend(),
            BackTrack::Empty => {
                panic!("This must be unreachable.Check the transversal order.")
            }
        }
    }

    fn left_score(
        scoring_schema: &Box<dyn ScoringSchema<A>>,
        matrix: &Matrix<BackTrack>,
        i: usize,
        j: usize,
    ) -> f32 {
        // i, j-1
        match matrix[[i, j - 1]] {
            // left_gap + left_gap is and extension
            BackTrack::L(v) => v - scoring_schema.get_extend(),
            // not(left_gap) + left_gap is a new gap
            BackTrack::T(v) => v - scoring_schema.get_function(1),
            BackTrack::D(v) => v - scoring_schema.get_function(1),
            BackTrack::DT(v) => v - scoring_schema.get_function(1),
            // Max(v - extend_existing_gap, v - add_new_gap) = v - extend_gap
            // because extend_existing_gap <= add_new_gap
            BackTrack::DL(v) => v - scoring_schema.get_extend(),
            BackTrack::TL(v) => v - scoring_schema.get_extend(),
            BackTrack::All(v) => v - scoring_schema.get_extend(),
            BackTrack::Empty => {
                panic!("This must be unreachable.Check the transversal order.")
            }
        }
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use super::{AlignmentSequence, BackTrack};
    use crate::{
        bioseq::{Aac, Protein},
        matrix::Matrix,
    };

    #[test]
    fn backtrack_from_scores() {
        // top i
        let test_cases = [
            (BackTrack::T(1.0), [1.0, 0.0, 0.0]),
            (BackTrack::T(1.0), [1.0, -1.0, 0.0]),
            (BackTrack::T(1.0), [1.0, 0.0, -1.0]),
            (BackTrack::D(1.0), [0.0, 1.0, 0.0]),
            (BackTrack::D(1.0), [-1.0, 1.0, 0.0]),
            (BackTrack::D(1.0), [0.0, 1.0, -1.0]),
            (BackTrack::L(1.0), [0.0, 0.0, 1.0]),
            (BackTrack::L(1.0), [-1.0, 0.0, 1.0]),
            (BackTrack::L(1.0), [0.0, -1.0, 1.0]),
            (BackTrack::DT(1.0), [1.0, 1.0, 0.0]),
            (BackTrack::DT(1.0), [1.0, 1.0, -1.0]),
            (BackTrack::DL(1.0), [0.0, 1.0, 1.0]),
            (BackTrack::DL(1.0), [-1.0, 1.0, 1.0]),
            (BackTrack::TL(1.0), [1.0, 0.0, 1.0]),
            (BackTrack::TL(1.0), [1.0, -1.0, 1.0]),
            (BackTrack::All(1.0), [1.0, 1.0, 1.0]),
            (BackTrack::All(0.0), [0.0, 0.0, 0.0]),
            (BackTrack::All(-1.0), [-1.0, -1.0, -1.0]),
        ];

        for (expected, [top, diagonal, left]) in test_cases {
            assert_eq!(
                expected,
                BackTrack::make_backtrack(top, diagonal, left).0,
                "Failed at case [{},{},{}]",
                top,
                diagonal,
                left
            )
        }
    }

    #[test]
    fn matrix_backtrack_one_path_any_start() {
        let container = vec![
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::L(0.0),
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
        ];
        let matrix = Matrix::from_vec(container, 4, 4);

        let test_cases = [
            ([3, 3], vec![vec![[3, 3], [2, 2], [1, 1], [0, 0]]]), // start at last cell
            // start any other cell
            ([2, 3], vec![vec![[2, 3], [2, 2], [1, 1], [0, 0]]]),
            ([3, 2], vec![vec![[3, 2], [2, 2], [1, 1], [0, 0]]]),
        ];

        for (start, expected_path) in test_cases {
            let [init_row, init_col] = start;
            let actual_path =
                BackTrack::backtracking(&matrix, init_row, init_col, f32::NEG_INFINITY);
            assert_eq!(
                expected_path, actual_path,
                "Failed with starting cell [{init_row}, {init_col}]."
            )
        }
    }

    #[test]
    #[should_panic(
        expected = "Empty at [2, 2]. Any implementation must remove all Empty from the matrix."
    )]
    fn matrix_backtrack_fails_bad_cell_in_path() {
        let container = vec![
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::L(0.0),
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::T(0.0),
            BackTrack::Empty,
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
        ];
        let matrix = Matrix::from_vec(container, 4, 4);
        BackTrack::backtracking(&matrix, 3, 3, f32::NEG_INFINITY);
    }

    #[test]
    fn matrix_backtrack_multiple_paths() {
        let test_cases = [
            (
                vec![
                    BackTrack::D(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::L(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::DT(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                ],
                HashSet::from([
                    vec![
                        [6, 3],
                        [5, 2],
                        [4, 1],
                        [4, 0],
                        [3, 0],
                        [2, 0],
                        [1, 0],
                        [0, 0],
                    ],
                    vec![[6, 3], [5, 2], [4, 2], [3, 2], [2, 2], [1, 1], [0, 0]],
                ]),
            ),
            (
                vec![
                    BackTrack::D(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::DL(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::DL(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                ],
                HashSet::from([
                    vec![[6, 3], [5, 2], [4, 1], [3, 1], [2, 1], [1, 0], [0, 0]],
                    vec![
                        [6, 3],
                        [5, 2],
                        [5, 1],
                        [4, 0],
                        [3, 0],
                        [2, 0],
                        [1, 0],
                        [0, 0],
                    ],
                ]),
            ),
            (
                vec![
                    BackTrack::D(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::L(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::L(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::TL(0.0),
                ],
                HashSet::from([
                    vec![
                        [6, 3],
                        [5, 3],
                        [4, 3],
                        [3, 2],
                        [2, 1],
                        [2, 0],
                        [1, 0],
                        [0, 0],
                    ],
                    vec![
                        [6, 3],
                        [6, 2],
                        [5, 1],
                        [4, 1],
                        [3, 0],
                        [2, 0],
                        [1, 0],
                        [0, 0],
                    ],
                ]),
            ),
            (
                vec![
                    BackTrack::D(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::L(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::L(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::All(0.0),
                    BackTrack::T(0.0),
                    BackTrack::T(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                    BackTrack::D(0.0),
                ],
                HashSet::from([
                    vec![[6, 3], [5, 2], [4, 1], [3, 1], [2, 0], [1, 0], [0, 0]],
                    vec![[6, 3], [5, 2], [4, 2], [3, 2], [2, 2], [1, 1], [0, 0]],
                    vec![
                        [6, 3],
                        [5, 2],
                        [5, 1],
                        [4, 0],
                        [3, 0],
                        [2, 0],
                        [1, 0],
                        [0, 0],
                    ],
                ]),
            ),
        ];

        let mut counter = 0;
        for (container, expected_paths) in test_cases {
            counter += 1;
            let rows = 7;
            let cols = 4;
            let matrix = Matrix::from_vec(container, rows, cols);
            let [init_row, init_col] = [rows - 1, cols - 1];
            let actual_paths: HashSet<Vec<[usize; 2]>> =
                BackTrack::backtracking(&matrix, init_row, init_col, f32::NEG_INFINITY)
                    .into_iter()
                    .collect();

            // Didn't missed any path
            let diff_missing: HashSet<_> =
                expected_paths.difference(&actual_paths).collect();
            // Didn't get an invalid path
            let diff_extra: HashSet<_> =
                actual_paths.difference(&expected_paths).collect();

            assert!(
                diff_missing.is_empty() && diff_extra.is_empty(),
                "Error at test case {}.\n- Extra paths: {:?}.\n- Missing paths: {:?}.",
                counter,
                diff_extra,
                diff_missing
            );
        }
    }
    #[test]
    fn matrix_backtrack_nested_multiple_paths() {
        let container = vec![
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::L(0.0),
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::D(0.0),
            BackTrack::D(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::T(0.0),
            BackTrack::L(0.0),
            BackTrack::T(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::T(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::DL(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::T(0.0),
            BackTrack::T(0.0),
            BackTrack::T(0.0),
            BackTrack::D(0.0),
            BackTrack::L(0.0),
            BackTrack::TL(0.0),
        ];

        let expected_paths = HashSet::from([
            vec![[6, 3], [5, 3], [4, 3], [3, 2], [2, 1], [1, 0], [0, 0]],
            vec![
                [6, 3],
                [5, 3],
                [4, 3],
                [4, 2],
                [4, 1],
                [3, 0],
                [2, 0],
                [1, 0],
                [0, 0],
            ],
            vec![
                [6, 3],
                [6, 2],
                [6, 1],
                [5, 0],
                [4, 0],
                [3, 0],
                [2, 0],
                [1, 0],
                [0, 0],
            ],
        ]);

        let matrix = Matrix::from_vec(container, 7, 4);
        let [init_row, init_col] = [6, 3];
        let actual_path: HashSet<Vec<[usize; 2]>> =
            BackTrack::backtracking(&matrix, init_row, init_col, f32::NEG_INFINITY)
                .into_iter()
                .collect();

        // Didn't missed any path
        let diff_missing: HashSet<_> = expected_paths.difference(&actual_path).collect();
        // Didn't get an invalid path
        let diff_extra: HashSet<_> = actual_path.difference(&expected_paths).collect();

        assert!(
            diff_missing.is_empty() && diff_extra.is_empty(),
            "- Extra paths: {:?}.\n- Missing paths: {:?}.",
            diff_extra,
            diff_missing
        );
    }

    #[test]
    fn matrix_backtrack_one_path_cutoff() {
        let container = vec![
            BackTrack::D(1.0),
            BackTrack::L(1.0),
            BackTrack::L(1.0),
            BackTrack::L(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::D(1.0),
            BackTrack::L(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::D(1.0),
            BackTrack::L(1.0),
            BackTrack::T(0.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::L(1.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::D(10.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::L(1.0),
            BackTrack::T(1.0),
            BackTrack::T(1.0),
        ];
        let matrix = Matrix::from_vec(container, 7, 4);

        let [init_row, init_col] = [5, 2];
        let actual_path = BackTrack::backtracking(&matrix, init_row, init_col, 0.0);
        let expected_path = vec![vec![[5, 2], [4, 1], [3, 0]]];

        assert_eq!(
            expected_path, actual_path,
            "Failed with starting cell [{init_row}, {init_col}]."
        )
    }

    #[test]
    fn matrix_backtrack_multiple_paths_cutoff() {
        let container = vec![
            BackTrack::D(1.0),
            BackTrack::L(0.0),
            BackTrack::L(1.0),
            BackTrack::L(1.0),
            BackTrack::T(0.0),
            BackTrack::D(1.0),
            BackTrack::D(1.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::D(10.0),
            BackTrack::T(0.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::L(1.0),
            BackTrack::DL(1.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::T(1.0),
            BackTrack::T(10.0),
            BackTrack::T(1.0),
            BackTrack::D(1.0),
            BackTrack::L(1.0),
            BackTrack::TL(1.0),
        ];
        let matrix = Matrix::from_vec(container, 7, 4);

        let expected_paths = HashSet::from([
            vec![[5, 3], [4, 3], [3, 2], [2, 1], [1, 0]],
            vec![[5, 3], [4, 3], [4, 2], [4, 1], [3, 0]],
            vec![[2, 3], [1, 2], [0, 1]],
        ]);

        let mut actual_paths: HashSet<Vec<[usize; 2]>> =
            BackTrack::backtracking(&matrix, 5, 3, 0.0)
                .into_iter()
                .collect();

        actual_paths.extend(BackTrack::backtracking(&matrix, 2, 3, 0.0));

        // Didn't missed any path
        let diff_missing: HashSet<_> = expected_paths.difference(&actual_paths).collect();
        // Didn't get an invalid path
        let diff_extra: HashSet<_> = actual_paths.difference(&expected_paths).collect();

        assert!(
            diff_missing.is_empty() && diff_extra.is_empty(),
            "- Extra paths: {:?}.\n- Missing paths: {:?}.",
            diff_extra,
            diff_missing
        );
    }

    #[test]
    fn alignment_sequence_no_gap() {
        let sequence_left = Protein::new("MVLSPADKT").unwrap();
        let sequence_top = Protein::new("MVLSGEDKS").unwrap();

        let backtrack_path = vec![
            [9, 9],
            [8, 8],
            [7, 7],
            [6, 6],
            [5, 5],
            [4, 4],
            [3, 3],
            [2, 2],
            [1, 1],
            [0, 0],
        ];
        let expected_alignment = AlignmentSequence {
            pairs: vec![
                [Some(Aac::M), Some(Aac::M)],
                [Some(Aac::V), Some(Aac::V)],
                [Some(Aac::L), Some(Aac::L)],
                [Some(Aac::S), Some(Aac::S)],
                [Some(Aac::P), Some(Aac::G)],
                [Some(Aac::A), Some(Aac::E)],
                [Some(Aac::D), Some(Aac::D)],
                [Some(Aac::K), Some(Aac::K)],
                [Some(Aac::T), Some(Aac::S)],
            ],
        };
        let actual_alignment =
            AlignmentSequence::new(backtrack_path, &sequence_left, &sequence_top);

        assert_eq!(expected_alignment.pairs.len(), actual_alignment.pairs.len());
        for p in 0..expected_alignment.pairs.len() {
            assert!(
                (expected_alignment.pairs[p][0] == actual_alignment.pairs[p][0])
                    && (expected_alignment.pairs[p][1] == actual_alignment.pairs[p][1]),
                "Error at index {}. Expected: {:?}. Got: {:?}.",
                p,
                expected_alignment.pairs,
                actual_alignment.pairs
            )
        }
    }

    #[test]
    fn alignment_sequence_left_gap() {
        let sequence_left = Protein::new("KVGAHAGEYA").unwrap();
        let sequence_top = Protein::new("KIGGHGAEYGA").unwrap();
        let backtrack_path = vec![
            [10, 11],
            [9, 10],
            [9, 9],
            [8, 8],
            [7, 7],
            [6, 6],
            [5, 5],
            [4, 4],
            [3, 3],
            [2, 2],
            [1, 1],
            [0, 0],
        ];

        let expected_alignment = AlignmentSequence {
            pairs: vec![
                [Some(Aac::K), Some(Aac::K)],
                [Some(Aac::V), Some(Aac::I)],
                [Some(Aac::G), Some(Aac::G)],
                [Some(Aac::A), Some(Aac::G)],
                [Some(Aac::H), Some(Aac::H)],
                [Some(Aac::A), Some(Aac::G)],
                [Some(Aac::G), Some(Aac::A)],
                [Some(Aac::E), Some(Aac::E)],
                [Some(Aac::Y), Some(Aac::Y)],
                [None, Some(Aac::G)],
                [Some(Aac::A), Some(Aac::A)],
            ],
        };
        let actual_alignment =
            AlignmentSequence::new(backtrack_path, &sequence_left, &sequence_top);

        assert_eq!(expected_alignment.pairs.len(), actual_alignment.pairs.len());
        for p in 0..expected_alignment.pairs.len() {
            assert!(
                (expected_alignment.pairs[p][0] == actual_alignment.pairs[p][0])
                    && (expected_alignment.pairs[p][1] == actual_alignment.pairs[p][1]),
                "Error at index {}. Expected: {:?}. Got: {:?}.",
                p,
                expected_alignment.pairs,
                actual_alignment.pairs
            )
        }
    }

    #[test]
    fn alignment_sequence_top_gap() {
        let sequence_left = Protein::new("KIGGHGAEYGA").unwrap();
        let sequence_top = Protein::new("KVGAHAGEYA").unwrap();
        let backtrack_path = vec![
            [11, 10],
            [10, 9],
            [9, 9],
            [8, 8],
            [7, 7],
            [6, 6],
            [5, 5],
            [4, 4],
            [3, 3],
            [2, 2],
            [1, 1],
            [0, 0],
        ];

        let expected_alignment = AlignmentSequence {
            pairs: vec![
                [Some(Aac::K), Some(Aac::K)],
                [Some(Aac::I), Some(Aac::V)],
                [Some(Aac::G), Some(Aac::G)],
                [Some(Aac::G), Some(Aac::A)],
                [Some(Aac::H), Some(Aac::H)],
                [Some(Aac::G), Some(Aac::A)],
                [Some(Aac::A), Some(Aac::G)],
                [Some(Aac::E), Some(Aac::E)],
                [Some(Aac::Y), Some(Aac::Y)],
                [Some(Aac::G), None],
                [Some(Aac::A), Some(Aac::A)],
            ],
        };
        let actual_alignment =
            AlignmentSequence::new(backtrack_path, &sequence_left, &sequence_top);

        assert_eq!(expected_alignment.pairs.len(), actual_alignment.pairs.len());
        for p in 0..expected_alignment.pairs.len() {
            assert!(
                (expected_alignment.pairs[p][0] == actual_alignment.pairs[p][0])
                    && (expected_alignment.pairs[p][1] == actual_alignment.pairs[p][1]),
                "Error at index {}. Expected: {:?}. Got: {:?}.",
                p,
                expected_alignment.pairs,
                actual_alignment.pairs
            )
        }
    }
}
