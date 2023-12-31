//! common data structures and functions used for multiple align algorithms

use crate::matrix::Matrix;
use std::mem::replace;

/// Represent values for backtracking
#[derive(Clone)]
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
    pub fn make_backtrack(top: f32, diagonal: f32, left: f32) -> BackTrack {
        let max_score = [top, diagonal, left].into_iter().reduce(f32::max).unwrap();
        let mut which_back: u8 = 0b000;
        for (value, indicator) in [(top, 0b001), (diagonal, 0b010), (left, 0b100)] {
            if value == max_score {
                which_back |= indicator
            }
        }

        // [Empty, Top, Diag., Left, Diag-top, Diag-left, Top-left, Any]
        // [0b000, 0b001, 0b010, 0b100, 0b011, 0b110, 0b101, 0b111]
        match which_back {
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
        }
    }

    /// from an entry matrix cell tracks all the paths
    pub fn backtracking(
        matrix: &Matrix<BackTrack>,
        init_row: usize,
        init_col: usize,
    ) -> Vec<Vec<[usize; 2]>> {
        let mut paths: Vec<Vec<[usize; 2]>> = Vec::new();
        let mut pending_stack: Vec<Vec<[usize; 2]>> = Vec::new();
        let current_path: Vec<[usize; 2]> = vec![[init_row, init_col]];
        Self::find_paths(matrix, current_path, &mut pending_stack, &mut paths);
        paths
    }

    fn find_paths(
        matrix: &Matrix<BackTrack>,
        mut current_path: Vec<[usize; 2]>,
        pending_stack: &mut Vec<Vec<[usize; 2]>>,
        paths: &mut Vec<Vec<[usize; 2]>>,
    ) {
        let [row, col] = *current_path.last().unwrap();
        if (row == 0) && (col == 0) {
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
            match matrix[[row, col]] {
                BackTrack::T(_) => current_path.push([row - 1, col]),
                BackTrack::D(_) => current_path.push([row - 1, col - 1]),
                BackTrack::L(_) => current_path.push([row, col - 1]),
                BackTrack::DT(_) => {
                    let mut branch = current_path.clone();
                    branch.push([row - 1, col]);
                    pending_stack.push(branch);
                    current_path.push([row - 1, col - 1]);
                }
                BackTrack::DL(_) => {
                    let mut branch = current_path.clone();
                    branch.push([row, col - 1]);
                    pending_stack.push(branch);
                    current_path.push([row - 1, col - 1]);
                }
                BackTrack::TL(_) => {
                    let mut branch = current_path.clone();
                    branch.push([row, col - 1]);
                    pending_stack.push(branch);
                    current_path.push([row - 1, col]);
                }
                BackTrack::All(_) => {
                    let mut branch = current_path.clone();
                    branch.push([row, col - 1]);
                    pending_stack.push(branch);

                    let mut branch = current_path.clone();
                    branch.push([row - 1, col]);
                    pending_stack.push(branch);

                    current_path.push([row - 1, col - 1])
                }
                BackTrack::Empty => panic!(
                    "Empty at [{row}, {col}]. Any implementation must remove all Empty from the matrix."
                ),
            };
        };
        Self::find_paths(matrix, current_path, pending_stack, paths);
    }
}

#[cfg(test)]
mod test {
    use std::collections::{HashMap, HashSet};

    use super::BackTrack;
    use crate::matrix::Matrix;

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
                BackTrack::make_backtrack(top, diagonal, left),
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
            let actual_path = BackTrack::backtracking(&matrix, init_row, init_col);
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
        BackTrack::backtracking(&matrix, 3, 3);
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
                BackTrack::backtracking(&matrix, init_row, init_col)
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
            BackTrack::backtracking(&matrix, init_row, init_col)
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
}
