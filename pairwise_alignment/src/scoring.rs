//! Alignment scoring schemas

use crate::bioseq::Aac;

pub trait HasName {
    fn name() -> &'static str;
}

/// Amino acid scoring Schema
pub trait AaSchema {
    /// Returns the score value s_{ij}, with i ≡ code_1 and j ≡ code_2.
    fn get_score(id: u8, code_1: &Aac, code_2: &Aac) -> &'static i8;
    fn get_table(id: u8) -> &'static [(u16, i8)];
}

/// Trait for Amino acid symmetric score tables
/// Schemas that are symmetric only needs the upper triangular matrix
/// Schema keys must be ordered or the binary search will be meaningless
///
/// The logic is that if A_{n,m} is a symmetric matrix, then a_{i,j} = a_{j, i}. So,
/// We only need to store the values in the diagonal plus values above or bellow it.
/// The choice of the upper triangle is arbitrary.
trait AaSymTable {
    fn search(schema: &'static [(u16, i8)], code_1: &Aac, code_2: &Aac) -> &'static i8 {
        let mut param = [code_1, code_2];
        param.sort();
        let key: u16 = Aac::duple_pairing(param[0], param[1]);
        match schema.binary_search_by(|(k, _)| k.cmp(&key)) {
            Ok(index) => &schema[index].1,
            Err(_) => panic!("Index not found. The scoring schema may be incomplete."),
        }
    }
}

/// BLOSUM substitution matrices
pub struct Blosum {}

impl HasName for Blosum {
    fn name() -> &'static str {
        "BLOSUM"
    }
}

impl AaSchema for Blosum {
    fn get_score(id: u8, code_1: &Aac, code_2: &Aac) -> &'static i8 {
        Blosum::search(Blosum::get_table(id), code_1, code_2)
    }

    fn get_table(id: u8) -> &'static [(u16, i8)] {
        match id {
            62 => BLOSUM62,
            _ => panic!("BLOSUM{id} is not implemented."),
        }
    }
}

impl AaSymTable for Blosum {}

static BLOSUM62: &[(u16, i8)] = &[
    (0, 4),
    (2, 0),
    (4, 9),
    (5, -2),
    (8, -3),
    (9, -1),
    (12, 6),
    (13, -4),
    (14, -2),
    (18, 2),
    (19, -2),
    (20, 0),
    (24, 5),
    (25, -3),
    (26, -3),
    (27, -2),
    (32, -3),
    (33, -1),
    (34, -3),
    (35, -1),
    (40, 6),
    (41, -2),
    (42, -1),
    (43, -1),
    (44, -1),
    (50, -3),
    (51, 0),
    (52, -3),
    (53, -3),
    (54, -1),
    (60, 6),
    (61, -1),
    (62, -3),
    (63, -1),
    (64, -1),
    (65, -1),
    (72, -2),
    (73, 0),
    (74, 1),
    (75, -4),
    (76, -1),
    (77, -2),
    (84, 8),
    (85, -4),
    (86, -3),
    (87, -3),
    (88, -3),
    (89, -3),
    (90, -1),
    (98, -3),
    (99, -2),
    (100, 0),
    (101, -2),
    (102, 1),
    (103, -3),
    (104, -1),
    (112, 4),
    (113, -1),
    (114, -4),
    (115, 0),
    (116, 0),
    (117, -1),
    (118, -3),
    (119, -1),
    (128, -3),
    (129, -3),
    (130, -3),
    (131, -3),
    (132, -1),
    (133, 0),
    (134, -3),
    (135, 1),
    (144, 5),
    (145, 2),
    (146, -2),
    (147, 0),
    (148, -4),
    (149, 2),
    (150, -2),
    (151, -1),
    (152, 0),
    (162, -2),
    (163, 1),
    (164, 1),
    (165, -2),
    (166, -3),
    (167, 0),
    (168, 0),
    (169, -1),
    (170, 0),
    (180, 4),
    (181, -1),
    (182, -3),
    (183, -2),
    (184, -2),
    (185, -3),
    (186, 0),
    (187, -1),
    (188, -1),
    (189, -3),
    (200, 2),
    (201, 0),
    (202, -3),
    (203, 0),
    (204, -2),
    (205, -2),
    (206, -1),
    (207, -3),
    (208, -2),
    (209, -2),
    (220, 5),
    (221, -3),
    (222, -1),
    (223, -3),
    (224, 0),
    (225, 0),
    (226, -2),
    (227, -2),
    (228, -4),
    (229, -2),
    (242, -2),
    (243, -3),
    (244, 1),
    (245, -3),
    (246, -1),
    (247, -2),
    (248, -1),
    (249, -3),
    (250, -3),
    (264, 6),
    (265, -2),
    (266, -2),
    (267, 2),
    (268, -2),
    (269, -2),
    (270, -3),
    (271, 1),
    (272, -2),
    (288, -2),
    (289, 0),
    (290, -2),
    (291, 0),
    (292, -1),
    (293, -3),
    (294, -2),
    (295, 3),
    (312, 7),
    (313, 0),
    (314, -1),
    (315, -2),
    (316, -1),
    (317, 3),
    (318, -2),
    (319, -3),
    (338, -1),
    (339, 0),
    (340, -1),
    (341, -1),
    (342, -2),
    (343, -3),
    (344, 2),
    (364, 5),
    (365, -2),
    (366, 1),
    (367, -1),
    (368, 1),
    (369, -3),
    (370, -1),
    (392, 1),
    (393, -1),
    (394, 0),
    (395, 1),
    (396, -2),
    (397, -2),
    (420, 5),
    (421, 0),
    (422, -1),
    (423, -3),
    (424, -1),
    (425, -1),
    (450, -1),
    (451, -1),
    (452, -2),
    (453, -4),
    (454, -1),
    (480, 4),
    (481, -1),
    (482, -2),
    (483, -4),
    (484, -2),
    (512, 1),
    (513, -3),
    (514, -2),
    (515, -3),
    (544, 5),
    (545, -2),
    (546, -3),
    (547, -1),
    (578, 0),
    (579, -3),
    (580, -2),
    (612, 4),
    (613, -2),
    (614, -2),
    (648, -3),
    (649, -2),
    (684, 11),
    (685, -1),
    (722, 2),
    (760, 7),
];

#[cfg(test)]
mod test {
    //! Table schema keys must be ORDERED and UNIQUE

    use super::*;

    /// General function to validate all tables.
    /// First checks that the table keys are ordered and unique.
    /// Then, verify the table size is the wanted.
    ///
    /// # Arguments
    ///
    /// * `table` - the table to be validated.
    /// * `expected_size` - the number of entries the table must have.
    fn validate_table_format(table: &[(u16, i8)], expected_size: usize) {
        // if the table is ordered and its keys are unique, then for all keys k_{i}, k_{i+1}
        // the table must satisfy k_{i} < k_{i+1}.
        for i in 0..(table.len() - 1) {
            let a = table[i].0;
            let j = i + 1;
            let b = table[j].0;
            assert!(
                a < b,
                "Invalid table. At position ({i}, {j}) found ({a}, {b})",
            )
        }

        let current_size = table.len();
        assert_eq!(
            expected_size,
            table.len(),
            "The table is incomplete. Expect {expected_size} entries and found {current_size}.",

        )
    }

    #[test]
    #[should_panic = "The table is incomplete. Expect 5 entries and found 2."]
    fn bad_table_size() {
        let table: [(u16, i8); 2] = [(5, 10), (8, 42)];
        validate_table_format(&table, 5);
    }

    #[test]
    #[should_panic = "Invalid table. At position (1, 2) found (1, 1)"]
    fn bad_table_keys_order_with_contiguous_repeats() {
        let table: [(u16, i8); 5] = [(0, 12), (1, 21), (1, 30), (5, 10), (8, 42)];
        validate_table_format(&table, 5);
    }

    #[test]
    #[should_panic = "Invalid table. At position (1, 2) found (5, 1)"]
    fn bad_table_keys_order_with_uncontiguous_repeats() {
        let table: [(u16, i8); 4] = [(1, 0), (5, 10), (1, 30), (8, 42)];
        validate_table_format(&table, 4);
    }

    #[test]
    #[should_panic = "Invalid table. At position (2, 3) found (8, 5)"]
    fn bad_table_keys_order_without_repeats() {
        let table: [(u16, i8); 4] = [(1, 0), (2, 7), (8, 10), (5, 42)];
        validate_table_format(&table, 4);
    }

    /// General function to validate all Amino acid scoring tables.
    ///
    /// # Arguments
    ///
    /// * 'tabe_id' - The identifier of the table. For example, for BLOSUM62 the identifier is 62.
    /// * `expected_size` - the number of entries the table must have.
    /// * 'score_cases' - An array of test cases for some values of the matrix.
    fn check_aminoacid_tables<T>(
        table_id: u8,
        expected_size: usize,
        score_cases: &[(i8, &Aac, &Aac)],
    ) where
        T: AaSymTable + AaSchema + HasName,
    {
        let name = T::name();
        let table = T::get_table(table_id);
        validate_table_format(table, expected_size);

        for (expected_score, aa_code1, aa_code2) in score_cases {
            let current_score = T::get_score(table_id, aa_code1, aa_code2);
            assert_eq!(
                *expected_score, *current_score,
                "Failed {} ID {} for ({:?}, {:?}): expected {} and found {}.",
                name, table_id, aa_code1, aa_code2, expected_score, current_score
            );
        }
    }

    #[test]
    #[should_panic = "Index not found. The scoring schema may be incomplete."]
    fn bad_table_missing_entry() {
        // (F,R) maps to 204
        let fake_table: &'static [(u16, i8); 4] = &[(0, 1), (31, 2), (200, 3), (760, 4)];
        struct Fake {}
        impl AaSymTable for Fake {}
        Fake::search(fake_table, &Aac::G, &Aac::R);
    }

    #[test]
    #[should_panic = "Failed Fake ID 89 for (G, R): expected 50 and found 42."]
    fn bad_table_mapping() {
        struct Fake {}
        impl AaSymTable for Fake {}
        impl HasName for Fake {
            fn name() -> &'static str {
                "Fake"
            }
        }
        impl AaSchema for Fake {
            fn get_score(id: u8, code_1: &Aac, code_2: &Aac) -> &'static i8 {
                Fake::search(Fake::get_table(id), code_1, code_2)
            }

            fn get_table(id: u8) -> &'static [(u16, i8)] {
                let fake_table: &'static [(u16, i8); 4] =
                    &[(0, -10), (31, 127), (204, 42), (760, -4)];
                match id {
                    89 => fake_table,
                    _ => panic!("The test is bad configured."),
                }
            }
        }

        // The test in the strict sense start here
        let table_id = 89;
        let expected_size = 4;
        let score_cases: [(i8, &Aac, &Aac); 2] =
            [(-4, &Aac::Y, &Aac::Y), (50, &Aac::G, &Aac::R)];

        check_aminoacid_tables::<Fake>(table_id, expected_size, &score_cases)
    }

    #[test]
    fn check_blosum62() {
        let table_id = 62;
        let expected_size = 210;
        let static_table = BLOSUM62;

        let score_cases = [
            // Diagonal extremes and mid
            (4, &Aac::A, &Aac::A),
            (7, &Aac::Y, &Aac::Y),
            (6, &Aac::N, &Aac::N),
            // Check symmetria
            (-3, &Aac::P, &Aac::I),
            (-3, &Aac::I, &Aac::P),
            (-2, &Aac::Q, &Aac::G),
            (-2, &Aac::G, &Aac::Q),
        ];

        // the following lines must be present in all blosum checks
        assert_eq!(static_table, Blosum::get_table(table_id));
        check_aminoacid_tables::<Blosum>(table_id, expected_size, &score_cases);
    }

    #[test]
    #[should_panic = "BLOSUM13 is not implemented."]
    fn blosum_index_doesnt_exist() {
        Blosum::get_table(13);
    }
}
