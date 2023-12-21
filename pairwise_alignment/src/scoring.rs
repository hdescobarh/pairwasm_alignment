//! Alignment scoring schemas.

use crate::bioseq::Aac;

/// Allows access to scoring schema name.
pub trait HasName {
    fn name() -> &'static str;
}

/// Types implementing this trait are considered as an **amino acid** scoring schema.
pub trait AaSchema {
    /// Returns the score value sᵢⱼ, with i ≡ code_1 and j ≡ code_2, from the schema with identifier id.
    fn get_score(id: u8, code_1: &Aac, code_2: &Aac) -> &'static i8;
    /// Returns the complete scoring schema with identifier id.
    fn get_table(id: u8) -> &'static [(u16, i8)];
}

/// Amino acid schemas implementing this trait are able to search scores
/// looking only in the upper triangle and main diagonal.
/// Schema keys must be ordered or the binary search will be meaningless.
///
/// The logic is that if Aₙₘ is a symmetric matrix, then aᵢⱼ = aⱼᵢ. So,
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
            45 => BLOSUM45,
            _ => panic!("BLOSUM{id} is not implemented."),
        }
    }
}

impl AaSymTable for Blosum {}

pub struct Pam {}

impl HasName for Pam {
    fn name() -> &'static str {
        "PAM"
    }
}

impl AaSchema for Pam {
    fn get_score(id: u8, code_1: &Aac, code_2: &Aac) -> &'static i8 {
        Pam::search(Pam::get_table(id), code_1, code_2)
    }

    fn get_table(id: u8) -> &'static [(u16, i8)] {
        match id {
            160 => PAM160,
            _ => panic!("PAM{id} is not implemented."),
        }
    }
}

impl AaSymTable for Pam {}

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

static BLOSUM45: &[(u16, i8)] = &[
    (0, 5),
    (2, -1),
    (4, 12),
    (5, -2),
    (8, -3),
    (9, -1),
    (12, 7),
    (13, -3),
    (14, -2),
    (18, 2),
    (19, -2),
    (20, 0),
    (24, 6),
    (25, -4),
    (26, -3),
    (27, -2),
    (32, -3),
    (33, -1),
    (34, -3),
    (35, -1),
    (40, 8),
    (41, -2),
    (42, 0),
    (43, -3),
    (44, -1),
    (50, -3),
    (51, 0),
    (52, -4),
    (53, -3),
    (54, -1),
    (60, 7),
    (61, -2),
    (62, -3),
    (63, 0),
    (64, -2),
    (65, -1),
    (72, -2),
    (73, 0),
    (74, 1),
    (75, -3),
    (76, -2),
    (77, -1),
    (84, 10),
    (85, -4),
    (86, -3),
    (87, -2),
    (88, -3),
    (89, -2),
    (90, -1),
    (98, -3),
    (99, -2),
    (100, 1),
    (101, -2),
    (102, 2),
    (103, -4),
    (104, -1),
    (112, 5),
    (113, -1),
    (114, -3),
    (115, 0),
    (116, 0),
    (117, -1),
    (118, -3),
    (119, -2),
    (128, -3),
    (129, -2),
    (130, -2),
    (131, -2),
    (132, 0),
    (133, 0),
    (134, -3),
    (135, 1),
    (144, 5),
    (145, 2),
    (146, 0),
    (147, 0),
    (148, -3),
    (149, 2),
    (150, -1),
    (151, -1),
    (152, 0),
    (162, -3),
    (163, 2),
    (164, 1),
    (165, -2),
    (166, -4),
    (167, 0),
    (168, 0),
    (169, -1),
    (170, 0),
    (180, 5),
    (181, -1),
    (182, -2),
    (183, -2),
    (184, -2),
    (185, -2),
    (186, 0),
    (187, -1),
    (188, -1),
    (189, -2),
    (200, 2),
    (201, 0),
    (202, -2),
    (203, 1),
    (204, -2),
    (205, -2),
    (206, -1),
    (207, -3),
    (208, -5),
    (209, -2),
    (220, 6),
    (221, -3),
    (222, -1),
    (223, -2),
    (224, 0),
    (225, 0),
    (226, -1),
    (227, -3),
    (228, -4),
    (229, -3),
    (242, -2),
    (243, -3),
    (244, 1),
    (245, -3),
    (246, -1),
    (247, -2),
    (248, 0),
    (249, -3),
    (250, -2),
    (264, 6),
    (265, -2),
    (266, -2),
    (267, 3),
    (268, -2),
    (269, -2),
    (270, -3),
    (271, 1),
    (272, -2),
    (288, -2),
    (289, 0),
    (290, -2),
    (291, -1),
    (292, -1),
    (293, -3),
    (294, -2),
    (295, 3),
    (312, 9),
    (313, 0),
    (314, -1),
    (315, -3),
    (316, -1),
    (317, 3),
    (318, -3),
    (319, -3),
    (338, -1),
    (339, 0),
    (340, -2),
    (341, -1),
    (342, -2),
    (343, -2),
    (344, 2),
    (364, 6),
    (365, -2),
    (366, 1),
    (367, -1),
    (368, 1),
    (369, -2),
    (370, 0),
    (392, 1),
    (393, -1),
    (394, 0),
    (395, 1),
    (396, -2),
    (397, -1),
    (420, 7),
    (421, 0),
    (422, -1),
    (423, -3),
    (424, -2),
    (425, 0),
    (450, -1),
    (451, -1),
    (452, -3),
    (453, -4),
    (454, 0),
    (480, 4),
    (481, -1),
    (482, -3),
    (483, -3),
    (484, -2),
    (512, 2),
    (513, -2),
    (514, -2),
    (515, -3),
    (544, 5),
    (545, -1),
    (546, -2),
    (547, -1),
    (578, 0),
    (579, -4),
    (580, -1),
    (612, 5),
    (613, -3),
    (614, -2),
    (648, -3),
    (649, -1),
    (684, 15),
    (685, -1),
    (722, 3),
    (760, 8),
];

static PAM160: &[(u16, i8)] = &[
    (0, 2),
    (2, -2),
    (4, 9),
    (5, 0),
    (8, -5),
    (9, 0),
    (12, 4),
    (13, -5),
    (14, -3),
    (18, 3),
    (19, -5),
    (20, 1),
    (24, 4),
    (25, -6),
    (26, -3),
    (27, -2),
    (32, -5),
    (33, 0),
    (34, -3),
    (35, -1),
    (40, 7),
    (41, 0),
    (42, 0),
    (43, -2),
    (44, -2),
    (50, -4),
    (51, 0),
    (52, -3),
    (53, -5),
    (54, -2),
    (60, 4),
    (61, -2),
    (62, -2),
    (63, 0),
    (64, -6),
    (65, -1),
    (72, -3),
    (73, 0),
    (74, -1),
    (75, -4),
    (76, -5),
    (77, 0),
    (84, 6),
    (85, -3),
    (86, -5),
    (87, -3),
    (88, -3),
    (89, -4),
    (90, 1),
    (98, -3),
    (99, -2),
    (100, 1),
    (101, -2),
    (102, 2),
    (103, -3),
    (104, -1),
    (112, 5),
    (113, -1),
    (114, -4),
    (115, 0),
    (116, 1),
    (117, -2),
    (118, -5),
    (119, -2),
    (128, -2),
    (129, -2),
    (130, -3),
    (131, -3),
    (132, -1),
    (133, 1),
    (134, -3),
    (135, 1),
    (144, 4),
    (145, 2),
    (146, -3),
    (147, 0),
    (148, -4),
    (149, 2),
    (150, -2),
    (151, 0),
    (152, 1),
    (162, -3),
    (163, 2),
    (164, 2),
    (165, -1),
    (166, -5),
    (167, -2),
    (168, 0),
    (169, -2),
    (170, 0),
    (180, 5),
    (181, 0),
    (182, -2),
    (183, -1),
    (184, -2),
    (185, -4),
    (186, 0),
    (187, -1),
    (188, -2),
    (189, -5),
    (200, 3),
    (201, 1),
    (202, -2),
    (203, 2),
    (204, -3),
    (205, -3),
    (206, -1),
    (207, -3),
    (208, -7),
    (209, -3),
    (220, 7),
    (221, -3),
    (222, -2),
    (223, -2),
    (224, 1),
    (225, 1),
    (226, -3),
    (227, -2),
    (228, -6),
    (229, 0),
    (242, -2),
    (243, -3),
    (244, 0),
    (245, -2),
    (246, -1),
    (247, -1),
    (248, -2),
    (249, -7),
    (250, -4),
    (264, 3),
    (265, -2),
    (266, -2),
    (267, 3),
    (268, -2),
    (269, -2),
    (270, -2),
    (271, -1),
    (272, -4),
    (288, -1),
    (289, -1),
    (290, -3),
    (291, -1),
    (292, 0),
    (293, -2),
    (294, -7),
    (295, 5),
    (312, 5),
    (313, 0),
    (314, -1),
    (315, -3),
    (316, 0),
    (317, 3),
    (318, -3),
    (319, -5),
    (338, 0),
    (339, -1),
    (340, -2),
    (341, -2),
    (342, -3),
    (343, -5),
    (344, 0),
    (364, 5),
    (365, -1),
    (366, 1),
    (367, -1),
    (368, 1),
    (369, -4),
    (370, -2),
    (392, 1),
    (393, 1),
    (394, 0),
    (395, 1),
    (396, -2),
    (397, -4),
    (420, 6),
    (421, -1),
    (422, 0),
    (423, -2),
    (424, -4),
    (425, -2),
    (450, -1),
    (451, -1),
    (452, -2),
    (453, -4),
    (454, -3),
    (480, 2),
    (481, -1),
    (482, -2),
    (483, -5),
    (484, -2),
    (512, 1),
    (513, -3),
    (514, -5),
    (515, -5),
    (544, 3),
    (545, -1),
    (546, 1),
    (547, -4),
    (578, 0),
    (579, -2),
    (580, -4),
    (612, 4),
    (613, -6),
    (614, -3),
    (648, -6),
    (649, -3),
    (684, 12),
    (685, -3),
    (722, -1),
    (760, 8),
];

#[cfg(test)]
mod test {
    //! Schema keys must be ORDERED and UNIQUE.

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
        // if the table is ordered and its keys are unique, then for all keys kᵢ, kᵢ₊₁
        // the table must satisfy kᵢ < kᵢ₊₁.
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

    /// General function to validate all amino acid scoring tables.
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
    #[should_panic = "BLOSUM13 is not implemented."]
    fn blosum_index_doesnt_exist() {
        Blosum::get_table(13);
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
    fn check_blosum45() {
        let table_id = 45;
        let expected_size = 210;
        let static_table = BLOSUM45;

        let score_cases = [
            // Diagonal extremes and mid
            (5, &Aac::A, &Aac::A),
            (8, &Aac::Y, &Aac::Y),
            (6, &Aac::N, &Aac::N),
            // Check symmetria
            (-2, &Aac::P, &Aac::I),
            (-2, &Aac::I, &Aac::P),
            (-2, &Aac::Q, &Aac::G),
            (-2, &Aac::G, &Aac::Q),
        ];

        // the following lines must be present in all blosum checks
        assert_eq!(static_table, Blosum::get_table(table_id));
        check_aminoacid_tables::<Blosum>(table_id, expected_size, &score_cases);
    }
}
