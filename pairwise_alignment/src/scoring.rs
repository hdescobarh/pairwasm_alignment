//! Alignment scoring schemas

use crate::bioseq::Aac;

/// Schemas that are symmetric only needs the upper triangular matrix
/// Schema keys must be ordered or the binary search will be meaningless
///
/// The logic is that if A_{n,m} is a symmetric matrix, then a_{i,j} = a_{j, i}. So,
/// We only need to store the values in the diagonal plus values above or bellow it.
/// The choice of the upper triangle is arbitrary.
pub fn read_symmetric_protein_schema(
    schema: &[(u16, i8)],
    code_1: &Aac,
    code_2: &Aac,
) -> i8 {
    let mut param = [code_1, code_2];
    param.sort();
    let key: u16 = Aac::duple_pairing(param[0], param[1]);
    match schema.binary_search_by(|(k, _)| k.cmp(&key)) {
        Ok(index) => schema[index].1,
        Err(_) => panic!("Index not found. The scoring schema may be incomplete."),
    }
}

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

    /// Checks the table keys are ordered and unique
    fn validate_table(table: &[(u16, i8)]) {
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
    }

    #[test]
    fn check_blosum62() {
        let schema = BLOSUM62;

        // ensure the lenght is the correct
        let expected_entries = 210;
        assert_eq!(
            expected_entries,
            schema.len(),
            "The table is incomplete. Expect {} entries and found {}.",
            expected_entries,
            schema.len()
        );

        // table ordering and uniqueness validation
        validate_table(schema);

        // Diagonal extremes and mid
        assert_eq!(4, read_symmetric_protein_schema(schema, &Aac::A, &Aac::A));
        assert_eq!(7, read_symmetric_protein_schema(schema, &Aac::Y, &Aac::Y));
        assert_eq!(6, read_symmetric_protein_schema(schema, &Aac::N, &Aac::N));
        // Check symmetria
        assert_eq!(-3, read_symmetric_protein_schema(schema, &Aac::P, &Aac::I));
        assert_eq!(-3, read_symmetric_protein_schema(schema, &Aac::I, &Aac::P));
        assert_eq!(-2, read_symmetric_protein_schema(schema, &Aac::Q, &Aac::G));
        assert_eq!(-2, read_symmetric_protein_schema(schema, &Aac::G, &Aac::Q));
    }
}
