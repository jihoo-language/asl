mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn same_complete_character_factors_distinct_de_rham_final_types() {
    assert_eq!(
        stdout(example("rank-as-odd-stationary")),
        concat!(
            "rvASOddCharacterL[First, 1] = [[1, 0], [1, -1], [0, 0], [6, 3], [3, -3], [0, 0], [18, 0], [0, 0], [54, 27], [81, -81], [0, 0], [486, 243], [729, 0]]\n",
            "rvASOddLPolynomial[First] = [1, 3, 3, 9, 18, 18, 63, 81, 108, 486, 648, 1215, 3483, 3645, 5832, 13122, 8748, 19683, 45927, 39366, 118098, 177147, 177147, 531441, 531441]\n",
            "rvASOddLPolynomial[Second] = [1, 3, 3, 9, 18, 18, 63, 81, 108, 486, 648, 1215, 3483, 3645, 5832, 13122, 8748, 19683, 45927, 39366, 118098, 177147, 177147, 531441, 531441]\n",
            "rvASOddPoints[First, 12] = 535573\n",
            "rvASOddPoints[Second, 12] = 535573\n",
            "rvASOddaNumber[First] = 7\n",
            "rvASOddaNumber[Second] = 7\n",
            "rvASOddEO[First] = [0, 0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]\n",
            "rvASOddEO[Second] = [0, 0, 0, 1, 2, 2, 3, 4, 4, 4, 5, 5]\n",
        )
    );
}

#[test]
fn signed_constant_twists_at_all_angles_and_coefficient_indices() {
    let mut source = String::from(
        "rvasodd Base[3] {(13,1),(11,2),(10,1),(1,1)};\n\
         rvasodd Twist[3] {(13,1),(11,-1),(10,1),(1,1),(0,-1)};\n",
    );
    for character in 1..3 {
        for k in 0..13 {
            let phase = -(character * k);
            writeln!(source, "assert rvASOddCoeff[Twist,{character},{k}] == rvCycloMul(rvCyclo[3,{phase}],rvASOddCoeff[Base,{character},{k}]);").unwrap();
        }
    }
    source.push_str("rvASOddPoints[Base,12]; rvASOddPoints[Twist,12];");
    assert_eq!(
        stdout(interpret(&source)),
        "rvASOddPoints[Base, 12] = 535573\nrvASOddPoints[Twist, 12] = 535573\n",
    );
}

#[test]
fn exponential_sum_respects_statement_budget() {
    let error = diagnostic(interpret("rvasodd C[5] {(21,1)}; rvASOddCharacterL[C,1];"));
    assert!(
        error.contains("statement exceeds the computation limit"),
        "{error}"
    );
}
