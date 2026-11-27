mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn curves_with_the_same_zeta_numerator_have_incomparable_final_types() {
    assert_eq!(
        stdout(example("rank-as-odd-eo-closure")),
        concat!(
            "rvASOddLPolynomial[First] = [1, 3, 3, 9, 18, 18, 63, 81, 108, 486, 648, 1215, 3483, 3645, 5832, 13122, 8748, 19683, 45927, 39366, 118098, 177147, 177147, 531441, 531441]\n",
            "rvASOddLPolynomial[Second] = [1, 3, 3, 9, 18, 18, 63, 81, 108, 486, 648, 1215, 3483, 3645, 5832, 13122, 8748, 19683, 45927, 39366, 118098, 177147, 177147, 531441, 531441]\n",
            "EOType(12){0, 0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}\n",
            "EOType(12){0, 0, 0, 1, 2, 2, 3, 4, 4, 4, 5, 5}\n",
            "eoDim = 29\n",
            "eoDim = 30\n",
            "eoClosureWitness = none\n",
            "eoClosureWitness = none\n",
        )
    );
}

#[test]
fn odd_bridge_returns_the_same_nominal_type() {
    let source = r"
        rvasodd Odd[5] {(1,1)};
        assert rvASOddEOType[Odd] == eoType{};
        assert eoClosure(rvASOddEOType[Odd],eoType{}) == S.hi;
        eoWeyl(rvASOddEOType[Odd]);
    ";
    assert_eq!(stdout(interpret(source)), "eoWeyl = []\n");
    let error = diagnostic(interpret("show rvASOddEOType[Unknown];"));
    assert!(
        error.contains("unknown odd Artin–Schreier cover"),
        "{error}"
    );
}
