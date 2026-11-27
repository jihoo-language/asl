mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_sharp_bound_is_shared_by_curves_with_four_different_first_slopes() {
    assert_eq!(
        stdout(example("rank-eo-first-slope")),
        concat!(
            "EOType(8){0, 1, 1, 2, 2, 3, 3, 4}\n",
            "eoFirstSlope = 1/4\n",
            "eoSlopeCycles = 1\n",
            "rvASNewton[Pure] = [1/2: 16]\n",
            "rvASNewton[WithSeven] = [1/3: 3, 1/2: 10, 2/3: 3]\n",
            "rvASNewton[WithEleven] = [3/8: 8, 5/8: 8]\n",
            "rvASNewton[WithFifteen] = [1/4: 4, 1/2: 8, 3/4: 4]\n",
        )
    );
}

#[test]
fn full_final_flags_have_the_periodic_cycle_counts() {
    let program = r"
        let ordinary = eoType{1,2,3};
        let superspecial = eoType{0,0};
        let two_fifths = eoType{0,0,1,1,2};
        let one_third = eoType{0,0,1,2,2,3};
        eoFirstSlope(ordinary); eoSlopeCycles(ordinary);
        eoFirstSlope(superspecial); eoSlopeCycles(superspecial);
        eoFirstSlope(two_fifths); eoSlopeCycles(two_fifths);
        eoFirstSlope(one_third); eoSlopeCycles(one_third);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "eoFirstSlope = 0\n",
            "eoSlopeCycles = 3\n",
            "eoFirstSlope = 1/2\n",
            "eoSlopeCycles = 2\n",
            "eoFirstSlope = 2/5\n",
            "eoSlopeCycles = 1\n",
            "eoFirstSlope = 1/3\n",
            "eoSlopeCycles = 2\n",
        )
    );
}

#[test]
fn slope_queries_require_a_nonempty_final_type() {
    for operation in ["eoFirstSlope", "eoSlopeCycles"] {
        let error = diagnostic(interpret(&format!("{operation}(eoType{{}});")));
        assert!(
            error.contains("first-slope invariants require positive genus"),
            "{error}"
        );
        let error = diagnostic(interpret(&format!("{operation}(S.lo);")));
        assert!(error.contains("requires an EO final type"), "{error}");
    }
}
