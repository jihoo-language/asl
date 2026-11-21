mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cubic_s_pair_has_an_exact_certificate_and_cusp_has_the_expected_hilbert_function() {
    assert_eq!(
        stdout(example("rank-ideals")),
        concat!(
            "rvgroebner[Intersection] = [x2^2 + x0*x1, x1*x2 + x0^2, x1^2*x2 + x0*x2^2, x2^4 + x1^3*x2]\n",
            "rvhilbert[Intersection, 0] = 1\n",
            "rvhilbert[Intersection, 1] = 3\n",
            "rvhilbert[Intersection, 2] = 4\n",
            "rvhilbert[Intersection, 1000] = 4\n",
            "rvgroebner[Cuspidal] = [x1^3 + x0^2*x2]\n",
            "rvhilbert[Cuspidal, 1000] = 3000\n",
        )
    );
}

#[test]
fn zero_and_unit_ideals_keep_generator_positions_and_polynomial_types() {
    let program = r"
        rvideal Empty[[]] {};
        rvideal Unit[[]] { zero[RVPoly([], 0)], rvmon[[]]{}, rvmon[[]]{} };
        let one = rvmon[[]]{};
        assert rvInIdeal[Empty](one) == S.lo;
        assert rvNormal[Empty](one) == one;
        assert rvInIdeal[Unit](one) == S.hi;
        assert rvWitness[Unit, 0](one) == zero[RVPoly([], 0)];
        assert rvWitness[Unit, 1](one) == one;
        assert rvWitness[Unit, 2](one) == zero[RVPoly([], 0)];
        rvgroebner[Empty];
        rvgroebner[Unit];
        rvhilbert[Empty, 0];
        rvhilbert[Empty, 10];
        rvhilbert[Unit, 0];
        rvhilbert[Unit, 10];

        rvfitting EmptyMinor[One, [0], 0];
        rvfitting AbsentMinor[One, [0], 1];
        assert rvInIdeal[EmptyMinor](rvmon[[0]]{10}) == S.hi;
        assert rvInIdeal[AbsentMinor](rvmon[[0]]{10}) == S.lo;
        let x = rvmon[[0]]{1};
        rvideal Mixed[[0]] { zero[RVPoly([0], 0)], rvMul(x,x), x, rvMul(x,x) };
        let normal = \f : RVPoly([0], 2) => rvNormal[Mixed](f);
        assert normal(rvMul(x,x)) == zero[RVPoly([0], 2)];
        assert normal(zero[RVPoly([0], 2)]) == zero[RVPoly([0], 2)];
        assert rvWitness[Mixed, 2](rvMul(x,x)) == x;
        assert rvWitness[Mixed, 3](rvMul(x,x)) == zero[RVPoly([0], 0)];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvgroebner[Empty] = []\n",
            "rvgroebner[Unit] = [1]\n",
            "rvhilbert[Empty, 0] = 1\n",
            "rvhilbert[Empty, 10] = 0\n",
            "rvhilbert[Unit, 0] = 0\n",
            "rvhilbert[Unit, 10] = 0\n",
        )
    );
}

#[test]
fn hilbert_counts_remain_exact_beyond_machine_integers() {
    let program = r"
        rvideal Empty[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]] {};
        rvideal Hyperplane[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]] {
            rvmon[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]]{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
        };
        rvhilbert[Empty, 1000000];
        rvhilbert[Hyperplane, 1000000];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvhilbert[Empty, 1000000] = 764808144178764238728047439565706886810089560734756525268145657040109274325001\n",
            "rvhilbert[Hyperplane, 1000000] = 11471950083430212127738795511552930008201220392715457147164977380940925001\n",
        )
    );
}

#[test]
fn invalid_ideal_operations_report_local_diagnostics() {
    let setup = "rvideal J[[0,1]] { rvmon[[0,1]]{2,0} };";
    for (program, reason) in [
        ("rvgroebner[Unknown];", "unknown rank ideal"),
        ("rvideal J[[0,1]] {};", "already defined"),
        (
            "rvideal K[[0]] { rvmon[[1]]{1} };",
            "different primitive selection",
        ),
        ("rvideal K[[0]] { S.hi };", "requires an RVPoly"),
        (
            "show rvNormal[J](rvmon[[0]]{1});",
            "different primitive selection",
        ),
        (
            "show rvWitness[J,1](rvmon[[0,1]]{2,0});",
            "index 1 is out of range",
        ),
        ("show rvWitness[J,0](rvmon[[0,1]]{1,0});", "degree greater"),
        (
            "show rvWitness[J,0](rvmon[[0,1]]{0,2});",
            "polynomial in the ideal",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{program}")));
        assert!(error.contains(reason), "{program}\n{error}");
    }
}
