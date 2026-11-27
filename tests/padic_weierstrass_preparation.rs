mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_preparation_unramified_roots_and_ramified_slopes() {
    let output = stdout(example("padic-weierstrass-preparation"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        &lines[..8],
        [
            "rvTatePrepare = coordinate=X, mu=0, precision=4, modulus=81, degree=2, polynomial=[55, 72, 1], unit=[28, 12]",
            "rvTateFactors = coordinate=X, mu=0, precision=4, modulus=81, factors=[([1, 0, 1], 1, [55, 72, 1])]",
            "rvTateUnramified = coordinate=X, mu=0, precision=4, modulus=81, factor=0, defining_polynomial=[1, 0, 1], roots=[[45, 28], [45, 53]]",
            "rvTatePrepare = coordinate=X, mu=0, precision=4, modulus=81, degree=2, polynomial=[78, 0, 1], unit=[1, 3]",
            "rvTateNewton = coordinate=X, mu=0, precision=4, modulus=81, valuations=[(1, 2, 2)], unresolved=None",
            "rvTateUnramified = coordinate=X, mu=0, precision=4, modulus=81, factor=0, residual=[0, 1], multiplicity=2, unresolved",
            "rvTatePrepare = coordinate=RVMazurDisc(O,3,0,0,0,1), mu=1, precision=2, modulus=121, degree=1, polynomial=[6, 1], unit=[24, 110]",
            "rvTateUnramified = coordinate=RVMazurDisc(O,3,0,0,0,1), mu=1, precision=2, modulus=121, factor=0, defining_polynomial=[6, 1], roots=[[115]]",
        ]
    );
    assert_eq!(
        lines[8],
        "rvTatePrepare = coordinate=RVMazurDisc(O,3,0,0,55,2), mu=2, precision=1, modulus=11, degree=1, polynomial=[1, 1], unit=[2]",
    );
}

#[test]
fn typed_series_and_coprime_primary_lifts_preserve_unit_content() {
    let program = "
        let first = rvTatePoly[3,5](4,1);
        let repeated = rvTatePoly[3,5](3,0,1);
        let unramified = rvTatePoly[3,5](1,3,1);
        let p = rvTateMul(first,rvTateMul(repeated,unramified));
        let unit = rvTatePoly[3,5](1,3,0,9);
        let h = rvTateMul(rvTatePoly[3,5](3),rvTateMul(unit,p));
        assert add(h,neg(h)) == zero[RVTate(3,5)];
        assert p == rvTateMul(rvTateMul(first,unramified),repeated);
        rvTatePrepare(h);
        rvTateFactors(h);
        rvTateNewton(h);
        rvTateCoeff[4000](h);
        rvTateUnramified[0](h);
        rvTateUnramified[1](h);
        rvTatePrepare(rvTatePoly[3,5](486,-243));";
    let output = stdout(interpret(program));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].ends_with("unit=[1, 3, 0, 9]"), "{output}");
    assert_eq!(
        &lines[1..],
        [
            "rvTateFactors = coordinate=X, mu=1, precision=4, modulus=81, factors=[([0, 1], 2, [3, 0, 1]), ([1, 1], 1, [4, 1]), ([1, 0, 1], 1, [1, 3, 1])]",
            "rvTateNewton = coordinate=X, mu=1, precision=4, modulus=81, valuations=[(0, 1, 3), (1, 2, 2)], unresolved=None",
            "rvTateCoeff[4000] = coordinate=X, precision=5, modulus=243, residue=0",
            "rvTateUnramified = coordinate=X, mu=1, precision=4, modulus=81, factor=0, residual=[0, 1], multiplicity=2, unresolved",
            "rvTateUnramified = coordinate=X, mu=1, precision=4, modulus=81, factor=1, defining_polynomial=[1, 1], roots=[[77]]",
            "rvTatePrepare = coordinate=X, precision=5, unresolved",
        ]
    );
}

#[test]
fn cubic_over_f2_and_valuation_floors_agree_with_independent_oracles() {
    // The root coordinates and Frobenius order were obtained by enumerating
    // the 32^3 coefficient triples in the displayed cubic quotient.
    let program = "
        let cubic = rvTatePoly[2,5](5,3,4,1);
        rvTateUnramified[0](cubic);
        rvTateNewton(rvTatePoly[3,3](0,0,0,0,0,9,1));
        rvTateNewton(rvTatePoly[3,3](0,3,3,1));
        rvTateNewton(rvTatePoly[3,3](0,0,3,1));
        rvTateNewton(rvTatePoly[3,4](27,9,3,1));
        rvTatePrepare(rvTatePoly[2,5](2,4,24));";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvTateUnramified = coordinate=X, mu=0, precision=5, modulus=32, factor=0, defining_polynomial=[1, 1, 0, 1], roots=[[16, 9, 10], [22, 2, 19], [22, 21, 3]]\n",
            "rvTateNewton = coordinate=X, mu=0, precision=3, modulus=27, valuations=[], unresolved=Some((1, 2, 6))\n",
            "rvTateNewton = coordinate=X, mu=0, precision=3, modulus=27, valuations=[(1, 2, 2)], unresolved=Some((2, 1, 1))\n",
            "rvTateNewton = coordinate=X, mu=0, precision=3, modulus=27, valuations=[], unresolved=Some((1, 1, 3))\n",
            "rvTateNewton = coordinate=X, mu=0, precision=4, modulus=81, valuations=[(1, 1, 3)], unresolved=None\n",
            "rvTatePrepare = coordinate=X, mu=1, precision=4, modulus=16, degree=0, polynomial=[1], unit=[1, 2, 12]\n",
        )
    );
}

#[test]
fn malformed_grades_factor_indices_and_other_types_are_diagnosed() {
    for (source, expected) in [
        ("show rvTatePoly[4,2](1);", "characteristic must be prime"),
        ("show rvTatePoly[3,0](1);", "precision must be positive"),
        ("show rvTatePoly[3,100](1);", "modulus p^s exceeds u64"),
        (
            "show rvTateMul(rvTatePoly[3,4](0,1),rvTatePoly[3,2](1));",
            "RVTate(3,4)",
        ),
        (
            "rvTatePrepare(S.hi);",
            "an RVTate(p,s) series or an Amice weight disc",
        ),
        (
            "rvTateUnramified[1](rvTatePoly[3,4](1,0,1));",
            "primary factor index 1 is out of range",
        ),
    ] {
        let message = diagnostic(interpret(source));
        assert!(message.contains(expected), "{message}");
    }
}
