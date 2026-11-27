mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn genuine_noncommuting_ordered_ratios_certify_weil_factors() {
    let text = stdout(example("hyperelliptic-hasse-witt"));
    for expected in [
        "rvHWStatus[C] = genus=2, ordinary=yes",
        "rvHWPoints[C, 1] = 4, trace=0",
        "rvHWPoints[C, 2] = 6, trace=4",
        "rvHWWeil[C] = [1, 0, -2, 0, 9]",
        "RVHWMatrix(3, 2, 3){[[16, 4], [13, 11]]}",
        "rvHWConverges[C, 3] mod 9 = yes",
        "rvHWCertify[C, 3] mod 27 = unit=[1, 0, 16], complement=[1, 0, 9], remainder=[0, 0]",
        "rvHWCertify[C, 7] mod 2187 = unit=[1, 0, 1015], complement=[1, 0, 1170], remainder=[0, 0]",
        "rvHWPoints[G, 3] = 22, trace=6",
        "rvHWWeil[G] = [1, 0, 2, -2, 6, 0, 27]",
        "rvHWCertify[G, 3] mod 27 = unit=[1, 12, 11, 22]",
        "rvHWWeil[E] = [1, -4, 7]",
        "ratio=10, unit=39, trace=4",
    ] {
        assert!(text.contains(expected), "{expected}\n{text}");
    }
}

#[test]
fn actual_high_precision_truncation_and_a_different_integral_lift() {
    let text = stdout(interpret(
        r"
        rvhyperelliptic C[3]{0,1,1,1,2,1};
        rvhyperelliptic D[3]{3,1,1,1,2,1};
        assert rvHW[D,3,3] == rvHWMat[3,2,3]{11,2,20,4};
        assert rvHW[C,3,3] == rvHWMat[3,2,3]{11,2,11,22};
        assert rvHWReduce[1](rvHW[D,3,3]) == rvHW[C,3,1];
        assert rvHW[C,6,6] == rvHWMat[3,2,6]{449,159,528,437};
        assert rvHW[C,5,6] == rvHWMat[3,2,6]{148,589,148,161};
        assert rvHWUnit[C,6] == rvHWMat[3,2,6]{691,706,202,38};
        rvHWConverges[C,6];
        rvHWCertify[C,6];
        rvHWCertify[D,3];
        rvHWWeil[D];
        ",
    ));
    assert!(text.contains("rvHWConverges[C, 6] mod 243 = yes"));
    assert!(text.contains("rvHWCertify[C, 6] mod 729 = unit=["));
    assert!(text.contains("rvHWCertify[D, 3] mod 27 = unit=[1, 0, 16]"));
    assert!(text.ends_with("rvHWWeil[D] = [1, 0, -2, 0, 9]\n"));
}

#[test]
fn smooth_nonordinary_curves_still_have_actual_raw_coefficients() {
    let source = r"
        rvhyperelliptic S[3]{0,-1,0,1};
        rvHWStatus[S];
        assert rvHW[S,1,1] == zero[RVHWMatrix(3,1,1)];
        assert rvHWEntry[1,1](rvHW[S,2,2]) == rvHWScalar[3,2,6];
        rvHWVal(rvHWEntry[1,1](rvHW[S,2,2]));
        rvHWVal(rvHWEntry[1,1](rvHW[S,1,1]));
        rvHWWeil[S];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvHWStatus[S] = genus=1, ordinary=no\n",
            "rvHWVal = 1\n",
            "rvHWVal = >=1\n",
            "rvHWWeil[S] = [1, 0, 3]\n",
        )
    );
    for call in ["show rvHWUnit[S,1];", "rvHWCertify[S,2];"] {
        let result = diagnostic(interpret(&format!("{source}{call}")));
        assert!(result.contains("smooth nonordinary"), "{result}");
    }
}

#[test]
fn characteristic_coefficients_and_matrix_arithmetic_support_the_residue_ring() {
    let text = stdout(interpret(
        r"
        let a = rvHWMat[3,3,2]{1,0,0,0,1,0,0,0,1};
        assert rvHWChar[1](a) == rvHWScalar[3,2,-3];
        assert rvHWChar[2](a) == rvHWScalar[3,2,3];
        assert rvHWChar[3](a) == rvHWScalar[3,2,-1];
        assert rvHWCompose(rvHWInv(a),a) == a;
        assert add(neg(a),a) == zero[RVHWMatrix(3,3,2)];
        assert rvHWReduce[1](rvHWScalar[7,3,-1]) == rvHWScalar[7,1,6];
        let v = rvHWMat[4294967291,1,2]{-1};
        assert rvHWCompose(v,rvHWInv(v)) == rvHWMat[4294967291,1,2]{1};
        let big = rvHWScalar[4294967291,2,-1];
        assert add(rvHWInv(big),rvHWScalar[4294967291,2,1]) == zero[RVHWScalar(4294967291,2)];
        card (\x : RVHWMatrix(3,2,1) => S.hi);
        card (\x : RVHWScalar(3,2) => S.hi);
        ",
    ));
    assert_eq!(text, "card = 81 / 81\ncard = 9 / 9\n");
}

#[test]
fn explicit_prime_polynomial_and_size_checks_are_distinct() {
    for (input, expected) in [
        ("rvhyperelliptic C[2]{0,1};", "odd prime"),
        ("rvhyperelliptic C[9]{0,1};", "not prime"),
        ("rvhyperelliptic C[3]{0,0,0,1};", "repeated roots"),
        ("rvhyperelliptic C[3]{1,1,1};", "monic odd-degree"),
        ("rvhyperelliptic C[3]{1,1,0,4};", "monic odd-degree"),
        ("show rvHWMat[3,2,1]{1,2};", "wrong number of entries"),
        (
            "show rvHWInv(rvHWMat[3,2,2]{3,0,0,1});",
            "nonunit determinant",
        ),
        ("show rvHWInv(rvHWScalar[3,2,6]);", "requires a unit"),
        (
            "show rvHWEntry[0,1](rvHWMat[3,1,1]{1});",
            "entry is out of range",
        ),
        (
            "show rvHWChar[2](rvHWMat[3,1,1]{1});",
            "coefficient is out of range",
        ),
        (
            "show rvHWReduce[3](rvHWScalar[3,2,4]);",
            "no greater precision",
        ),
        ("show rvHWScalar[3,0,1];", "precision must be positive"),
        ("show rvHWScalar[3,81,1];", "modulus p^s exceeds u64"),
        (
            "rvhyperelliptic C[3]{0,1,0,1}; show rvHW[C,80,1];",
            "truncation p^r exceeds u64",
        ),
        (
            "rvhyperelliptic C[3]{0,1,0,1}; rvHWPoints[C,0];",
            "degree must be in",
        ),
    ] {
        let result = diagnostic(interpret(input));
        assert!(result.contains(expected), "{input} -> {result}");
    }
}

#[test]
fn the_linear_equation_gives_the_empty_genus_zero_matrix() {
    let text = stdout(interpret(
        r"
        rvhyperelliptic L[3]{2,1};
        rvHWStatus[L];
        rvHWPoints[L,2];
        rvHWWeil[L];
        assert rvHW[L,0,3] == zero[RVHWMatrix(3,0,3)];
        assert rvHWChar[0](rvHWUnit[L,2]) == rvHWScalar[3,2,1];
        rvHWCertify[L,2];
        ",
    ));
    assert_eq!(
        text,
        concat!(
            "rvHWStatus[L] = genus=0, ordinary=yes\n",
            "rvHWPoints[L, 2] = 10, trace=0\n",
            "rvHWWeil[L] = [1]\n",
            "rvHWCertify[L, 2] mod 9 = unit=[1], complement=[1], remainder=[]\n",
        )
    );
}
