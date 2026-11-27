mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn oriented_prime_cycles_match_the_f49_canonical_frobenius() {
    assert_eq!(
        stdout(example("cm-artin-reciprocity")),
        concat!(
            "rvCMArtinPrime[O,2,0] = form=[2, 1, 3], class=[2, 1, 3], order=3\n",
            "rvCMArtinOrbit[O,2,0] = period=3, classes=[[1, 1, 6], [2, 1, 3], [2, -1, 3]]\n",
            "rvCMArtinOrbit[O,2,1] = period=3, classes=[[1, 1, 6], [2, -1, 3], [2, 1, 3]]\n",
            "rvCMArtinFactors[H48,F] = p=7, residue_degree=1, artin_order=2, factor_degree=2, factor_count=1, factors=[[[5], [5], [1]]]\n",
            "rvCMArtinFactors[H48,L] = p=7, residue_degree=2, artin_order=2, factor_degree=1, factor_count=2, factors=[[[0, 6], [1, 0]], [[5, 1], [1, 0]]]\n",
            "rvCMArtinLift[M,W] = p=7, residue_degree=2, precision=3, artin_order=2, exact_period=2, j=[287, 71]\n",
        )
    );
}

#[test]
fn binary_factorization_splits_an_order_four_cycle_over_f4() {
    let output = stdout(interpret(
        "rvcmclasspoly H7[-7]; rvcmclasspolylift H63[H7,3];
         rvcmorder O[-63]; rvunram F[2,1]{0,1}; rvunram L[2,1]{1,1,1};
         let p=rvCMArtinClass[O,2,0];
         assert p==rvCMClass[O,2,1,8];
         assert rvCMArtinTranslate[O,2,0,2](rvCMPrincipal[O])==rvCMClass[O,4,1,4];
         assert rvCMArtinTranslate[O,2,1,4](p)==p;
         rvCMArtinFactors[H63,F]; rvCMArtinFactors[H63,L];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvCMArtinFactors[H63,F] = p=2, residue_degree=1, artin_order=4, factor_degree=4, factor_count=1, factors=[[[1], [0], [0], [1], [1]]]\n",
            "rvCMArtinFactors[H63,L] = p=2, residue_degree=2, artin_order=4, factor_degree=2, factor_count=2, factors=[[[0, 1], [0, 1], [1, 0]], [[1, 1], [1, 1], [1, 0]]]\n",
        )
    );
}

#[test]
fn odd_prime_orientation_and_actual_cubic_and_split_factors() {
    let output = stdout(interpret(
        "rvcmorder O[-23];
         let a=rvCMClass[O,2,1,3];
         assert rvCMArtinClass[O,3,0]==rvCMInverse(a);
         assert rvCMArtinClass[O,13,0]==a;
         assert rvCMArtinClass[O,13,1]==rvCMInverse(a);
         rvCMArtinPrime[O,13,0];
         rvcmclasspoly H3[-3]; rvcmclasspolylift H12[H3,2];
         rvcmclasspolylift H48[H12,2]; rvcmclasspolylift H108[H12,3];
         rvunram F[7,1]{0,1}; rvunram G[13,1]{0,1};
         rvCMArtinFactors[H108,F]; rvCMArtinFactors[H48,G];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvCMArtinPrime[O,13,0] = form=[13, -9, 2], class=[2, 1, 3], order=3\n",
            "rvCMArtinFactors[H108,F] = p=7, residue_degree=1, artin_order=3, factor_degree=3, factor_count=1, factors=[[[1], [3], [4], [1]]]\n",
            "rvCMArtinFactors[H48,G] = p=13, residue_degree=1, artin_order=1, factor_degree=1, factor_count=2, factors=[[[7], [1]], [[9], [1]]]\n",
        )
    );
}

#[test]
fn prime_forms_and_signed_actions_use_the_full_integer_boundary() {
    let min = i128::MIN;
    let output = stdout(interpret(&format!(
        "rvcmorder Wide[{min}]; rvcmorder O[-23];
         let p=rvCMClass[Wide,3,-2,14178431955039102644307275309657008811];
         assert rvCMArtinClass[Wide,3,0]==p;
         assert rvCMArtinClass[Wide,3,1]==rvCMInverse(p);
         assert rvCMArtinTranslate[O,2,0,{min}](rvCMPrincipal[O])==rvCMClass[O,2,1,3];"
    )));
    assert!(output.is_empty(), "{output}");
}

#[test]
fn unsplit_ramified_conductor_and_type_errors() {
    let setup = "rvcmorder O[-23]; rvcmorder Twin[-23]; rvcmorder R[-48]; rvcmorder I[-3];
        rvcmclasspoly H3[-3]; rvunram W[7,2]{0,1}; rvunram F[5,1]{0,1};";
    for (query, expected) in [
        ("show rvCMArtinClass[I,2,0];", "prime must split"),
        ("show rvCMArtinClass[O,5,0];", "prime must split"),
        (
            "show rvCMArtinClass[O,23,0];",
            "prime must not divide the order discriminant",
        ),
        (
            "show rvCMArtinClass[R,2,0];",
            "prime must not divide the order discriminant",
        ),
        ("show rvCMArtinClass[O,9,0];", "9 is not prime"),
        ("show rvCMArtinClass[O,2,2];", "orientation must be 0 or 1"),
        (
            "show rvCMArtinTranslate[O,2,0,1](rvCMPrincipal[Twin]);",
            "expected type RVCMClass(O)",
        ),
        (
            "rvCMArtinOrbit[O,2,0](rvCMPrincipal[Twin]);",
            "expected type RVCMClass(O)",
        ),
        (
            "rvCMArtinFactors[H3,W];",
            "requires a precision-one RVUnram field",
        ),
        ("rvCMArtinFactors[H3,F];", "prime must split"),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {query}")));
        assert!(error.contains(expected), "{error}");
    }
}
