mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exceptional_units_and_commuting_conductor_diamond() {
    assert_eq!(
        stdout(example("ring-class-polynomial-norms")),
        concat!(
            "rvCMClassPolyData[H48] = discriminant=-48, degree=2\n",
            "rvCMClassPolyCoeff[H48,0] = 6549518250000\n",
            "rvCMClassPolyCoeff[H48,1] = -2835810000\n",
            "rvCMClassPolyCoeff[H16,0] = -287496\n",
            "rvCMClassPolyCoeff[H27,0] = 12288000\n",
            "rvCMClassPolyData[H108b] = discriminant=-108, degree=3\n",
            "rvCMClassPolyCoeff[H108a,0] = -1879994705688000000000\n",
            "rvCMClassPolyCoeff[H108b,0] = -1879994705688000000000\n",
            "rvCMClassPolyData[H432] = discriminant=-432, degree=6\n",
        )
    );
}

#[test]
fn norms_of_quadratic_and_cubic_parents_match_independent_cm_products() {
    let mut source = "rvcmclasspoly H3[-3];
        rvcmclasspolylift H12[H3,2];
        rvcmclasspolylift H48[H12,2];
        rvcmclasspolylift H192[H48,2];
        rvcmclasspolylift H27[H3,3];
        rvcmclasspolylift H108[H12,3];
        rvcmclasspolylift H432[H108,2];"
        .to_string();
    // Independent products of j(tau) over the reduced forms, computed
    // with the analytic q-product at 220 decimal digits.
    let expected: [(&str, &[&str]); 2] = [
        (
            "H192",
            &[
                "-1080060886113159937649308593750000",
                "826335556188178615474500000000",
                "15705521635909735050750000",
                "-8041801037378436000",
                "1",
            ],
        ),
        (
            "H432",
            &[
                "42889619864187195342544128412237640625000000000000",
                "3869372376492639837782614434923625000000000000",
                "34904627315764077727184412247908187500000000",
                "1007059405271040783775694468925000000000",
                "280179539493990596285512318134750000",
                "-22804995243537595825782822000",
                "1",
            ],
        ),
    ];
    let mut expected_output = String::new();
    for (name, coefficients) in expected {
        for (index, coefficient) in coefficients.iter().enumerate() {
            writeln!(source, "rvCMClassPolyCoeff[{name},{index}];").unwrap();
            writeln!(
                expected_output,
                "rvCMClassPolyCoeff[{name},{index}] = {coefficient}"
            )
            .unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), expected_output);
}

#[test]
fn all_nine_fundamental_seeds_are_available_unconditionally() {
    let mut source = String::new();
    let mut expected_output = String::new();
    for (discriminant, constant) in [
        (-3, 0),
        (-4, -1728),
        (-7, 3375),
        (-8, -8000),
        (-11, 32768),
        (-19, 884736),
        (-43, 884736000),
        (-67, 147197952000),
        (-163, 262537412640768000i128),
    ] {
        let name = format!("H{}", -discriminant);
        writeln!(
            source,
            "rvcmclasspoly {name}[{discriminant}]; rvCMClassPolyCoeff[{name},0];"
        )
        .unwrap();
        writeln!(expected_output, "rvCMClassPolyCoeff[{name},0] = {constant}").unwrap();
    }
    assert_eq!(stdout(interpret(&source)), expected_output);
}

#[test]
fn genuinely_required_ascent_and_invalid_inputs() {
    let setup = "rvcmclasspoly H3[-3]; rvcmclasspolylift H12[H3,2];
        rvcmclasspolylift H108[H12,3]; rvunram F[5,2]{0,1};";
    for (command, expected) in [
        (
            "rvcmclasspoly Bad[-12];",
            "requires a class-number-one fundamental discriminant",
        ),
        ("rvcmclasspoly H3[-4];", "is already defined"),
        (
            "rvcmclasspolylift Missing[Unknown,2];",
            "unknown CM class polynomial",
        ),
        (
            "rvcmclasspolylift Missing[H108,2];",
            "ascending CM class polynomial for discriminant -27 must be registered",
        ),
        ("rvcmclasspolylift Bad[H3,9];", "9 is not prime"),
        (
            "rvcmclasspolylift Bad[H3,11];",
            "exact classical modular polynomial exceeds the computation limit",
        ),
        (
            "show rvCMClassPolyEval[H3](rvUnramElem[F]{1});",
            "requires a precision-one RVUnram field",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}
