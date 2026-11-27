mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-ks-nomizu-massey.asl")
        .split_once("// Rational cup products")
        .unwrap()
        .0
}

#[test]
fn example_computes_a_closed_nonformal_triple_and_exact_product() {
    assert_eq!(
        stdout(example("hyperelliptic-ks-nomizu-massey")),
        concat!(
            "rvKSNomizuCohomEqual = witness_valid=true\n",
            "rvKSNomizuClosed = true\n",
            "rvKSNomizuCoeff[1,1] = 1/6\n",
            "rvKSNomizuMasseyCertificate = triple=<eta,eta,Df>, degree=2, excludes_zero=true\n",
            "rvKSNomizuData[A,2,2] = cochains=2080, closed=1328, exact=16, cohomology=1312\n",
        )
    );
}

#[test]
fn mixed_rational_cochains_have_graded_cup_and_differential_signs() {
    let source = format!(
        r"{}
        let eta=rvKSNomizuEta[A,2];
        let x=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,1]);
        let y=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,64]);
        let targetx=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,2]);
        let targety=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,128]);
        let unit=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,0]);
        let top=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,18446744073709551615]);
        assert rvKSNomizuCup(unit,x)==x;
        assert rvKSNomizuCup(eta,eta)==zero[RVKSNomizu(A,2,2)];
        assert rvKSNomizuCup(x,eta)==neg(rvKSNomizuCup(eta,x));
        let xy=rvKSNomizuCup(x,y);
        assert rvKSNomizuCup(xy,eta)==rvKSNomizuCup(eta,xy);
        assert rvKSNomizuCup(xy,eta)==rvKSNomizuCup(x,rvKSNomizuCup(y,eta));
        assert rvKSNomizuD(x)==rvKSNomizuScale[6,1](rvKSNomizuCup(eta,targetx));
        assert rvKSNomizuD(y)==rvKSNomizuScale[6,1](rvKSNomizuCup(eta,targety));
        let a=add(rvKSNomizuScale[2,3](x),rvKSNomizuScale[-5,7](eta));
        let b=add(rvKSNomizuScale[-9,11](y),rvKSNomizuScale[1,3](eta));
        let c=add(rvKSNomizuScale[2,5](targetx),rvKSNomizuScale[3,13](targety));
        let ab=rvKSNomizuCup(a,b);
        assert ab==neg(rvKSNomizuCup(b,a));
        assert rvKSNomizuCup(ab,c)==rvKSNomizuCup(a,rvKSNomizuCup(b,c));
        assert rvKSNomizuD(ab)==add(rvKSNomizuCup(rvKSNomizuD(a),b),neg(rvKSNomizuCup(a,rvKSNomizuD(b))));
        let ac=rvKSNomizuCup(a,c);
        assert rvKSNomizuD(rvKSNomizuCup(ac,b))==add(rvKSNomizuCup(rvKSNomizuD(ac),b),rvKSNomizuCup(ac,rvKSNomizuD(b)));
        assert rvKSNomizuD(rvKSNomizuD(ab))==zero[RVKSNomizu(A,2,4)];
        assert rvKSNomizuD(rvKSNomizuCup(eta,x))==zero[RVKSNomizu(A,2,3)];
        let oriented=rvKSNomizuCup(eta,top);
        assert rvKSNomizuCup(top,eta)==oriented;
        assert rvKSNomizuD(oriented)==zero[RVKSNomizu(A,2,66)];
        assert rvKSNomizuCup(oriented,targetx)==zero[RVKSNomizu(A,2,66)];
        assert rvKSNomizuD(zero[RVKSNomizu(A,2,-1)])==zero[RVKSNomizu(A,2,0)];
        let productFiber=rvKSNomizuFiber[A,0](rvKSExteriorBasis[P,1]);
        let productEta=rvKSNomizuEta[A,0];
        let productZero=zero[RVKSNomizu(A,0,1)];
        assert rvKSNomizuD(productFiber)==zero[RVKSNomizu(A,0,2)];
        assert rvKSNomizuCup(productFiber,productEta)==neg(rvKSNomizuCup(productEta,productFiber));
        assert rvKSNomizuMassey(productEta,productEta,productEta,productZero,productZero)==zero[RVKSNomizu(A,0,2)];
        rvKSNomizuCoeff[65,1](rvKSNomizuCup(xy,eta));
        rvKSNomizuCoeff[18446744073709551615,1](oriented);
        rvKSNomizuClosed(a);
        rvKSNomizuClosed(rvKSNomizuD(a));
        ",
        setup()
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvKSNomizuCoeff[65,1] = 1\n",
            "rvKSNomizuCoeff[18446744073709551615,1] = 1\n",
            "rvKSNomizuClosed = false\nrvKSNomizuClosed = true\n",
        )
    );
}

#[test]
fn explicit_coboundaries_and_general_massey_signs_are_checked() {
    let source = format!(
        r"{}
        let eta=rvKSNomizuEta[A,2];
        let x=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,1]);
        let targetx=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,2]);
        let targety=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,128]);
        let inert=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,4]);
        let hx=rvKSNomizuScale[1,6](x);
        let exact=rvKSNomizuCup(eta,targetx);
        let rep=rvKSNomizuCup(eta,hx);
        let other=add(rep,exact);
        rvKSNomizuCohomEqual(other,rep,hx);
        rvKSNomizuCohomEqual(other,rep,zero[RVKSNomizu(A,2,1)]);
        let adjusted=rvKSNomizuMassey(eta,eta,targetx,eta,add(hx,targety));
        assert adjusted==add(add(rep,rvKSNomizuCup(eta,targety)),exact);
        let even=rvKSNomizuCup(eta,inert);
        let closed2=rvKSNomizuCup(inert,targety);
        let higher=rvKSNomizuMassey(even,eta,targetx,closed2,hx);
        assert higher==add(rvKSNomizuCup(even,hx),neg(rvKSNomizuCup(closed2,targetx)));
        rvKSNomizuCoeff[134,0](higher);
        rvKSNomizuClosed(higher);
        ",
        setup()
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvKSNomizuCohomEqual = witness_valid=true\n",
            "rvKSNomizuCohomEqual = witness_valid=false\n",
            "rvKSNomizuCoeff[134,0] = -1\n",
            "rvKSNomizuClosed = true\n",
        )
    );
}

#[test]
fn actual_integral_spin_factors_and_signed_covers_give_exact_certificates() {
    let mut source = setup().to_owned();
    let cases = [
        ("A", 2, 1u64, 2u64, 1, 3, "1/6"),
        ("B3", -2, 64, 128, 1, 1, "-1/2"),
        ("Mixed", 3, 16, 32, -1, 2, "-1/6"),
    ];
    for (name, m, a, b, num, den, _) in cases {
        let index = name;
        let (signed, absolute) = if m < 0 { (-1, -m) } else { (1, m) };
        writeln!(source, r"
            let s{index}=rvKSNomizuScale[{num},{den}](rvKSNomizuFiber[{name},{m}](rvKSExteriorBasis[P,{a}]));
            let t{index}=rvKSNomizuFiber[{name},{m}](rvKSExteriorBasis[P,{b}]);
            let eta{index}=rvKSNomizuEta[{name},{m}];
            let primitive{index}=rvKSNomizuScale[{signed},{absolute}](s{index});
            assert rvKSNomizuD(primitive{index})==rvKSNomizuCup(eta{index},t{index});
            let rep{index}=rvKSNomizuEtaMassey(s{index});
            assert rep{index}==rvKSNomizuMassey(eta{index},eta{index},t{index},zero[RVKSNomizu({name},{m},1)],primitive{index});
            rvKSNomizuMasseyCertificate(s{index});
            rvKSNomizuCoeff[{a},1](rep{index});
            ").unwrap();
    }
    source.push_str(r"
        let extreme=rvKSNomizuFiber[A,-170141183460469231731687303715884105728](rvKSExteriorBasis[P,1]);
        rvKSNomizuCoeff[1,1](rvKSNomizuEtaMassey(extreme));
        rvKSNomizuMasseyCertificate(extreme);
    ");
    let actual = stdout(interpret(&source));
    let mut lines = actual.lines();
    for (_, _, mask, _, _, _, coefficient) in cases {
        assert_eq!(
            lines.next().unwrap(),
            "rvKSNomizuMasseyCertificate = triple=<eta,eta,Df>, degree=2, excludes_zero=true"
        );
        assert_eq!(
            lines.next().unwrap(),
            format!("rvKSNomizuCoeff[{mask},1] = {coefficient}")
        );
    }
    assert_eq!(
        lines.next().unwrap(),
        "rvKSNomizuCoeff[1,1] = -1/170141183460469231731687303715884105728"
    );
    assert!(lines.next().unwrap().ends_with("excludes_zero=true"));
    assert!(lines.next().is_none());
}

fn choose(n: usize, k: usize) -> u128 {
    (0..k).fold(1, |value, i| value * (n - i) as u128 / (i + 1) as u128)
}

fn report(line: &str) -> Vec<u128> {
    line.split_once(" = ")
        .unwrap()
        .1
        .split(", ")
        .map(|field| field.split_once('=').unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn all_cochain_dimensions_match_wang_and_identity_cover() {
    let mut source = setup().to_owned();
    let modes = [("A", 2), ("B3", -1), ("Mixed", 3), ("Mixed", 0)];
    for (name, m) in modes {
        for degree in 0..=65 {
            writeln!(
                source,
                "rvKSNomizuData[{name},{m},{degree}]; rvKSAlexanderJets[{name},{m},{degree},1,Q];"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    for (_, m) in modes {
        let mut prior_closed = 0;
        for degree in 0..=65 {
            let ch = report(lines.next().unwrap());
            let wang = report(lines.next().unwrap());
            assert_eq!(ch[0], choose(65, degree));
            assert_eq!(ch[3], wang[0]);
            assert_eq!(ch[1] - ch[2], ch[3]);
            let prior = degree.checked_sub(1).map_or(0, |k| choose(65, k));
            assert_eq!(ch[2], prior - prior_closed);
            prior_closed = ch[1];
            if m == 0 {
                assert_eq!(ch, [ch[0], ch[0], 0, ch[0]]);
            }
            if matches!(degree, 0 | 65) {
                assert_eq!(ch[3], 1);
            }
        }
    }
    assert!(lines.next().is_none());
}

#[test]
fn type_nullhomotopy_and_nonformality_claims_have_diagnostics() {
    let prefix = format!(
        r"{}
        let eta=rvKSNomizuEta[A,2];
        let x=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,1]);
        let targetx=rvKSNomizuFiber[A,2](rvKSExteriorBasis[P,2]);
        let hx=rvKSNomizuScale[1,6](x);
        let zz=zero[RVKSNomizu(A,2,1)];
        ",
        setup()
    );
    for (command, expected) in [
        (
            "rvKSNomizuMasseyCertificate(rvKSNomizuFiber[A,0](rvKSExteriorBasis[P,1]));",
            "requires nonzero monodromy exponent",
        ),
        (
            "show rvKSNomizuEtaMassey(add(x,eta));",
            "requires a pure fiber one-form",
        ),
        (
            "rvKSNomizuMasseyCertificate(targetx);",
            "requires D(f) to be nonzero",
        ),
        (
            "rvKSNomizuMasseyCertificate(rvKSNomizuCup(eta,targetx));",
            "expected type RVKSNomizu(A,2,1)",
        ),
        (
            "show rvKSNomizuMassey(eta,eta,x,zz,hx);",
            "requires three closed inputs",
        ),
        (
            "show rvKSNomizuMassey(eta,eta,targetx,hx,hx);",
            "first nullhomotopy has the wrong differential",
        ),
        (
            "show rvKSNomizuMassey(eta,eta,targetx,zz,zz);",
            "second nullhomotopy has the wrong differential",
        ),
        (
            "show rvKSNomizuCup(x,rvKSNomizuEta[A,0]);",
            "expected type RVKSNomizu(A,2,1)",
        ),
        (
            "show add(x,rvKSNomizuEta[B3,2]);",
            "expected type RVKSNomizu(A,2,1)",
        ),
        (
            "show rvKSNomizuScale[1,0](x);",
            "denominator must be positive",
        ),
        (
            "rvKSNomizuCohomEqual(x,targetx,zero[RVKSNomizu(A,2,0)]);",
            "requires closed cochains",
        ),
        ("rvKSNomizuData[A,2,66];", "degree must be between 0 and 65"),
        ("rvKSNomizuCoeff[1,2](x);", "eta exponent must be 0 or 1"),
    ] {
        let actual = diagnostic(interpret(&format!("{prefix}{command}")));
        assert!(actual.contains(expected), "{actual}");
    }
}
