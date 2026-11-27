mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const SURFACE: &str = "rvhodgewick B[21]{0,1}{[2]};";

#[test]
fn all_internal_dimensions_through_arity_four_match_truncated_polynomial_oracle() {
    let mut source = String::from(SURFACE);
    let mut grades = Vec::new();
    for arity in 0..=4 {
        for weight in (-2 * arity - 1)..=(3 - arity) {
            writeln!(source, "rvVerbitskyHHDims[B,1,{arity},{weight}];").unwrap();
            grades.push((arity, weight));
        }
    }
    let output = stdout(interpret(&source));
    let dims = output
        .lines()
        .map(|line| {
            line.split_once(", HH=")
                .unwrap()
                .1
                .trim_end_matches(')')
                .parse::<usize>()
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(grades.len(), dims.len());
    let support = grades
        .into_iter()
        .zip(dims)
        .filter_map(|(grade, dimension)| {
            if dimension == 0 {
                None
            } else {
                Some((grade, dimension))
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(
        support,
        [
            ((0, 0), 1),
            ((0, 1), 1),
            ((0, 2), 1),
            ((1, 0), 1),
            ((1, 1), 1),
            ((2, -3), 1),
            ((2, -2), 1),
            ((3, -3), 1),
            ((3, -2), 1),
            ((4, -6), 1),
            ((4, -5), 1),
        ]
    );
}

#[test]
fn carry_is_associative_deformation_and_its_higher_cups_are_nonzero() {
    let output = stdout(example("hyperelliptic-kummer-verbitsky-hochschild"));
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 6);
    assert_eq!(lines[3], "rvVerbitskyHHCoeff = 1");
    assert!(lines[4].starts_with("RVVerbitskyHH(B,1,3,-3)[([1, 2, 1], 1, 1), "));
    assert!(lines[5].starts_with("RVVerbitskyHH(B,1,4,-6)[([1, 2, 1, 2], 0, 1), "));
}

#[test]
fn direct_coboundary_braces_bv_and_gerstenhaber_identities() {
    let output = stdout(interpret(&format!(
        "{SURFACE}
        let x=rvVerbitskyHHBasis[B,1,0,1,0];
        let d=rvVerbitskyHHClass(add(rvVerbitskyHHTerm[B,1,1,0,1]{{1}},rvVerbitskyHHTerm[B,1,1,0,2,2]{{2}}));
        let t=rvVerbitskyHHBasis[B,1,2,-3,0];
        let dt=rvVerbitskyHHCup(d,t);
        let xt=rvVerbitskyHHCup(x,t);
        assert rvVerbitskyHHBV(dt)==rvVerbitskyHHScale[5,1](t);
        assert rvVerbitskyHHBV(rvVerbitskyHHBV(dt))==zero[RVVerbitskyHH(B,1,1,-3)];
        assert rvVerbitskyHHBracket(d,xt)==rvVerbitskyHHScale[-2,1](xt);
        assert rvVerbitskyHHBracket(d,dt)==rvVerbitskyHHScale[-3,1](dt);
        assert rvVerbitskyHHBracket(d,xt)==add(rvVerbitskyHHCup(rvVerbitskyHHBracket(d,x),t),rvVerbitskyHHCup(x,rvVerbitskyHHBracket(d,t)));
        assert rvVerbitskyHHBracket(d,t)==neg(add(rvVerbitskyHHBV(dt),neg(rvVerbitskyHHCup(rvVerbitskyHHBV(d),t))));
        assert rvVerbitskyHHBracket(x,x)==zero[RVVerbitskyHH(B,1,-1,2)];
        assert rvVerbitskyHHBracket(x,d)==neg(x);
        let tx=rvVerbitskyHHRepresentative(t);
        let xc=rvVerbitskyHHRepresentative(x);
        assert rvVerbitskyHHBrace(tx,xc)==zero[RVVerbitskyHochschildCochain(B,1,1,-2)];
        assert rvVerbitskyHHBrace(rvVerbitskyHHRepresentative(d),xc)==xc;
        let f=rvVerbitskyHHTerm[B,1,1,-1,0]{{1}};
        let delta=rvVerbitskyHHDiff(f);
        assert rvVerbitskyHHDiff(delta)==zero[RVVerbitskyHochschildCochain(B,1,3,-1)];
        assert rvVerbitskyHHClass(delta)==zero[RVVerbitskyHH(B,1,2,-1)];
        rvVerbitskyHHCoeff[1]{{1,1}}(delta);
        let u=rvVerbitskyHHTerm[B,1,2,-1,2]{{1,2}};
        let dbv=rvVerbitskyHHDiff(rvVerbitskyHHBV(u));
        assert add(dbv,rvVerbitskyHHBV(rvVerbitskyHHDiff(u)))==zero[RVVerbitskyHochschildCochain(B,1,2,-1)];
        rvVerbitskyHHCoeff[1]{{1,1}}(dbv);
        let j=rvVerbitskyTangent[B,1,-2]{{(0,rvVerbitskyScale[-2,21](rvVerbitskyBasis[B,1,0,0])),(1,rvVerbitskyBasis[B,1,2,0])}};
        assert rvVerbitskyHHFromTangent(j)==zero[RVVerbitskyHH(B,1,2,-2)];
        "
    )));
    assert_eq!(output, "rvVerbitskyHHCoeff = 2\nrvVerbitskyHHCoeff = -3\n");
}

#[test]
fn non_diagonal_frobenius_pairing_detects_bivector_bv_sign_and_jacobi() {
    let output = stdout(interpret(
        r"
        rvhodgewick C[20]{0,1}{[2,1],[1,-2]};
        rvVerbitskyIntegral(rvVerbitskyBasis[C,1,4,0]);
        rvVerbitskyHHDims[C,1,2,0];
        let p=rvVerbitskyHHClass(add(rvVerbitskyHHTerm[C,1,2,0,3]{1,2},rvVerbitskyHHTerm[C,1,2,0,3,-1]{2,1}));
        let x=rvVerbitskyHHBasis[C,1,0,1,0];
        let d=rvVerbitskyHHClass(add(add(rvVerbitskyHHTerm[C,1,1,0,1]{1},rvVerbitskyHHTerm[C,1,1,0,2]{2}),rvVerbitskyHHTerm[C,1,1,0,3,2]{3}));
        let v=rvVerbitskyHHBracket(p,x);
        assert v==rvVerbitskyHHClass(rvVerbitskyHHTerm[C,1,1,1,3,2]{2});
        assert v==add(rvVerbitskyHHCup(rvVerbitskyHHBV(p),x),neg(rvVerbitskyHHBV(rvVerbitskyHHCup(p,x))));
        assert rvVerbitskyHHBV(rvVerbitskyHHBV(p))==zero[RVVerbitskyHH(C,1,0,0)];
        assert rvVerbitskyHHBracket(d,v)==add(rvVerbitskyHHBracket(rvVerbitskyHHBracket(d,p),x),rvVerbitskyHHBracket(p,rvVerbitskyHHBracket(d,x)));
        rvVerbitskyHHCoeff[1]{1}(rvVerbitskyHHBV(p));
        rvVerbitskyHHCoeff[2]{1}(rvVerbitskyHHBV(p));
        rvVerbitskyHHCoeff[1]{2}(rvVerbitskyHHBV(p));
        rvVerbitskyHHCoeff[2]{2}(rvVerbitskyHHBV(p));
        ",
    ));
    assert_eq!(output, "rvVerbitskyIntegral = 1\nrvVerbitskyHHDims[C, 1, 2, 0] = (C=4, Z=4, B=3, HH=1)\nrvVerbitskyHHCoeff = -2/5\nrvVerbitskyHHCoeff = 4/5\nrvVerbitskyHHCoeff = 4/5\nrvVerbitskyHHCoeff = 2/5\n");
}

#[test]
fn reports_coordinate_normalization_type_and_work_errors() {
    for (query, reason) in [
        ("let x=rvVerbitskyHHTerm[B,1,1,0,0]{0};", "positive weight"),
        ("let x=rvVerbitskyHHTerm[B,1,2,0,0]{1};", "wrong arity"),
        ("let x=rvVerbitskyHHTerm[B,1,1,0,2]{1};", "wrong internal shift"),
        ("let x=rvVerbitskyHHTerm[B,1,1,0,0]{3};", "outside the ring basis"),
        ("let x=rvVerbitskyHHClass(rvVerbitskyHHTerm[B,1,1,-1,0]{1});", "is not a cocycle"),
        ("let x=rvVerbitskyHHBasis[B,1,3,-3,1];", "basis index out of range"),
        ("rvVerbitskyHHDims[B,1,1000000000000,-1000000000000];", "work or storage limit"),
        ("let x=rvVerbitskyHHBrace(rvVerbitskyHHBasis[B,1,1,0,0],rvVerbitskyHHBasis[B,1,0,1,0]);", "with cochains for braces"),
        ("let x=rvVerbitskyHHScale[1,0](rvVerbitskyHHBasis[B,1,0,0,0]);", "denominator must be positive"),
        ("rvhodgewick C[21]{0,1}{[2]}; let x=rvVerbitskyHHCup(rvVerbitskyHHBasis[B,1,0,0,0],rvVerbitskyHHBasis[C,1,0,0,0]);", "same marking"),
    ] {
        assert!(diagnostic(interpret(&format!("{SURFACE} {query}"))).contains(reason), "{query}");
    }
}
