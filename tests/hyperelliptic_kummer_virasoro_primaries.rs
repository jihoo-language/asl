mod common;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zeroW=rvUnramElem[W]{};
    let oneW=rvUnramElem[W]{1};
    let a=rvUnramElem[W]{0,1};
    let r0=rvJacPoint[L](zeroW,zeroW);
    let r1=rvJacPoint[L](oneW,zeroW);
    let r2=rvJacPoint[L](neg(oneW),zeroW);
    let ra=rvJacPoint[L](a,zeroW);
    rvjactorsionbasis B[L,2](neg(oneW)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),ra)};
    let v=rvNakajimaVac[B];
";

#[test]
fn highest_weight_example_fixes_free_column_order_and_primitive_sign() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-virasoro-primaries")),
        "rvNakajimaPrimaryDim[B, 2, 4] = 253\nrvNakajimaIsPrimary = true\nrvNakajimaIsPrimary = false\n"
    );
}

#[test]
fn grades_recover_the_lowest_weight_dimensions() {
    let mut source = MARKED.to_owned();
    let mut expected = String::new();
    let gradings: &[&[usize]] = &[
        &[1],
        &[1, 22, 1],
        &[1, 22, 253, 22, 1],
        &[1, 22, 275, 2256, 275, 22, 1],
    ];
    for (length, dimensions) in gradings.iter().enumerate() {
        for (half_degree, dimension) in dimensions.iter().enumerate() {
            let degree = half_degree * 2;
            source.push_str(&format!("rvNakajimaPrimaryDim[B,{length},{degree}];"));
            expected.push_str(&format!(
                "rvNakajimaPrimaryDim[B, {length}, {degree}] = {dimension}\n"
            ));
        }
    }
    assert_eq!(
        gradings
            .iter()
            .map(|grades| grades.iter().sum::<usize>())
            .collect::<Vec<_>>(),
        [1, 24, 299, 2852]
    );
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn each_positive_constraint_is_necessary() {
    assert_eq!(
        stdout(interpret(&format!(r"{MARKED}
            rvNakajimaIsPrimary(v);
            let x = rvNakajimaMonomial[B]{{(2,0,1)}};
            let y = rvNakajimaMonomial[B]{{(1,0,1),(1,23,1)}};
            assert rvVirasoro[2](x) == zero[RVNakajima(B,0,-2)];
            assert rvVirasoro[1](x) == rvNakajimaScale[2](rvNakajimaCreate[1,0](v));
            assert rvVirasoro[1](y) == zero[RVNakajima(B,1,2)];
            assert rvVirasoro[2](y) == v;
            rvNakajimaIsPrimary(x);
            rvNakajimaIsPrimary(y);
            let twiceY = rvNakajimaScale[2](y);
            let b = rvNakajimaMonomial[B]{{(1,1,2)}};
            rvNakajimaIsPrimary(add(twiceY,b));
            rvNakajimaIsPrimary(zero[RVNakajima(B,-5,7)]);
        "))),
        "rvNakajimaIsPrimary = true\nrvNakajimaIsPrimary = false\nrvNakajimaIsPrimary = false\nrvNakajimaIsPrimary = true\nrvNakajimaIsPrimary = true\n"
    );
}

#[test]
fn generator_values_cover_nontrivial_grades_and_nominal_markings() {
    let mut source = format!(
        r"{MARKED}
        assert rvNakajimaPrimaryGen[B,1,0,0] == rvNakajimaCreate[1,0](v);
        assert rvNakajimaPrimaryGen[B,1,2,0] == rvNakajimaCreate[1,1](v);
        assert rvNakajimaPrimaryGen[B,1,2,21] == rvNakajimaCreate[1,22](v);
        assert rvNakajimaPrimaryGen[B,1,4,0] == rvNakajimaCreate[1,23](v);
        assert rvNakajimaPrimaryGen[B,3,0,0] == rvNakajimaMonomial[B]{{(1,0,3)}};
        assert rvNakajimaPrimaryGen[B,3,12,0] == rvNakajimaMonomial[B]{{(1,23,3)}};
    "
    );
    for (degree, index) in [
        (2, 21),
        (4, 171),
        (6, 0),
        (6, 1128),
        (6, 2255),
        (8, 255),
        (10, 21),
    ] {
        source.push_str(&format!(
            r"
            let p{degree}x{index} = rvNakajimaPrimaryGen[B,3,{degree},{index}];
            assert rvVirasoro[1](p{degree}x{index}) == zero[RVNakajima(B,2,{})];
            assert rvVirasoro[2](p{degree}x{index}) == zero[RVNakajima(B,1,{})];
            assert rvVirasoro[3](p{degree}x{index}) == zero[RVNakajima(B,0,{})];
        ",
            degree - 2,
            degree - 4,
            degree - 6
        ));
    }
    assert_eq!(stdout(interpret(&source)), "");

    let wrong_marking = format!(
        r"{MARKED}
        rvjactorsionbasis H[L,2](neg(oneW)){{r0,add(add(r0,r1),r2),r1,add(add(r0,r1),ra)}};
        assert rvNakajimaPrimaryGen[H,0,0,0] == v;
    "
    );
    let err = diagnostic(interpret(&wrong_marking));
    assert!(
        err.contains("RVNakajima(H,0,0)") && err.contains("RVNakajima(B,0,0)"),
        "{err}"
    );
}

#[test]
fn generator_clears_nontrivial_rational_denominators() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let p = rvNakajimaPrimaryGen[B,4,4,275];
            let free = rvNakajimaScale[-3](rvNakajimaMonomial[B]{{(2,0,2)}});
            let first = rvNakajimaScale[4](rvNakajimaMonomial[B]{{(1,0,1),(3,0,1)}});
            let second = rvNakajimaScale[12](rvNakajimaMonomial[B]{{(1,0,2),(1,1,1),(1,2,1)}});
            assert p == add(add(first,second),free);
            assert rvVirasoro[1](p) == zero[RVNakajima(B,3,2)];
            assert rvVirasoro[2](p) == zero[RVNakajima(B,2,0)];
        "
        ))),
        ""
    );
}

#[test]
fn sparse_checks_support_huge_grades_and_wide_images() {
    let length = 1_000_000_000_000_000_000_000_000i128;
    let exponent = i128::MAX - 2;
    let paired = exponent + 1;
    assert_eq!(
        stdout(interpret(&format!(r"{MARKED}
            rvNakajimaPrimaryDim[B,{length},0];
            rvNakajimaPrimaryDim[B,{length},2];
            assert rvNakajimaPrimaryGen[B,{length},0,0] == rvNakajimaMonomial[B]{{(1,0,{length})}};
            rvNakajimaIsPrimary(rvNakajimaMonomial[B]{{(1,0,{length})}});
            rvNakajimaIsPrimary(rvNakajimaScale[{}](rvNakajimaMonomial[B]{{({length},0,1)}}));
            rvNakajimaIsPrimary(rvNakajimaMonomial[B]{{(1,0,{length}),(1,23,1)}});
            let left = rvNakajimaScale[2](rvNakajimaMonomial[B]{{(1,0,{paired}),(1,23,1)}});
            let right = rvNakajimaScale[{paired}](rvNakajimaMonomial[B]{{(1,0,{exponent}),(1,1,2)}});
            rvNakajimaIsPrimary(add(left,right));
        ", i128::MAX))),
        format!("rvNakajimaPrimaryDim[B, {length}, 0] = 1\nrvNakajimaPrimaryDim[B, {length}, 2] = 22\nrvNakajimaIsPrimary = true\nrvNakajimaIsPrimary = false\nrvNakajimaIsPrimary = false\nrvNakajimaIsPrimary = true\n")
    );
}

#[test]
fn unsupported_grades_and_limited_enumeration_are_diagnosed() {
    assert_eq!(
        stdout(interpret(&format!(r"{MARKED}
            rvNakajimaPrimaryDim[B,-1,0];
            rvNakajimaPrimaryDim[B,1,3];
            rvNakajimaPrimaryDim[B,1,6];
        "))),
        "rvNakajimaPrimaryDim[B, -1, 0] = 0\nrvNakajimaPrimaryDim[B, 1, 3] = 0\nrvNakajimaPrimaryDim[B, 1, 6] = 0\n"
    );
    for (input, expected) in [
        (
            "show rvNakajimaPrimaryGen[B,2,2,22];",
            "primary generator index is out of range",
        ),
        (
            "show rvNakajimaPrimaryGen[B,1,3,0];",
            "primary generator index is out of range",
        ),
        (
            "rvNakajimaPrimaryDim[B,1000000000,2000000000];",
            "graded enumeration storage limit",
        ),
        ("rvNakajimaIsPrimary(o);", "a Nakajima state"),
        (
            "rvNakajimaPrimaryDim[C,0,0];",
            "unknown Jacobian torsion basis",
        ),
    ] {
        let err = diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(err.contains(expected), "{err}");
    }
}
