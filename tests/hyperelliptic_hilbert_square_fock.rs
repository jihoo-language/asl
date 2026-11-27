mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[W]{};
    let o=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
    let e1=rvJacPoint[C](z,z);
    let f1=rvJacPoint[C](o,z);
    let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
    let f2=add(add(e1,f1),rvJacPoint[C](i,z));
    rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};
";

#[test]
fn example_covers_the_incidence_sign_point_class_and_unit_denominator() {
    assert_eq!(
        stdout(example("hyperelliptic-hilbert-square-fock")),
        concat!(
            "rvHilbSquareIntegral = 1\n",
            "rvHilbSquareIntegral = -1\n",
            "rvNakajimaPair = -2\n",
            "rvHilbSquareIntegral = -2\n",
            "rvHilbSquareNakajimaDen = 2\n",
        )
    );
}

fn mode_one_pair(i: usize, j: usize) -> String {
    if i == j {
        format!("rvNakajimaMonomial[B]{{(1,{i},2)}}")
    } else {
        format!("rvNakajimaMonomial[B]{{(1,{i},1),(1,{j},1)}}")
    }
}

#[test]
fn every_weight_two_fock_basis_vector_round_trips() {
    let mut monomials = Vec::new();
    for i in 0..24 {
        for j in i..24 {
            monomials.push(mode_one_pair(i, j));
        }
        monomials.push(format!("rvNakajimaMonomial[B]{{(2,{i},1)}}"));
    }
    assert_eq!(monomials.len(), 324);
    for chunk in monomials.chunks(64) {
        let mut source = MARKED.to_owned();
        for (index, monomial) in chunk.iter().enumerate() {
            writeln!(
                source,
                "let x{index}={monomial}; assert rvHilbSquareNakajima[1](rvNakajimaHilbSquare(x{index}))==x{index};"
            )
            .unwrap();
        }
        assert_eq!(stdout(interpret(&source)), "");
    }
}

#[test]
fn all_cup_coordinate_axes_round_trip_after_scaling() {
    for (grade, dimension) in [(0, 1), (2, 23), (4, 276), (6, 23), (8, 1)] {
        for start in (0..dimension).step_by(64) {
            let mut source = MARKED.to_owned();
            for index in start..dimension.min(start + 64) {
                writeln!(
                    source,
                    "let c{index}=rvHilbSquareBasis[B,{grade},{index}]; assert rvNakajimaHilbSquare(rvHilbSquareNakajima[2](c{index}))==rvHilbSquareScale[2,1](c{index});"
                )
                .unwrap();
            }
            assert_eq!(stdout(interpret(&source)), "");
        }
    }
}

#[test]
fn wick_pairing_equals_the_cup_trace_in_every_partition_and_degree() {
    let setup = r"
        let u=rvNakajimaMonomial[B]{(1,0,2)};
        let top=rvNakajimaMonomial[B]{(1,23,2)};
        let eu=rvNakajimaMonomial[B]{(1,1,1),(1,0,1)};
        let ft=rvNakajimaMonomial[B]{(1,2,1),(1,23,1)};
        let p=rvNakajimaMonomial[B]{(1,0,1),(1,23,1)};
        let ee=rvNakajimaMonomial[B]{(1,1,2)};
        let ef=rvNakajimaMonomial[B]{(1,1,1),(1,2,1)};
        let ff=rvNakajimaMonomial[B]{(1,2,2)};
        let d=rvNakajimaMonomial[B]{(2,0,1)};
        let dt=rvNakajimaMonomial[B]{(2,23,1)};
        let de=rvNakajimaMonomial[B]{(2,1,1)};
        let df=rvNakajimaMonomial[B]{(2,2,1)};
        let mixed1=add(add(ee,rvNakajimaScale[-5](de)),rvNakajimaScale[3](p));
        let mixed2=add(rvNakajimaScale[-7](df),rvNakajimaScale[11](ff));
    ";
    let mut source = format!("{MARKED}{setup}");
    let mut expected = String::new();
    for (left, right, result) in [
        ("u", "top", 2),
        ("eu", "ft", -1),
        ("d", "dt", -2),
        ("eu", "dt", 0),
        ("d", "ft", 0),
        ("p", "p", 1),
        ("p", "ef", 0),
        ("ee", "ff", 2),
        ("de", "df", 2),
        ("ee", "df", 0),
        ("de", "p", 0),
        ("mixed1", "mixed2", 92),
    ] {
        writeln!(
            source,
            "rvNakajimaPair({left},{right}); rvHilbSquareIntegral(rvHilbSquareCup(rvNakajimaHilbSquare({left}),rvNakajimaHilbSquare({right})));"
        )
        .unwrap();
        writeln!(expected, "rvNakajimaPair = {result}").unwrap();
        writeln!(expected, "rvHilbSquareIntegral = {result}").unwrap();
    }
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn denominator_is_minimal_and_scale_preserves_marking_and_formal_grades() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let e=rvHilbSquareBasis[B,2,0];
            let d=rvHilbSquareBasis[B,2,22];
            let x=add(rvHilbSquareScale[1,3](e),rvHilbSquareScale[1,5](d));
            let y=rvHilbSquareNakajima[-30](x);
            let fe=rvNakajimaMonomial[B]{{(1,0,1),(1,1,1)}};
            let fd=rvNakajimaMonomial[B]{{(2,0,1)}};
            assert y==add(rvNakajimaScale[-10](fe),rvNakajimaScale[-3](fd));
            assert rvNakajimaHilbSquare(y)==rvHilbSquareScale[-30,1](x);
            assert rvHilbSquareNakajima[0](x)==zero[RVNakajima(B,2,2)];
            assert rvHilbSquareNakajima[3](zero[RVHilbSquare(B,10)])==zero[RVNakajima(B,2,10)];
            assert rvNakajimaHilbSquare(zero[RVNakajima(B,2,3)])==zero[RVHilbSquare(B,3)];
            let negative=rvHilbSquareLower(d,rvHilbSquareBasis[B,0,0]);
            assert rvHilbSquareNakajima[13](negative)==zero[RVNakajima(B,2,-2)];
            assert rvNakajimaHilbSquare(zero[RVNakajima(B,2,-2)])==negative;
            rvHilbSquareNakajimaDen(negative);
            rvHilbSquareNakajimaDen(x);
            rvHilbSquareNakajimaDen(rvHilbSquareCup(d,d));
            rvHilbSquareNakajimaDen(rvHilbSquareC2[B]);
            rvHilbSquareNakajimaDen(zero[RVHilbSquare(B,10)]);
            "
        ))),
        concat!(
            "rvHilbSquareNakajimaDen = 1\n",
            "rvHilbSquareNakajimaDen = 30\n",
            "rvHilbSquareNakajimaDen = 1\n",
            "rvHilbSquareNakajimaDen = 1\n",
            "rvHilbSquareNakajimaDen = 1\n",
        )
    );
}

#[test]
fn errors_check_weight_nominal_markings_integrality_and_range() {
    for (source, expected) in [
        (
            "show rvNakajimaHilbSquare(rvNakajimaVac[B]);",
            "requires a Nakajima state of length two",
        ),
        (
            "show rvNakajimaHilbSquare(z);",
            "a Nakajima state",
        ),
        (
            "rvHilbSquareNakajimaDen(rvNakajimaVac[B]);",
            "an RVHilbSquare cohomology class",
        ),
        (
            "show rvHilbSquareNakajima[1](rvHilbSquareBasis[B,0,0]);",
            "nonintegral Fock coefficients",
        ),
        (
            "let x=add(rvHilbSquareScale[1,3](rvHilbSquareBasis[B,2,0]),rvHilbSquareScale[1,5](rvHilbSquareBasis[B,2,22])); show rvHilbSquareNakajima[15](x);",
            "nonintegral Fock coefficients",
        ),
        (
            "rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; assert rvNakajimaHilbSquare(zero[RVNakajima(D,2,0)])==zero[RVHilbSquare(B,0)];",
            "RVHilbSquare(D, 0)",
        ),
        (
            "rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; assert rvHilbSquareNakajima[1](zero[RVHilbSquare(D,0)])==zero[RVNakajima(B,2,0)];",
            "RVNakajima(D,2,0)",
        ),
        (
            "let u=rvHilbSquareBasis[B,0,0]; show rvHilbSquareNakajima[170141183460469231731687303715884105726](rvHilbSquareScale[4,1](u));",
            "exceeds the signed 128-bit coefficient range",
        ),
        (
            "let a=rvHilbSquareScale[1,170141183460469231731687303715884105727](rvHilbSquareBasis[B,2,0]); let b=rvHilbSquareScale[1,2](rvHilbSquareBasis[B,2,1]); rvHilbSquareNakajimaDen(add(a,b));",
            "denominator exceeds the signed 128-bit range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{MARKED}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
