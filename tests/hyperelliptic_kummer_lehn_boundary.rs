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
fn half_boundary_example_matches_the_hilbert_cube_fujiki_value() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-lehn-boundary")),
        "rvKummerFujiki[3, 6, 0, 1] = -960\n"
    );
}

#[test]
fn top_intersection_for_four_points_has_positive_sign() {
    let mut source = format!("{MARKED} let d0=rvNakajimaMonomial[B]{{(1,0,4)}};");
    for power in 1..=8 {
        writeln!(source, "let d{power}=rvNakajimaBoundary(d{});", power - 1).unwrap();
    }
    source.push_str(
        r"
        assert d8==rvNakajimaScale[3265920](rvNakajimaMonomial[B]{(1,23,4)});
        let zeroDiv=zero[RVKummerDiv(B)];
        rvKummerFujiki[4,8,0,1](zeroDiv,zeroDiv);
    ",
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvKummerFujiki[4, 8, 0, 1] = 136080\n"
    );
}

#[test]
fn all_weight_two_fock_axes_agree_with_geometric_half_boundary_cup() {
    let mut axes = Vec::new();
    for i in 0..24 {
        for j in i..24 {
            let factors = if i == j {
                format!("(1,{i},2)")
            } else {
                format!("(1,{i},1),(1,{j},1)")
            };
            axes.push(format!("rvNakajimaMonomial[B]{{{factors}}}"));
        }
        axes.push(format!("rvNakajimaMonomial[B]{{(2,{i},1)}}"));
    }
    assert_eq!(axes.len(), 324);
    for chunk in axes.chunks(64) {
        let mut source = format!("{MARKED} let delta=rvHilbSquareBasis[B,2,22];");
        for (index, monomial) in chunk.iter().enumerate() {
            writeln!(source,
                "let x{index}={monomial}; assert rvNakajimaHilbSquare(rvNakajimaBoundary(x{index}))==rvHilbSquareCup(delta,rvNakajimaHilbSquare(x{index}));"
            ).unwrap();
        }
        assert_eq!(stdout(interpret(&source)), "");
    }
}

#[test]
fn joins_and_splits_include_equal_lengths_and_repeated_copies() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let one=rvNakajimaMonomial[B]{{(1,0,3)}};
        assert rvNakajimaBoundary(one)==rvNakajimaScale[3](rvNakajimaMonomial[B]{{(1,0,1),(2,0,1)}});
        let pair=rvNakajimaMonomial[B]{{(1,0,1),(2,23,1)}};
        let joined=rvNakajimaScale[2](rvNakajimaMonomial[B]{{(3,23,1)}});
        let split=neg(rvNakajimaMonomial[B]{{(1,0,1),(1,23,2)}});
        assert rvNakajimaBoundary(pair)==add(joined,split);
        let twin=rvNakajimaMonomial[B]{{(2,23,2)}};
        assert rvNakajimaBoundary(twin)==rvNakajimaScale[-2](rvNakajimaMonomial[B]{{(1,23,2),(2,23,1)}});
        let p=rvNakajimaMonomial[B]{{(4,23,1)}};
        let splitDistinct=rvNakajimaScale[-4](rvNakajimaMonomial[B]{{(1,23,1),(3,23,1)}});
        let splitEqual=rvNakajimaScale[-2](twin);
        assert rvNakajimaBoundary(p)==add(splitDistinct,splitEqual);
        let pMiddle=rvNakajimaMonomial[B]{{(4,1,1)}};
        let first=rvNakajimaMonomial[B]{{(1,1,1),(3,23,1)}};
        let second=rvNakajimaMonomial[B]{{(1,23,1),(3,1,1)}};
        let third=rvNakajimaMonomial[B]{{(2,1,1),(2,23,1)}};
        assert rvNakajimaBoundary(pMiddle)==rvNakajimaScale[-4](add(add(first,second),third));
    "
        ))),
        ""
    );
}

#[test]
fn joining_and_splitting_are_adjoint_for_wick_pairing() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let x=rvNakajimaMonomial[B]{{(1,0,1),(2,0,1)}};
        let y=rvNakajimaMonomial[B]{{(3,23,1)}};
        rvNakajimaPair(rvNakajimaBoundary(x),y);
        rvNakajimaPair(x,rvNakajimaBoundary(y));
        let r=rvNakajimaMonomial[B]{{(1,1,1),(2,0,1)}};
        let s=rvNakajimaMonomial[B]{{(1,0,1),(1,2,1),(1,23,1)}};
        rvNakajimaPair(rvNakajimaBoundary(r),s);
        rvNakajimaPair(r,rvNakajimaBoundary(s));
    "
        ))),
        concat!(
            "rvNakajimaPair = 6\n",
            "rvNakajimaPair = 6\n",
            "rvNakajimaPair = 2\n",
            "rvNakajimaPair = 2\n"
        )
    );
}

#[test]
fn sparse_multiplicity_wide_cancellation_and_formal_grades() {
    let n = 1_000_000_000_000_000_000i128;
    let copies = n - 2;
    let join = n * (n - 1) / 2;
    let max = i128::MAX;
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let x=rvNakajimaMonomial[B]{{(1,0,{n})}};
        let expected=rvNakajimaScale[{join}](rvNakajimaMonomial[B]{{(1,0,{copies}),(2,0,1)}});
        assert rvNakajimaBoundary(x)==expected;
        let a=rvNakajimaMonomial[B]{{(1,1,2)}};
        let b=rvNakajimaMonomial[B]{{(1,2,2)}};
        let wide=add(rvNakajimaScale[{max}](a),rvNakajimaScale[-{max}](b));
        assert rvNakajimaBoundary(wide)==zero[RVNakajima(B,2,6)];
        let v=rvNakajimaVac[B];
        assert rvNakajimaBoundary(v)==zero[RVNakajima(B,0,2)];
        assert rvNakajimaBoundary(rvNakajimaCreate[1,1](v))==zero[RVNakajima(B,1,4)];
        assert rvNakajimaBoundary(zero[RVNakajima(B,-2,5)])==zero[RVNakajima(B,-2,7)];
        assert rvNakajimaBoundary(zero[RVNakajima(B,{n},0)])==zero[RVNakajima(B,{n},2)];
    "
        ))),
        ""
    );
}

#[test]
fn oversized_splits_and_checked_coordinates_are_diagnosed() {
    for (input, expected) in [
        (
            "show rvNakajimaBoundary(rvNakajimaMonomial[B]{(1000000000000000000,23,1)});"
                .to_owned(),
            "statement exceeds the computation limit",
        ),
        (
            format!(
                "show rvNakajimaBoundary(rvNakajimaScale[{}](rvNakajimaMonomial[B]{{(1,0,3)}}));",
                i128::MAX
            ),
            "coordinate exceeds signed 128-bit range",
        ),
        (
            format!(
                "show rvNakajimaBoundary(zero[RVNakajima(B,2,{})]);",
                i128::MAX
            ),
            "grade exceeds signed 128-bit range",
        ),
        ("show rvNakajimaBoundary(o);".to_owned(), "a Nakajima state"),
    ] {
        let err = diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(err.contains(expected), "{err}");
    }
}
