mod common;

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
fn example_exhibits_the_rational_unit_and_character_above_the_rank() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-rational-hilbert")),
        "rvHilbCohDen = 2\nrvHilbCohIntegral = 1\nrvHilbCohIntegral = 3/4\n"
    );
}

#[test]
fn weight_zero_and_geometric_units_through_four_points() {
    let mut source = format!(
        r"{MARKED}
        let v=rvHilbCohUnit[B,0];
        let point=rvHilbCohPoint[B,0];
        assert v==point;
        assert rvHilbCohCup(v,point)==v;
        assert rvHilbCohFock[1](v)==rvNakajimaVac[B];
        rvHilbCohIntegral(v);
    "
    );
    for (n, factorial, factors) in [
        (1, 1, "(1,1,1)"),
        (2, 2, "(2,1,1)"),
        (3, 6, "(1,0,1),(2,1,1)"),
        (4, 24, "(1,0,1),(3,1,1)"),
    ] {
        source.push_str(&format!(
            r"
            let u{n}=rvHilbCohUnit[B,{n}];
            let p{n}=rvHilbCohPoint[B,{n}];
            let x{n}=rvHilbCohQ[5,7](rvNakajimaMonomial[B]{{{factors}}});
            assert rvHilbCohCup(u{n},x{n})==x{n};
            assert rvHilbCohCup(x{n},u{n})==x{n};
            assert rvHilbCohCup(u{n},u{n})==u{n};
            assert rvHilbCohFock[{factorial}](u{n})==rvNakajimaMonomial[B]{{(1,0,{n})}};
            rvHilbCohIntegral(rvHilbCohCup(u{n},p{n}));
        "
        ));
    }
    assert_eq!(
        stdout(interpret(&source)),
        "rvHilbCohIntegral = 1\n".repeat(5)
    );
}

#[test]
fn nonintegral_hilbert_square_classes_transport_along_the_fock_bridge() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let x=rvHilbSquareScale[1,3](rvHilbSquareBasis[B,2,0]);
        let y=rvHilbSquareScale[-1,5](rvHilbSquareBasis[B,2,22]);
        let qx=rvHilbCohQ[1,3](rvHilbSquareNakajima[3](x));
        let qy=rvHilbCohQ[1,10](rvHilbSquareNakajima[10](y));
        let qproduct=rvHilbCohCup(qx,qy);
        rvHilbCohDen(qproduct);
        let bridge=rvHilbSquareScale[1,30](rvNakajimaHilbSquare(rvHilbCohFock[30](qproduct)));
        assert bridge==rvHilbSquareCup(x,y);
        let line=rvK3LatticeBasis[B,0];
        let c1=rvHilbCohChern[2,1](line);
        let c1Bridge=rvHilbSquareScale[1,2](rvNakajimaHilbSquare(rvHilbCohFock[2](c1)));
        assert c1Bridge==add(rvHilbSquareBasis[B,2,0],neg(rvHilbSquareBasis[B,2,22]));
    "
        ))),
        "rvHilbCohDen = 15\n"
    );
}

#[test]
fn line_character_and_unscaled_chern_have_the_expected_grades() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let line=rvK3LatticeBasis[B,0];
        assert rvHilbCohCh[0,0](line)==zero[RVHilbCoh(B,0,0)];
        assert rvHilbCohChern[0,0](line)==rvHilbCohUnit[B,0];
        assert rvHilbCohCh[1,0](line)==rvHilbCohUnit[B,1];
        assert rvHilbCohCh[1,1](line)==rvHilbCohChern[1,1](line);
        assert rvHilbCohCh[1,2](line)==neg(rvHilbCohPoint[B,1]);
        assert rvHilbCohChern[1,2](line)==zero[RVHilbCoh(B,1,4)];
        assert rvHilbCohCh[1,3](line)==zero[RVHilbCoh(B,1,6)];
        assert rvHilbCohCh[3,0](line)==rvHilbCohScale[3,1](rvHilbCohUnit[B,3]);
        assert rvHilbCohChern[3,4](line)==zero[RVHilbCoh(B,3,8)];
        assert rvHilbCohCh[3,7](line)==zero[RVHilbCoh(B,3,14)];
        rvHilbCohIntegral(rvHilbCohCh[1,2](line));
    "
        ))),
        "rvHilbCohIntegral = -1\n"
    );
}

#[test]
fn the_degree_zero_rank_uses_the_reduced_factorial() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let line=rvK3LatticeBasis[B,0];
        rvHilbCohDen(rvHilbCohUnit[B,34]);
        let rank=rvHilbCohCh[35,0](line);
        rvHilbCohDen(rank);
        let expected=rvHilbCohQ[1,295232799039604140847618609643520000000](rvNakajimaMonomial[B]{{(1,0,35)}});
        assert rank==expected;
    "
        ))),
        "rvHilbCohDen = 295232799039604140847618609643520000000\n".repeat(2)
    );
}

#[test]
fn rank_three_character_satisfies_newton_in_the_cup_algebra() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let line=rvK3LatticeBasis[B,0];
        let c1=rvHilbCohChern[3,1](line);
        let c2=rvHilbCohChern[3,2](line);
        let c3=rvHilbCohChern[3,3](line);
        let p2=rvHilbCohScale[2,1](rvHilbCohCh[3,2](line));
        let p3=rvHilbCohScale[6,1](rvHilbCohCh[3,3](line));
        let p4=rvHilbCohScale[24,1](rvHilbCohCh[3,4](line));
        assert p2==add(rvHilbCohCup(c1,c1),rvHilbCohScale[-2,1](c2));
        assert p3==add(add(rvHilbCohCup(c1,p2),neg(rvHilbCohCup(c2,c1))),rvHilbCohScale[3,1](c3));
        assert p4==add(add(rvHilbCohCup(c1,p3),neg(rvHilbCohCup(c2,p2))),rvHilbCohCup(c3,c1));
        let c1sq=rvHilbCohCup(c1,c1);
        let ch4=rvHilbCohCh[3,4](line);
        rvHilbCohIntegral(rvHilbCohCup(ch4,c1sq));
        rvHilbCohIntegral(rvHilbCohCh[3,6](line));
    "
        )))
        .lines()
        .filter(|line| *line != "rvHilbCohIntegral = 0")
        .count(),
        2
    );
}

#[test]
fn fractions_normalize_before_narrowing_and_odd_grades_stay_zero() {
    let m = i128::MAX;
    let before = m - 1;
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let raw=rvNakajimaMonomial[B]{{(1,23,2)}};
        let u=rvHilbCohUnit[B,2];
        let p=rvHilbCohPoint[B,2];
        let huge=rvHilbCohQ[{m},1](raw);
        assert rvHilbCohCup(u,huge)==huge;
        let a=rvHilbCohScale[{before},{m}](u);
        let b=rvHilbCohScale[{m},{before}](p);
        assert rvHilbCohCup(a,b)==p;
        let third=rvHilbCohScale[1,3](p);
        let quarter=rvHilbCohScale[1,4](p);
        assert add(third,quarter)==rvHilbCohQ[7,12](raw);
        assert add(third,neg(third))==zero[RVHilbCoh(B,2,8)];
        rvHilbCohDen(rvHilbCohScale[6,8](p));
        rvHilbCohIntegral(add(third,quarter));
        let formal=zero[RVHilbCoh(B,2,1)];
        assert rvHilbCohCup(formal,p)==zero[RVHilbCoh(B,2,9)];
        let enormous=rvHilbCohPoint[B,10000000000000000000000];
        rvHilbCohIntegral(enormous);
    "
        ))),
        "rvHilbCohDen = 4\nrvHilbCohIntegral = 7/12\nrvHilbCohIntegral = 1\n"
    );
}

#[test]
fn nominal_grades_and_exact_range_errors_are_reported() {
    let m = i128::MAX;
    for (source, expected) in [
        ("show rvHilbCohFock[1](rvHilbCohUnit[B,2]);".to_owned(), "nonintegral Fock coefficients"),
        ("show rvHilbCohScale[1,0](rvHilbCohUnit[B,0]);".to_owned(), "denominator must be positive"),
        ("rvHilbCohIntegral(rvHilbCohUnit[B,2]);".to_owned(), "requires the top cohomological degree"),
        ("rvjactorsionbasis H[C,2](neg(o)){e1,e2,f1,f2}; show rvHilbCohCup(rvHilbCohUnit[B,0],rvHilbCohUnit[H,0]);".to_owned(), "same named level basis and length"),
        ("show rvHilbCohCup(rvHilbCohUnit[B,0],rvHilbCohUnit[B,1]);".to_owned(), "same named level basis and length"),
        ("show rvHilbCohUnit[B,35];".to_owned(), "denominator exceeds the unsigned 128-bit range"),
        (format!("let x=rvHilbCohScale[{m},1](rvHilbCohPoint[B,2]); show add(x,x);"), "coefficient exceeds the signed 128-bit range"),
        (format!("show rvHilbCohScale[1,{m}](rvHilbCohScale[1,{m}](rvHilbCohPoint[B,2]));"), "denominator exceeds the unsigned 128-bit range"),
        (format!("show rvHilbCohCup(zero[RVHilbCoh(B,1,{m})],rvHilbCohPoint[B,1]);"), "grade exceeds signed 128-bit range"),
    ] {
        let error=diagnostic(interpret(&format!("{MARKED}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
