mod common;

use common::{diagnostic, example, interpret, stdout};

const BASIS: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zeroW=rvUnramElem[W]{};
    let oneW=rvUnramElem[W]{1};
    let a=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[L](zeroW,zeroW);
    let f1=rvJacPoint[L](oneW,zeroW);
    let e2=add(add(e1,f1),rvJacPoint[L](neg(oneW),zeroW));
    let f2=add(add(e1,f1),rvJacPoint[L](a,zeroW));
    rvjactorsionbasis B[L,2](neg(oneW)){e1,e2,f1,f2};
";

#[test]
fn parabolic_reciprocal_quadratic_and_salem_words_have_exact_spectra() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-lattice-dynamics")),
        concat!(
            "rvKummerDivPair = 4\n",
            "rvK3Entropy[Parabolic] = characteristic=[1, -17, 136, -680, 2380, -6188, 12376, -19448, 24310, -24310, 19448, -12376, 6188, -2380, 680, -136, 17, -1], cyclotomic=[(1, 17)], dynamical degree=1\n",
            "rvK3Entropy[Quadratic] = characteristic=[1, -31, 312, -1672, 5628, -12740, 19656, -19448, 8294, 8294, -19448, 19656, -12740, 5628, -1672, 312, -31, 1], cyclotomic=[(1, 14), (2, 1)], reciprocal quadratic=[1, -18, 1], dynamical degree=largest real root greater than 1 of [1, -18, 1]\n",
            "rvK3Entropy[Salem] = characteristic=[1, -81, 840, -4072, 11532, -20332, 21528, -10296, -3146, 3146, 10296, -21528, 20332, -11532, 4072, -840, 81, -1], cyclotomic=[(1, 13)], Salem factor=[1, -68, -122, -68, 1], dynamical degree=largest real root greater than 1 of [1, -68, -122, -68, 1]\n",
        ),
    );
}

#[test]
fn generator_formulas_squares_order_and_inverse_hold_on_the_marked_lattice() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{BASIS}
            let h=rvKummerHyperplane[B];
            let e0=rvKummerNodeDiv[B,0];
            let node1=rvKummerNodeDiv[B,1];
            let e3=rvKummerNodeDiv[B,3];
            let d=add(h,add(node1,rvKummerTropeDiv[B,5]));
            rvk3autom Id[B]{{}};
            assert rvK3AutoDiv[Id](d)==d;
            rvk3autom T[B]{{translate(10)}};
            assert rvK3AutoDiv[T](h)==h;
            assert rvK3AutoDiv[T](e3)==rvKummerNodeDiv[B,9];
            assert rvK3AutoDiv[T](rvK3AutoDiv[T](d))==d;
            rvk3autom S[B]{{switch(10)}};
            assert rvK3AutoDiv[S](e3)==rvKummerTropeDiv[B,9];
            assert rvK3AutoDiv[S](h)==rvKummerDiv[B]{{6,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2}};
            assert rvK3AutoDiv[S](rvK3AutoDiv[S](d))==d;
            rvk3autom P0[B]{{project(0)}};
            rvk3autom P1[B]{{project(1)}};
            assert rvK3AutoDiv[P0](e0)==add(rvKummerDivScale[2](h),rvKummerDivScale[-3](e0));
            assert rvK3AutoDiv[P0](node1)==node1;
            assert rvK3AutoDiv[P0](rvK3AutoDiv[P0](d))==d;
            rvk3autom P01[B]{{project(0),project(1)}};
            rvk3autom P10[B]{{project(1),project(0)}};
            assert rvK3AutoDiv[P01](h)==rvK3AutoDiv[P1](rvK3AutoDiv[P0](h));
            assert rvK3AutoDiv[P01](h)==rvKummerDiv[B]{{18,-8,-24,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
            assert rvK3AutoDiv[P10](h)==rvKummerDiv[B]{{18,-24,-8,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
            rvk3autom Forward[B]{{switch(5),translate(8),project(3),project(1)}};
            rvk3autom Reverse[B]{{project(1),project(3),translate(8),switch(5)}};
            assert rvK3AutoDiv[Reverse](rvK3AutoDiv[Forward](d))==d;
            rvKummerDivPair(h,d);
            rvKummerDivPair(rvK3AutoDiv[Forward](h),rvK3AutoDiv[Forward](d));
        "
        ))),
        "rvKummerDivPair = 6\nrvKummerDivPair = 6\n",
    );
}

#[test]
fn translation_conjugacy_reversal_and_mixed_words_match_exact_factors() {
    let output = stdout(interpret(&format!(
        r"{BASIS}
        rvk3autom A[B]{{project(0),project(1),project(2)}};
        rvk3autom I[B]{{project(2),project(1),project(0)}};
        rvk3autom T[B]{{translate(8),project(0),project(1),project(2),translate(8)}};
        rvk3autom U[B]{{project(8),project(9),project(10)}};
        rvk3autom S[B]{{switch(0),project(0),project(1)}};
        let d=rvKummerTropeDiv[B,11];
        assert rvK3AutoDiv[T](d)==rvK3AutoDiv[U](d);
        rvK3Entropy[A];
        rvK3Entropy[I];
        rvK3Entropy[T];
        rvK3Entropy[U];
        rvK3Entropy[S];
    "
    )));
    let lines: Vec<_> = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(lines.len(), 5);
    for line in &lines[1..4] {
        assert_eq!(*line, lines[0]);
    }
    assert!(lines[4].contains("cyclotomic=[(1, 4), (2, 9)], Salem factor=[1, -16, -2, -16, 1]"));
}

#[test]
fn exact_wide_intermediates_cover_long_words_and_cancellation() {
    let triples = vec!["project(0),project(1),project(2)"; 16].join(",");
    let value = i128::MAX - 1;
    let output = stdout(interpret(&format!(
        r"{BASIS}
        rvk3autom Long[B]{{{triples}}};
        rvk3autom P[B]{{project(0)}};
        let x=rvKummerDiv[B]{{{value},-{value},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
        assert rvK3AutoDiv[P](x)==x;
        rvK3Entropy[Long];
        "
    )));
    let (mut previous, mut current) = (2i128, 18i128);
    for _ in 1..16 {
        (previous, current) = (current, 18 * current - previous);
    }
    assert!(output.contains(&format!(
        "cyclotomic=[(1, 15)], reciprocal quadratic=[1, {}, 1]",
        -current
    )));
}

#[test]
fn malformed_words_wrong_markings_and_overflow_are_diagnosed() {
    let triples = vec!["project(0),project(1),project(2)"; 36].join(",");
    for (source, expected) in [
        ("rvk3autom Bad[B]{flip(1)};".into(), "expected translate, switch, or project"),
        ("rvk3autom Bad[B]{project(16)};".into(), "mask must be between 0 and 15"),
        ("rvk3autom Bad[Unknown]{};".into(), "unknown Jacobian torsion basis"),
        ("rvk3autom A[B]{}; rvk3autom A[B]{};".into(), "already defined"),
        ("rvK3Entropy[Unknown];".into(), "unknown K3 automorphism word"),
        ("rvjactorsionbasis Other[L,2](neg(oneW)){e1,e2,f1,f2}; rvk3autom A[B]{}; show rvK3AutoDiv[A](zero[RVKummerDiv(Other)]);".into(), "expected type RVKummerDiv(B)"),
        (format!("rvk3autom Huge[B]{{{triples}}}; show rvK3AutoDiv[Huge](rvKummerHyperplane[B]);"), "K3 automorphism coefficient exceeds signed 128-bit range"),
        (format!("rvk3autom Huge[B]{{{triples}}}; rvK3Entropy[Huge];"), "K3 automorphism coefficient exceeds signed 128-bit range"),
    ] {
        let problem = diagnostic(interpret(&format!("{BASIS}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
