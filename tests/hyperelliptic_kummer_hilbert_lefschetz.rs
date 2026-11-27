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
fn three_point_example_inverts_the_complete_degree_two_slice() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-hilbert-lefschetz")),
        "rvHilbCohDen = 6\n"
    );
}

#[test]
fn three_point_inverse_also_resolves_the_degree_four_slice() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let h=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{{(1,0,2),(1,1,1)}});
        let first=rvHilbCohQ[2,5](rvNakajimaMonomial[B]{{(1,0,1),(1,1,1),(1,2,1)}});
        let second=rvHilbCohQ[-3,7](rvNakajimaMonomial[B]{{(3,0,1)}});
        let p=add(first,second);
        let h2=rvHilbCohCup(h,h);
        let h2p=rvHilbCohCup(h2,p);
        assert rvHilbCohLefschetzInverse(h,h2p)==p;
    "
        ))),
        ""
    );
}

#[test]
fn normalized_lowering_matches_the_square_bridge_in_every_degree() {
    let mut source = format!(
        r"{MARKED}
        let hs=add(rvHilbSquareBasis[B,2,0],rvHilbSquareBasis[B,2,22]);
        let h=rvHilbCohQ[1,2](rvHilbSquareNakajima[2](hs));
    "
    );
    for (degree, index) in [(2, 1), (4, 0), (4, 22), (4, 275), (6, 4), (6, 22), (8, 0)] {
        source.push_str(&format!(r"
            let s{degree}x{index}=rvHilbSquareBasis[B,{degree},{index}];
            let q{degree}x{index}=rvHilbCohQ[1,2](rvHilbSquareNakajima[2](s{degree}x{index}));
            let lambda{degree}x{index}=rvHilbCohQ[1,8](rvHilbSquareNakajima[8](rvHilbSquareLambda(hs,s{degree}x{index})));
            assert rvHilbCohLower(h,q{degree}x{index})==lambda{degree}x{index};
        "));
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn primitive_strings_have_the_sl2_bracket_and_rational_scaling() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let h=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{{(1,0,2),(1,1,1)}});
        let a=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{{(1,0,2),(1,2,1)}});
        let p=add(h,rvHilbCohScale[-2,1](a));
        assert rvHilbCohLower(h,p)==zero[RVHilbCoh(B,3,0)];
        let hp=rvHilbCohCup(h,p);
        let hhp=rvHilbCohCup(h,hp);
        let lowered=rvHilbCohLower(h,hhp);
        assert lowered==rvHilbCohScale[6,1](hp);
        let scaledH=rvHilbCohScale[-2,3](h);
        assert rvHilbCohLower(scaledH,hhp)==rvHilbCohScale[-3,2](lowered);
        let first=rvHilbCohCup(h,rvHilbCohLower(h,hp));
        let second=rvHilbCohLower(h,hhp);
        assert add(first,neg(second))==rvHilbCohScale[-2,1](hp);
    "
        ))),
        ""
    );
}

#[test]
fn scalar_inverse_and_degenerate_types_avoid_graded_matrix_work() {
    let n = 1000000000000000000i128;
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let u=rvHilbCohUnit[B,0];
        let z2=zero[RVHilbCoh(B,0,2)];
        assert rvHilbCohLefschetzInverse(z2,u)==u;
        assert rvHilbCohLower(z2,u)==zero[RVHilbCoh(B,0,-2)];
        let h=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,0,1),(1,1,1)}});
        let cup2=rvHilbCohCup(h,h);
        let cup4=rvHilbCohCup(cup2,cup2);
        assert rvHilbCohLefschetzInverse(h,cup4)==rvHilbCohUnit[B,2];
        let hugeH=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,0,{repeats}),(1,1,1)}});
        let hugeA=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,0,{repeats}),(1,2,1)}});
        let hugePrimitive=add(hugeH,rvHilbCohScale[-2,1](hugeA));
        assert rvHilbCohLower(hugeH,hugePrimitive)==zero[RVHilbCoh(B,{n},0)];
        let middle=zero[RVHilbCoh(B,{n},{twice})];
        assert rvHilbCohLefschetzInverse(hugeH,middle)==middle;
        let hugeTop=zero[RVHilbCoh(B,{n},{top})];
        assert rvHilbCohLefschetzInverse(hugeH,hugeTop)==zero[RVHilbCoh(B,{n},0)];
        assert rvHilbCohLower(h,zero[RVHilbCoh(B,2,7)])==zero[RVHilbCoh(B,2,5)];
    ",
            repeats = n - 1,
            twice = 2 * n,
            top = 4 * n
        ))),
        ""
    );
}

#[test]
fn isotropic_maps_are_singular_even_for_zero_right_hand_sides() {
    let setup = r"
        let a=rvK3Lattice[B]{[0,2,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
        let c=rvHilbCohChern[2,1](a);
        let d=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{(2,0,1)});
        let h=add(c,d);
    ";
    for expression in [
        "rvHilbCohLefschetzInverse(h,zero[RVHilbCoh(B,2,6)])",
        "rvHilbCohLefschetzInverse(h,rvHilbCohPoint[B,2])",
        "rvHilbCohLefschetzInverse(h,zero[RVHilbCoh(B,2,8)])",
        "rvHilbCohLower(h,h)",
    ] {
        let error = diagnostic(interpret(&format!("{MARKED}{setup} show {expression};")));
        assert!(
            error.contains("Hilbert Lefschetz cup map is singular"),
            "{error}"
        );
    }
}

#[test]
fn wrong_degrees_and_markings_have_specific_diagnostics() {
    for (input, expected) in [
        ("show rvHilbCohLefschetzInverse(zero[RVHilbCoh(B,2,2)],rvHilbCohUnit[B,2]);", "requires a middle or upper cohomological degree"),
        ("show rvHilbCohLower(rvHilbCohUnit[B,2],rvHilbCohUnit[B,2]);", "requires a degree-two divisor"),
        ("rvjactorsionbasis H[C,2](neg(o)){e1,e2,f1,f2}; show rvHilbCohLower(zero[RVHilbCoh(B,2,2)],rvHilbCohPoint[H,2]);", "same named level basis and length"),
    ] {
        let error=diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(error.contains(expected), "{error}");
    }
}
