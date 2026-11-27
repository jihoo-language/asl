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

const THREE: &str = r"
    let a=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{(1,0,2),(1,1,1)});
    let b=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{(1,0,2),(1,2,1)});
    let c=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{(1,0,2),(1,3,1)});
    let delta=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{(1,0,1),(2,0,1)});
";

#[test]
fn example_rotates_a_class_outside_the_divisor_subring() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-hilbert-llv")),
        "rvHilbCohBeauville = -1\nrvHilbCohBeauville = -4\n"
    );
}

#[test]
fn marked_beauville_scales_surface_and_boundary_correctly() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}{THREE}
        let a1=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,1,1)}});
        let b1=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,2,1)}});
        rvHilbCohBeauville(a1,b1);
        assert rvHilbCohLLVLower(a1,b1)==rvHilbCohScale[-2,1](rvHilbCohUnit[B,1]);
        assert rvHilbCohRotate(a1,b1,a1)==add(neg(a1),rvHilbCohScale[2,1](b1));
        rvHilbCohBeauville(rvHilbCohScale[2,3](a),rvHilbCohScale[-3,5](a));
        rvHilbCohBeauville(rvHilbCohScale[2,1](a),rvHilbCohScale[2,1](a));
        rvHilbCohBeauville(a,delta);
        rvHilbCohBeauville(delta,delta);
        let d4=rvHilbCohQ[1,4](rvNakajimaMonomial[B]{{(1,0,2),(2,0,1)}});
        rvHilbCohBeauville(d4,d4);
        assert rvHilbCohLLVLower(d4,d4)==rvHilbCohScale[-48,1](rvHilbCohUnit[B,4]);
    "
        ))),
        concat!(
            "rvHilbCohBeauville = -1\n",
            "rvHilbCohBeauville = 4/5\n",
            "rvHilbCohBeauville = -8\n",
            "rvHilbCohBeauville = 0\n",
            "rvHilbCohBeauville = -4\n",
            "rvHilbCohBeauville = -6\n"
        )
    );
}

#[test]
fn three_point_lowering_and_rotation_cover_the_full_grading() {
    let mut source = format!("{MARKED}{THREE} let h0=rvHilbCohUnit[B,3];");
    for j in 1..=6 {
        source.push_str(&format!(
            r"
            let h{j}=rvHilbCohCup(a,h{});
            assert rvHilbCohLLVLower(a,h{j})==rvHilbCohScale[{},1](h{});
        ",
            j - 1,
            -2 * j * (7 - j),
            j - 1
        ));
    }
    for (degree, left, right) in [
        (4, "(2,0,1),(1,1,1)", "(2,0,1),(1,2,1)"),
        (6, "(3,1,1)", "(3,2,1)"),
        (8, "(2,23,1),(1,1,1)", "(2,23,1),(1,2,1)"),
        (10, "(1,23,2),(1,1,1)", "(1,23,2),(1,2,1)"),
    ] {
        source.push_str(&format!(r"
            let x{degree}=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{{left}}});
            let y{degree}=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{{right}}});
            assert rvHilbCohRotate(a,b,x{degree})==add(neg(x{degree}),rvHilbCohScale[2,1](y{degree}));
        "));
    }
    source.push_str(
        r"
        assert rvHilbCohRotate(a,b,h0)==zero[RVHilbCoh(B,3,0)];
        assert rvHilbCohRotate(a,b,rvHilbCohPoint[B,3])==zero[RVHilbCoh(B,3,12)];
    ",
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn isotropic_lowering_matches_the_square_in_every_degree_and_anchor() {
    let mut source = format!(
        r"{MARKED}
        let l=rvK3Lattice[B]{{[0,2,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        let delta=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{{(2,0,1)}});
        let iso=add(rvHilbCohChern[2,1](l),delta);
        let isos=rvHilbSquareScale[1,2](rvNakajimaHilbSquare(rvHilbCohFock[2](iso)));
        let a=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,0,1),(1,1,1)}});
        let asq=rvNakajimaHilbSquare(rvHilbCohFock[1](a));
        let tilt=add(iso,delta);
        rvHilbCohBeauville(iso,iso);
        rvHilbCohBeauville(tilt,tilt);
    "
    );
    for (degree, index) in [(2, 5), (4, 0), (4, 20), (6, 1), (8, 0)] {
        source.push_str(&format!(r"
            let s{degree}x{index}=rvHilbSquareBasis[B,{degree},{index}];
            let x{degree}c{index}=rvHilbCohQ[1,2](rvHilbSquareNakajima[2](s{degree}x{index}));
            let f{degree}c{index}=rvHilbCohLLVLower(iso,x{degree}c{index});
            let via{degree}c{index}=rvHilbCohQ[1,2](rvHilbSquareNakajima[2](rvHilbSquareLower(isos,s{degree}x{index})));
            assert f{degree}c{index}==via{degree}c{index};
            let anchored{degree}c{index}=add(rvHilbCohLLVLower(tilt,x{degree}c{index}),neg(rvHilbCohLLVLower(delta,x{degree}c{index})));
            assert f{degree}c{index}==anchored{degree}c{index};
            let r{degree}c{index}=rvHilbCohQ[1,2](rvHilbSquareNakajima[2](rvHilbSquareRotate(asq,isos,s{degree}x{index})));
            assert rvHilbCohRotate(a,iso,x{degree}c{index})==r{degree}c{index};
        "));
    }
    assert_eq!(
        stdout(interpret(&source)),
        "rvHilbCohBeauville = 0\nrvHilbCohBeauville = -2\n"
    );
}

#[test]
fn isotropic_translation_retries_when_the_first_anchor_is_isotropic() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let l=rvK3Lattice[B]{{[0,0,0,0,0,-2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        let d=rvHilbCohQ[1,2](rvNakajimaMonomial[B]{{(2,0,1)}});
        let iso=add(rvHilbCohChern[2,1](l),d);
        let a=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,0,1),(1,1,1)}});
        let first=add(iso,a);
        let second=add(iso,rvHilbCohScale[2,1](a));
        rvHilbCohBeauville(first,first);
        rvHilbCohBeauville(second,second);
        let dSquared=rvHilbCohCup(d,d);
        let f=rvHilbCohLLVLower(iso,dSquared);
        assert f==rvHilbCohScale[-4,1](iso);
    "
        ))),
        "rvHilbCohBeauville = 0\nrvHilbCohBeauville = -4\n"
    );
}

#[test]
fn three_point_lowering_commutes_and_rotation_is_a_derivation() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}{THREE}
        let middle=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(2,23,1),(1,1,1)}});
        let fa=rvHilbCohLLVLower(a,middle);
        let fb=rvHilbCohLLVLower(b,middle);
        assert rvHilbCohLLVLower(a,fb)==rvHilbCohLLVLower(b,fa);
        let extra=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(2,0,1),(1,1,1)}});
        let by=rvHilbCohCup(b,extra);
        let left=rvHilbCohRotate(a,b,by);
        let rb=rvHilbCohRotate(a,b,b);
        let rx=rvHilbCohRotate(a,b,extra);
        let term1=rvHilbCohCup(rb,extra);
        let term2=rvHilbCohCup(b,rx);
        assert left==add(term1,term2);
        assert rx==neg(rvHilbCohRotate(b,a,extra));
        let scale=rvHilbCohScale[-2,3](a);
        assert rvHilbCohLLVLower(scale,middle)==rvHilbCohScale[-2,3](fa);
    "
        ))),
        ""
    );
}

#[test]
fn three_point_isotropic_linearity_and_the_rotation_lie_bracket() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}{THREE}
        let isotropicLine=rvK3Lattice[B]{{[0,2,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        let iso=add(rvHilbCohChern[3,1](isotropicLine),delta);
        let tilted=add(iso,delta);
        let x=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(2,0,1),(1,1,1)}});
        let fIso=rvHilbCohLLVLower(iso,x);
        assert fIso==add(rvHilbCohLLVLower(tilted,x),neg(rvHilbCohLLVLower(delta,x)));
        let aplus=add(iso,a);
        assert rvHilbCohLLVLower(aplus,x)==add(fIso,rvHilbCohLLVLower(a,x));
        let xac=rvHilbCohRotate(a,c,x);
        let xab=rvHilbCohRotate(a,b,x);
        let fbxa=rvHilbCohRotate(a,b,xac);
        let faxb=rvHilbCohRotate(a,c,xab);
        let aa=rvHilbCohRotate(a,b,a);
        let ac=rvHilbCohRotate(a,b,c);
        let y1=rvHilbCohRotate(aa,c,x);
        let y2=rvHilbCohRotate(a,ac,x);
        assert add(fbxa,neg(faxb))==add(y1,y2);
    "
        ))),
        ""
    );
}

#[test]
fn nominal_divisor_requirements_and_formal_zeros() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let z2=zero[RVHilbCoh(B,0,2)];
        let z4=zero[RVHilbCoh(B,0,4)];
        rvHilbCohBeauville(z2,z2);
        assert rvHilbCohLLVLower(z2,z4)==z2;
        assert rvHilbCohRotate(z2,z2,z4)==z4;
    "
        ))),
        "rvHilbCohBeauville = 0\n"
    );
    for (input,expected) in [
        ("rvHilbCohBeauville(rvHilbCohUnit[B,2],zero[RVHilbCoh(B,2,2)]);", "requires a degree-two divisor"),
        ("show rvHilbCohLLVLower(zero[RVHilbCoh(B,2,2)],rvHilbCohPoint[B,3]);", "same named level basis and length"),
        ("rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; show rvHilbCohRotate(zero[RVHilbCoh(B,2,2)],zero[RVHilbCoh(D,2,2)],rvHilbCohPoint[B,2]);", "same named level basis and length"),
    ] {
        let err=diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(err.contains(expected),"{err}");
    }
}
