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
    let a=rvHilbSquareDiv[0](rvK3Lattice[B]{[0,2,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]});
    let b=rvHilbSquareDiv[0](rvK3Lattice[B]{[0,0,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]});
    let delta=rvHilbSquareBasis[B,2,22];
    let c=add(rvHilbSquareScale[1,3](a),rvHilbSquareScale[-2,5](b));
    let x0=rvHilbSquareScale[3,7](rvHilbSquareBasis[B,0,0]);
    let x2=add(c,rvHilbSquareScale[4,9](delta));
    let x4=add(rvHilbSquareCup(x2,add(a,delta)),rvHilbSquareScale[-5,3](rvHilbSquareC2[B]));
    let x6=rvHilbSquareCup(x4,add(b,delta));
    let x8=rvHilbSquareScale[-7,11](rvHilbSquareBasis[B,8,0]);
";

#[test]
fn raising_and_lowering_brackets_hold_in_every_degree() {
    let mut source = MARKED.to_owned();
    source.push_str(
        r"
        assert rvHilbSquareLower(a,x0)==zero[RVHilbSquare(B,-2)];
        assert rvHilbSquareLower(a,zero[RVHilbSquare(B,-2)])==zero[RVHilbSquare(B,-4)];
        assert rvHilbSquareCup(a,zero[RVHilbSquare(B,-2)])==zero[RVHilbSquare(B,0)];
        assert rvHilbSquareRotate(a,b,a)==rvHilbSquareScale[2,1](a);
        assert rvHilbSquareLower(delta,delta)==rvHilbSquareScale[-8,1](rvHilbSquareBasis[B,0,0]);
        let formalOdd=zero[RVHilbSquare(B,1)];
        let formalHigh=zero[RVHilbSquare(B,10)];
        ",
    );
    for (degree, class) in [
        (0, "x0"),
        (2, "x2"),
        (4, "x4"),
        (6, "x6"),
        (8, "x8"),
        (1, "formalOdd"),
        (10, "formalHigh"),
    ] {
        let weight = degree - 4;
        writeln!(
            source,
            r"
            assert add(rvHilbSquareCup(a,rvHilbSquareLower(b,{class})),
                       neg(rvHilbSquareLower(b,rvHilbSquareCup(a,{class}))))
                 ==add(rvHilbSquareScale[2,1](rvHilbSquareRotate(a,b,{class})),
                        rvHilbSquareScale[{},1]({class}));
            assert add(rvHilbSquareCup(delta,rvHilbSquareLambda(delta,{class})),
                       neg(rvHilbSquareLambda(delta,rvHilbSquareCup(delta,{class}))))
                 ==rvHilbSquareScale[{weight},1]({class});
            assert rvHilbSquareLower(a,rvHilbSquareLower(b,{class}))
                 ==rvHilbSquareLower(b,rvHilbSquareLower(a,{class}));
            assert rvHilbSquareLower(c,{class})
                 ==add(rvHilbSquareScale[1,3](rvHilbSquareLower(a,{class})),
                       rvHilbSquareScale[-2,5](rvHilbSquareLower(b,{class})));
            assert rvHilbSquareLambda(c,{class})
                 ==rvHilbSquareScale[-15,8](rvHilbSquareLower(c,{class}));
            ",
            2 * weight
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn rotations_are_skew_and_act_on_raising_and_lowering_generators() {
    let mut source = MARKED.to_owned();
    source.push_str(
        r"
        let v=add(b,delta);
        let rc=rvHilbSquareRotate(a,b,c);
        let rv=rvHilbSquareRotate(a,b,v);
        assert rvHilbSquareRotate(a,b,rvHilbSquareC2[B])==zero[RVHilbSquare(B,4)];
        ",
    );
    for degree in [0, 2, 4, 6, 8] {
        writeln!(
            source,
            r"
            assert rvHilbSquareRotate(a,b,x{degree})==neg(rvHilbSquareRotate(b,a,x{degree}));
            assert add(rvHilbSquareRotate(a,b,rvHilbSquareCup(c,x{degree})),
                       neg(rvHilbSquareCup(c,rvHilbSquareRotate(a,b,x{degree}))))
                 ==rvHilbSquareCup(rvHilbSquareRotate(a,b,c),x{degree});
            assert add(rvHilbSquareRotate(a,b,rvHilbSquareLower(c,x{degree})),
                       neg(rvHilbSquareLower(c,rvHilbSquareRotate(a,b,x{degree}))))
                 ==rvHilbSquareLower(rc,x{degree});
            assert add(rvHilbSquareRotate(a,b,rvHilbSquareRotate(c,v,x{degree})),
                       neg(rvHilbSquareRotate(c,v,rvHilbSquareRotate(a,b,x{degree}))))
                 ==add(rvHilbSquareRotate(rc,v,x{degree}),rvHilbSquareRotate(c,rv,x{degree}));
            "
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn primitive_projectors_split_off_the_lefschetz_summands() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            assert rvHilbSquarePrimitive(delta,x0)==x0;
            assert rvHilbSquarePrimitive(delta,x2)==c;
            assert rvHilbSquarePrimitive(delta,c)==c;
            assert rvHilbSquareLower(delta,c)==zero[RVHilbSquare(B,0)];
            let a2=rvHilbSquareCup(a,a);
            let delta2=rvHilbSquareCup(delta,delta);
            let deltaB=rvHilbSquareCup(delta,b);
            let beta=add(a2,add(rvHilbSquareScale[5,7](delta2),rvHilbSquareScale[-3,2](deltaB)));
            assert rvHilbSquarePrimitive(delta,beta)==a2;
            assert rvHilbSquarePrimitive(delta,a2)==a2;
            assert rvHilbSquareCup(a2,rvHilbSquareCup(b,b))
                ==rvHilbSquareScale[8,1](rvHilbSquareBasis[B,8,0]);
            assert rvHilbSquareCup(delta,a2)==zero[RVHilbSquare(B,6)];
            assert rvHilbSquareLower(delta,a2)==zero[RVHilbSquare(B,2)];
            let p=rvHilbSquarePrimitive(x2,x4);
            assert rvHilbSquarePrimitive(x2,p)==p;
            assert rvHilbSquareCup(x2,p)==zero[RVHilbSquare(B,6)];
            let pp=rvHilbSquarePrimitive(x2,c);
            assert rvHilbSquarePrimitive(x2,pp)==pp;
            assert rvHilbSquareLower(x2,pp)==zero[RVHilbSquare(B,0)];
            "
        ))),
        ""
    );
}

#[test]
fn llv_errors_preserve_the_marking_isotropic_and_signed_degree_rules() {
    for (source, expected) in [
        ("show rvHilbSquareLambda(a,x0);", "requires a nonisotropic generator"),
        ("show rvHilbSquarePrimitive(a,x4);", "requires a nonisotropic generator"),
        ("show rvHilbSquarePrimitive(delta,x6);", "requires cohomological degree zero, two, or four"),
        ("show rvHilbSquareLower(x4,x0);", "requires a generator of cohomological degree two"),
        ("show rvHilbSquareRotate(a,x0,x4);", "requires a generator of cohomological degree two"),
        ("show rvHilbSquareRotate(a,b,z);", "an RVHilbSquare cohomology class"),
        ("rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; show rvHilbSquareLower(a,zero[RVHilbSquare(D,0)]);", "requires the same named marking"),
        ("rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; show rvHilbSquareRotate(a,zero[RVHilbSquare(D,2)],x4);", "requires the same named marking"),
        ("show rvHilbSquareLower(delta,zero[RVHilbSquare(B,-170141183460469231731687303715884105728)]);", "lowering degree exceeds the signed 128-bit range"),
        ("show rvHilbSquareCup(delta,zero[RVHilbSquare(B,170141183460469231731687303715884105727)]);", "cup product degree exceeds the signed 128-bit range"),
    ] {
        let error = diagnostic(interpret(&format!("{MARKED}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn example_integrates_the_lowered_point_after_cupping_with_delta() {
    assert_eq!(
        stdout(example("hyperelliptic-hilbert-square-llv")),
        "rvHilbSquareIntegral = -8\nrvHilbSquareIntegral = 4\n"
    );
}
