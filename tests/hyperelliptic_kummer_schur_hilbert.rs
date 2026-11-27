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
fn example_exhibits_exterior_twists_and_the_weyl_rank() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-schur-hilbert")),
        "rvHilbExteriorEuler[2, 2] = 0\nrvHilbExteriorEuler[2, 1] = -2\n"
    );
}

#[test]
fn axes_vanishing_and_large_exact_weyl_dimension() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let a=rvK3LatticeBasis[B,0];
        assert rvHilbCohSchur[0,0]{{}}(a)==rvHilbCohUnit[B,0];
        assert rvHilbCohSchur[0,0]{{1}}(a)==zero[RVHilbCoh(B,0,0)];
        assert rvHilbCohSchur[3,0]{{0,0}}(a)==rvHilbCohUnit[B,3];
        assert rvHilbCohSchur[3,1]{{}}(a)==zero[RVHilbCoh(B,3,2)];
        assert rvHilbCohSchur[2,1]{{2,1,1}}(a)==zero[RVHilbCoh(B,2,2)];
        assert rvHilbCohSchur[2,5]{{2}}(a)==zero[RVHilbCoh(B,2,10)];
        assert rvHilbCohSchur[3,2]{{1}}(a)==rvHilbCohCh[3,2](a);
        assert rvHilbCohSchur[4,0]{{2,1,1}}(a)==rvHilbCohScale[15,1](rvHilbCohUnit[B,4]);
        assert rvHilbCohSchur[1,1]{{1000000000000}}(a)==rvHilbCohScale[1000000000000,1](rvHilbCohCh[1,1](a));
        assert rvHilbCohSchur[1,2]{{3}}(a)==rvHilbCohScale[-9,1](rvHilbCohPoint[B,1]);
        let veryLarge=rvHilbCohSchur[2,0]{{{}}}(a);
        assert veryLarge==rvHilbCohQ[{},1](rvNakajimaMonomial[B]{{(1,0,2)}});
    ",
            u128::MAX - 2,
            i128::MAX
        ))),
        ""
    );
}

#[test]
fn determinant_twists_and_nontrivial_high_cohomology() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let a=rvK3LatticeBasis[B,0];
        let c1=rvHilbCohChern[2,1](a);
        assert rvHilbCohSchur[2,1]{{5,4}}(a)==rvHilbCohScale[9,1](c1);
        let q2=rvHilbCohCup(c1,c1);
        let q4=rvHilbCohCup(q2,q2);
        let det=rvHilbCohSchur[2,4]{{1,1}}(a);
        assert det==rvHilbCohScale[1,24](q4);
        rvHilbCohIntegral(det);
    "
        ))),
        "rvHilbCohIntegral = 2\n"
    );
}

#[test]
fn complete_and_dual_jacobi_trudi_satisfy_pieri_in_positive_degree() {
    let mut source = format!("{MARKED} let a=rvK3LatticeBasis[B,0];");
    for (length, left, first, second) in [(3, "2", "3", "2,1"), (4, "1,1,1", "2,1,1", "1,1,1,1")] {
        source.push_str(&format!(r"
            let p0n{length}=rvHilbCohCup(rvHilbCohSchur[{length},0]{{1}}(a),rvHilbCohSchur[{length},2]{{{left}}}(a));
            let p1n{length}=rvHilbCohCup(rvHilbCohSchur[{length},1]{{1}}(a),rvHilbCohSchur[{length},1]{{{left}}}(a));
            let p2n{length}=rvHilbCohCup(rvHilbCohSchur[{length},2]{{1}}(a),rvHilbCohSchur[{length},0]{{{left}}}(a));
            let rhs{length}=add(rvHilbCohSchur[{length},2]{{{first}}}(a),rvHilbCohSchur[{length},2]{{{second}}}(a));
            assert add(add(p0n{length},p1n{length}),p2n{length})==rhs{length};
        "));
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn exterior_euler_agrees_with_hilbert_square_todd_even_with_a_twist() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let a=rvK3LatticeBasis[B,0];
        let trivial=zero[RVK3Lattice(B)];
        let point=rvHilbCohPoint[B,2];
        let td2=rvHilbCohQ[1,12](rvHilbSquareNakajima[1](rvHilbSquareC2[B]));
        let det2=rvHilbCohSchur[2,2]{{1,1}}(a);
        let det4=rvHilbCohSchur[2,4]{{1,1}}(a);
        let integralTop=add(add(det4,rvHilbCohCup(det2,td2)),rvHilbCohScale[3,1](point));
        assert integralTop==zero[RVHilbCoh(B,2,8)];
        rvHilbCohIntegral(integralTop);
        rvHilbExteriorEuler[2,2](a,trivial);
        let divisor=rvHilbCohQ[1,1](rvNakajimaMonomial[B]{{(1,0,1),(1,1,1)}});
        let d2=rvHilbCohCup(divisor,divisor);
        let d3=rvHilbCohCup(d2,divisor);
        let d4=rvHilbCohCup(d2,d2);
        let s1=rvHilbCohSchur[2,1]{{1}}(a);
        let s2=rvHilbCohSchur[2,2]{{1}}(a);
        let s3=rvHilbCohSchur[2,3]{{1}}(a);
        let s4=rvHilbCohSchur[2,4]{{1}}(a);
        let twist2=add(add(s2,rvHilbCohCup(s1,divisor)),d2);
        let term2=rvHilbCohScale[1,2](rvHilbCohCup(s2,d2));
        let term3=rvHilbCohScale[1,6](rvHilbCohCup(s1,d3));
        let term4=rvHilbCohScale[1,12](d4);
        let twist4=add(add(s4,rvHilbCohCup(s3,divisor)),add(add(term2,term3),term4));
        let twistTop=add(add(twist4,rvHilbCohCup(twist2,td2)),rvHilbCohScale[6,1](point));
        assert twistTop==rvHilbCohScale[-2,1](point);
        rvHilbCohIntegral(twistTop);
        rvHilbExteriorEuler[2,1](a,a);
    "
        ))),
        concat!(
            "rvHilbCohIntegral = 0\n",
            "rvHilbExteriorEuler[2, 2] = 0\n",
            "rvHilbCohIntegral = -2\n",
            "rvHilbExteriorEuler[2, 1] = -2\n"
        )
    );
}

#[test]
fn exterior_euler_uses_generalized_binomials_and_sparse_extremes() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let a=rvK3LatticeBasis[B,0];
        let trivial=zero[RVK3Lattice(B)];
        let chiZero=rvK3Lattice[B]{{[0,0,0,0,0,0],[2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        rvHilbExteriorEuler[0,0](a,a);
        rvHilbExteriorEuler[0,1](a,a);
        rvHilbExteriorEuler[2,2](a,a);
        rvHilbExteriorEuler[1000000000000000000000,0](a,trivial);
        rvHilbExteriorEuler[1000000000000000000000,0](a,chiZero);
    "
        ))),
        concat!(
            "rvHilbExteriorEuler[0, 0] = 1\n",
            "rvHilbExteriorEuler[0, 1] = 0\n",
            "rvHilbExteriorEuler[2, 2] = 3\n",
            "rvHilbExteriorEuler[1000000000000000000000, 0] = 1000000000000000000001\n",
            "rvHilbExteriorEuler[1000000000000000000000, 0] = 0\n"
        )
    );
}

#[test]
fn shape_budget_and_nominal_errors_are_reported() {
    for (input, expected) in [
        ("show rvHilbCohSchur[2,1]{1,2}(rvK3LatticeBasis[B,0]);", "Schur rows must be weakly decreasing"),
        ("show rvHilbCohSchur[4,2]{1000000000000000000}(rvK3LatticeBasis[B,0]);", "shape and degree exceed the computation limit"),
        ("show rvHilbCohSchur[2,0]{340282366920938463463374607431768211455}(rvK3LatticeBasis[B,0]);", "coefficient exceeds the signed 128-bit range"),
        ("rvjactorsionbasis H[C,2](neg(o)){e1,e2,f1,f2}; rvHilbExteriorEuler[2,1](rvK3LatticeBasis[B,0],rvK3LatticeBasis[H,0]);", "rvHilbExteriorEuler twist"),
    ] {
        let error=diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(error.contains(expected), "{error}");
    }
}
