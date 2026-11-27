mod common;

use common::{diagnostic, example, interpret, stdout};

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let a=rvUnramElem[F]{3};
    let h=rvUnramElem[F]{4};
    let e=rvUnramElem[F]{5};
    let s=rvUnramElem[F]{6};
    rvhyperelliptic C[F]{z,h,t,e,t,o};
    rvrichelot D[C]{[4],[2],[1]}{[],[1]}{[1],[],[1]};
    rvjacchart N[D](t);
";

const F9: &str = r"
    rvunram U[3,1]{1,0,1};
    let z=rvUnramElem[U]{};
    let o=rvUnramElem[U]{1};
    let alpha=rvUnramElem[U]{0,1};
    let b=rvUnramElem[U]{1,1};
    let c=rvUnramElem[U]{2,1};
    rvhyperelliptic C[U]{z,b,z,c,z,o};
    rvrichelot D[C]{[1],[],[1]}{[1,1],[],[1]}{[],[1]};
    rvrichelotdual E[D];
";

#[test]
fn charted_twist_only_image_has_the_same_coordinates_after_a_quadratic_lift() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-richelot")),
        "RVKummer(N)[[1], [2], [1], [2]]\n"
    );
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            rvunram W[7,1]{{4,0,1}};
            rvjacbase L[C,W];
            rvrichelot J[L]{{[4],[2],[1]}}{{[],[1]}}{{[1],[],[1]}};
            let two=rvUnramElem[W]{{2}};
            let one=rvUnramElem[W]{{1}};
            rvjacchart M[J](two);
            let lift=rvRichelotPoint[M](rvUnramElem[W]{{3}},rvUnramElem[W]{{0,1}});
            assert rvJacKummerPoint(lift)==rvKummer[M](one,two,one,two);
            let point=rvKummer[C](z,o,a,t);
            assert rvKummerRichelot[N](point)==rvKummer[N](o,t,o,t);
            "
        ))),
        ""
    );
}

#[test]
fn charted_map_commutes_with_cantor_on_the_curve_and_its_whole_twist_group() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            assert (\p:RVJac(C)=>rvKummerRichelot[N](rvJacKummerPoint(p)))
                ==(\p:RVJac(C)=>rvJacKummerPoint(rvRichelotPush[N](p)));
            rvhyperelliptic T[F]{{z,o,e,s,a,o}};
            rvrichelot H[T]{{[2],[3],[1]}}{{[],[1]}}{{[4],[],[1]}};
            rvjacchart M[H](a);
            let sourceBack=\u:RVKummer(T)=>rvKummer[C](rvKummerCoord[0](u),rvUnramMul(a,rvKummerCoord[1](u)),rvUnramMul(t,rvKummerCoord[2](u)),neg(rvKummerCoord[3](u)));
            let targetBack=\u:RVKummer(M)=>rvKummer[N](rvKummerCoord[0](u),neg(rvKummerCoord[1](u)),rvKummerCoord[2](u),neg(rvKummerCoord[3](u)));
            assert (\p:RVJac(T)=>rvKummerRichelot[N](sourceBack(rvJacKummerPoint(p))))
                ==(\p:RVJac(T)=>targetBack(rvJacKummerPoint(rvRichelotPush[M](p))));
            "
        ))),
        ""
    );
}

#[test]
fn kernels_doubling_and_differential_addition_descend_through_the_map() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            let project=\u:RVJac(C)=>rvJacKummerPoint(u);
            let push=\u:RVKummer(C)=>rvKummerRichelot[N](u);
            let id=zero[RVKummer(N)];
            let first=rvJacPoint[C](z,z);
            let second=add(rvJacPoint[C](o,z),rvJacPoint[C](h,z));
            assert push(zero[RVKummer(C)])==id;
            assert push(project(first))==id;
            assert push(project(second))==id;
            assert push(project(add(first,second)))==id;
            let p=rvKummer[C](z,o,a,t);
            let p2=rvKummerScale[2](p);
            let p3=rvKummerScale[3](p);
            assert push(p2)==rvKummerScale[2](push(p));
            assert push(p3)==rvKummerDiffAdd(push(p2),push(p),push(p));
            assert push(p3)==push(rvKummerDiffAdd(p2,p,p));
            "
        ))),
        ""
    );
}

#[test]
fn direct_target_dual_composes_to_doubling_on_the_entire_quartic() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F9}
            let special=rvKummer[C](z,o,o,o);
            assert rvKummerLiftable(special)==S.lo;
            let identify=\u:RVKummer(C)=>rvKummer[E](rvKummerCoord[0](u),rvKummerCoord[1](u),rvKummerCoord[2](u),rvKummerCoord[3](u));
            assert (\u:RVKummer(C)=>rvKummerRichelot[E](rvKummerRichelot[D](u)))
                ==(\u:RVKummer(C)=>identify(rvKummerScale[2](u)));
            rvjacchart N[D](o);
            assert (\p:RVJac(C)=>rvKummerRichelot[N](rvJacKummerPoint(p)))
                ==(\p:RVJac(C)=>rvJacKummerPoint(rvJacRechart[N](rvRichelotPush[D](p))));

            let k=rvUnramElem[U]{{1,2}};
            let t0=rvUnramElem[U]{{2,2}};
            let minusAlpha=rvUnramElem[U]{{0,2}};
            rvhyperelliptic T[U]{{z,t0,z,b,z,o}};
            rvrichelot H[T]{{[0,2],[],[1]}}{{[1,2],[],[1]}}{{[],[1]}};
            let q=rvJacPoint[T](t0,minusAlpha);
            let u=rvJacKummerPoint(rvRichelotPush[H](q));
            let expected=rvKummer[D](rvKummerCoord[0](u),rvUnramMul(k,rvKummerCoord[1](u)),rvUnramMul(alpha,rvKummerCoord[2](u)),rvUnramMul(b,rvKummerCoord[3](u)));
            assert rvKummerRichelot[D](special)==expected;
            "
        ))),
        ""
    );
}

#[test]
fn rejects_the_wrong_nominal_source_and_uncharted_nonmonic_targets() {
    for (source, expected) in [
        (
            "show rvKummerRichelot[D](zero[RVKummer(C)]);",
            "monic odd-degree",
        ),
        (
            "show rvKummerRichelot[C](zero[RVKummer(C)]);",
            "has no stored Richelot splitting",
        ),
        (
            "show rvKummerRichelot[N](zero[RVJac(C)]);",
            "expected type RVKummer(C)",
        ),
        (
            "rvhyperelliptic Twin[F]{z,h,t,e,t,o}; show rvKummerRichelot[N](zero[RVKummer(Twin)]);",
            "expected type RVKummer(C)",
        ),
        (
            "assert zero[RVKummer(C)]==rvKummerRichelot[N](zero[RVKummer(C)]);",
            "expected type RVKummer(C)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F7}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
