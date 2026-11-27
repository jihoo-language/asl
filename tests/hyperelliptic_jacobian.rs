mod common;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{0};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,o,o,o,t,o};
";

#[test]
fn genuine_frobenius_relation_uses_nonbase_class_of_order_194() {
    assert_eq!(
        stdout(example("hyperelliptic-jacobian")),
        concat!(
            "rvHWWWeil[C] = [1, 0, -2, 0, 9]\n",
            "rvJacCount[C] = 8\n",
            "rvJacEnumerate[C] = 8\n",
            "rvJacOrder = 2\n",
            "rvJacCount[L] = 776\n",
            "rvJacOrder = 194\n",
            "RVJac(L){u=[[1, 0, 0]], v=[]}\n",
        )
    );
}

#[test]
fn degree_three_reduction_and_independent_enumeration() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
        rvhyperelliptic G[F]{{o,o,z,z,z,z,z,o}};
        rvJacCount[G];
        rvJacEnumerate[G];
        let p=rvJacPoint[G](z,o);
        let q=rvJacPoint[G](o,z);
        assert rvJacV[0](p)==o;
        assert rvJacU[1](p)==o;
        assert rvJacScale[-5](p)==neg(rvJacScale[5](p));
        assert add(add(p,p),q)==add(p,add(q,p));
        assert add(p,neg(p))==zero[RVJac(G)];
        assert add(p,zero[RVJac(G)])==p;
        assert rvJacWeilApply[G](q)==zero[RVJac(G)];
    "
        ))),
        "rvJacCount[G] = 34\nrvJacEnumerate[G] = 34\n"
    );
}

const F9: &str = r"
    rvunram U[3,1]{1,0,1};
    rvunram V[3,1]{2,1,0,0,1};
    let a=rvUnramElem[U]{0,1};
    let b=rvUnramElem[U]{2,1};
    let o=rvUnramElem[U]{1};
    let z=rvUnramElem[U]{};
    rvhyperelliptic C[U]{b,o,z,z,z,o};
    rvjacbase D[C,V];
    let source=rvJacPoint[C](a,a);
    let alpha=rvUnramElem[V]{0,2,1,1};
    let point=rvJacPoint[D](alpha,alpha);
    let q=rvJacPoint[D](rvUnramElem[V]{2,1},rvUnramElem[V]{0,2,0,1});
";

#[test]
fn relative_frobenius_and_descent_on_nonprime_fields() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F9}
        assert rvJacExtend[D](source) == point;
        assert rvJacU[0](point) == neg(alpha);
        assert rvJacDescend[C](point) == source;
        assert rvJacDescend[C](add(point,point)) == add(source,source);
        assert rvJacFrob[C,1](point) == point;
        let conjugate=rvJacPoint[D](rvUnramElem[V]{{2,1,1,1}},rvUnramElem[V]{{0,0,1}});
        assert rvJacFrob[C,1](q) == conjugate;
        assert rvJacFrob[C,2](q) == q;
        let trace=rvJacDescend[C](add(q,conjugate));
        assert rvJacExtend[D](trace)==add(q,conjugate);
        assert rvJacWeilApply[C](q)==zero[RVJac(D)];
        let abar=rvUnramFrob[1](a);
        let bbar=rvUnramFrob[1](b);
        rvhyperelliptic Cbar[U]{{bbar,o,z,z,z,o}};
        rvjacbase Dbar[Cbar,V];
        let dbar=rvJacExtend[Dbar](rvJacPoint[Cbar](abar,abar));
        assert rvJacV[0](dbar)==rvUnramFrob[1](rvJacV[0](point));
        assert rvJacU[0](dbar)==rvUnramFrob[1](rvJacU[0](point));
    "
        ))),
        ""
    );
    assert!(
        diagnostic(interpret(&format!("{F9} show rvJacDescend[C](q);")))
            .contains("outside the requested subfield")
    );
}

#[test]
fn rejects_bad_mumford_pairs_and_nominal_endpoints() {
    for (source, expected) in [
        ("show rvjac[C]{[1],[2]}{};", "nonzero and monic"),
        ("show rvjac[C]{[],[1]}{[1]};", "must divide"),
        ("show rvJacPoint[C](z,o);", "must divide"),
        ("rvhyperelliptic D[F]{z,o,o,o,t,o}; show add(zero[RVJac(C)],zero[RVJac(D)]);", "expected type RVJac(C)"),
        ("rvhyperelliptic D[F]{o,o,o,o,t,o}; show rvJacExtend[D](zero[RVJac(C)]);", "not the declared base extension"),
        ("rvhyperelliptic D[F]{o,o,z,z,z,t}; show zero[RVJac(D)];", "monic odd-degree"),
        ("rvunram W[3,2]{0,1}; let w=rvUnramElem[W]{1}; rvhyperelliptic D[W]{w,w}; show zero[RVJac(D)];", "precision one"),
    ] {
        let problem=diagnostic(interpret(&format!("{F3}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}

#[test]
fn finite_functions_index_jacobian_values_in_canonical_order() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
        let fold=\d:RVJac(C)=>add(d,d);
        assert fold==(\d:RVJac(C)=>zero[RVJac(C)]);
        assert fold(rvJacPoint[C](o,z))==zero[RVJac(C)];
        let self=\d:RVJac(C)=>d;
        assert self(rvJacPoint[C](t,z))==rvJacPoint[C](t,z);
    "
        ))),
        ""
    );
}

#[test]
fn semireduced_composition_exposes_the_exact_norm_quotient() {
    assert_eq!(
        stdout(interpret(
            r"
        rvunram F[5,1]{0,1};
        let one=rvUnramElem[F]{1};
        let zero=rvUnramElem[F]{};
        rvhyperelliptic C[F]{one,one,zero,one};
        let p=rvJacPoint[C](zero,one);
        rvJacCompose(p,p);
        let semi=rvjac[C]{[],[],[1]}{[1],[3]};
        assert semi==add(p,p);
        assert semi==rvjac[C]{[1],[1]}{[2]};
    "
        )),
        "rvJacCompose = u=[[0], [0], [1]], v=[[1], [3]], norm quotient=[[1], [1]]\n"
    );
}
