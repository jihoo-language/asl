mod common;

use common::{diagnostic, example, interpret, stdout};

const F5: &str = r"
    rvunram F[5,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let three=rvUnramElem[F]{3};
    let four=rvUnramElem[F]{4};
    rvhyperelliptic C[F]{o,o,z,three,o,o};
    let p=rvJacPoint[C](z,o);
    let r=rvJacPoint[C](t,z);
    let s=rvJacPoint[C](three,t);
";

#[test]
fn nonsplit_trace_matches_frobenius_and_descent() {
    assert_eq!(
        stdout(example("hyperelliptic-quadratic-trace")),
        "RVJac(C){u=[[2], [2], [1]], v=[[1], [2]]}\n"
    );
}

#[test]
fn split_trace_handles_shared_support_both_signs_and_ramification() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            assert rvJacQuadTrace[C](o){{[1]}}{{}}{{[2]}}{{[3]}}==zero[RVJac(C)];
            assert rvJacQuadTrace[C](o){{[1],[1]}}{{[4]}}{{[4]}}{{[2]}}==add(p,s);
            assert rvJacQuadTrace[C](four){{[1],[1]}}{{[2]}}{{[4]}}{{[1]}}==add(p,s);
            assert rvJacQuadTrace[C](o){{[1],[1]}}{{[4]}}{{[4],[],[4]}}{{[2],[],[1],[1]}}==add(p,s);
            assert rvJacQuadTrace[C](o){{[],[1]}}{{}}{{[1]}}{{}}==add(p,p);
            assert rvJacQuadTrace[C](o){{[],[1]}}{{}}{{}}{{[1]}}==zero[RVJac(C)];
            assert rvJacQuadTrace[C](o){{[3],[1]}}{{}}{{}}{{}}==zero[RVJac(C)];

            assert rvJacQuadTrace[C](o){{[],[],[1]}}{{[],[3]}}{{[1],[2]}}{{}}
                ==add(add(p,r),add(p,s));
            assert rvJacQuadTrace[C](o){{[],[],[1]}}{{[],[3]}}{{[],[4]}}{{[1],[3]}}
                ==add(r,s);
        "
        ))),
        ""
    );
}

#[test]
fn repeated_support_survives_tracing_and_shared_smaller_multiplicities() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            let pp=add(p,p);
            assert rvJacQuadTrace[C](o){{[],[],[1]}}{{}}{{[1],[3]}}{{}}==add(pp,pp);
            assert rvJacQuadTrace[C](o){{[],[],[1]}}{{}}{{}}{{[1],[3]}}==zero[RVJac(C)];
            assert rvJacQuadTrace[C](o){{[],[1],[1]}}{{[],[4]}}{{[1]}}{{[],[3]}}
                ==add(pp,add(p,s));
            assert rvJacQuadTrace[C](o){{[],[1],[1]}}{{[],[4]}}{{[],[2]}}{{[1],[1]}}
                ==add(p,s);
            assert rvJacQuadTrace[C](o){{[],[1],[],[1]}}{{}}{{[1],[2]}}{{}}
                ==rvJacScale[2](add(add(p,r),s));
        "
        ))),
        ""
    );
}

#[test]
fn nonsplit_trace_uses_nonprime_coefficients_and_radicand() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{1,0,1};
            rvunram L[3,1]{2,0,1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=rvUnramElem[F]{2};
            let a=rvUnramElem[F]{0,1};
            rvhyperelliptic C[F]{z,o,o,o,t,o};
            rvjacbase D[C,L];
            let p=rvJacPoint[D](rvUnramElem[L]{0,1},rvUnramElem[L]{2,1,1,2});
            let norm=add(p,rvJacFrob[C,1](p));
            let trace=rvJacQuadTrace[C](add(o,a)){[],[1]}{[2]}{[0,1]}{[0,2]};
            assert trace==rvJacDescend[C](norm);
            assert rvJacExtend[D](trace)==norm;
        "
        )),
        ""
    );
}

#[test]
fn trace_requires_etale_radicand_valid_relative_pairs_and_nominal_field() {
    for (source, expected) in [
        (
            "show rvJacQuadTrace[C](z){[1]}{}{}{};",
            "radicand must be nonzero",
        ),
        (
            "show rvJacQuadTrace[C](o){}{}{}{};",
            "U must be nonzero and monic",
        ),
        (
            "show rvJacQuadTrace[C](o){[],[1]}{[],[1]}{}{};",
            "U must be nonzero and monic",
        ),
        (
            "show rvJacQuadTrace[C](o){[],[1]}{}{[3]}{[3]};",
            "U must divide f",
        ),
        (
            "show rvJacQuadTrace[C](o){[4],[1],[1]}{}{}{};",
            "U must divide f",
        ),
        (
            "rvunram G[5,1]{0,1}; show rvJacQuadTrace[C](rvUnramElem[G]{1}){[1]}{}{}{};",
            "expected type RVUnram(F)",
        ),
        (
            "rvhyperelliptic E[F]{o,o,o}; show rvJacQuadTrace[E](o){[1]}{}{}{};",
            "monic odd-degree",
        ),
        (
            "show rvJacQuadTrace[C](o){[1,0]}{}{}{};",
            "at most 1 coordinates",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F5}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
