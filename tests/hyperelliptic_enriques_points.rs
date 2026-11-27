mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/hyperelliptic-enriques-points.asl");

fn preamble() -> &'static str {
    EXAMPLE.split_once("\ncard ").unwrap().0
}

#[test]
fn all_rational_orbits_are_enumerated_and_quotiented_by_the_actual_switch() {
    assert_eq!(
        stdout(example("hyperelliptic-enriques-points")),
        concat!(
            "card = 22 / 22\n",
            "card = 16 / 22\n",
            "RVEnriques(Y){split, base=[[0, 0], [0, 0], [0, 0], [1, 0]], direction=[[0, 0], [0, 0], [1, 0], [0, 0]]}\n",
            "RVEnriques(Y){nonsplit, base=[[1, 0], [1, 0], [2, 1], [1, 2]]}\n",
        )
    );
}

#[test]
fn a_cubic_marking_translates_over_a_separate_quadratic_cover() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            rvunram W[3,1]{1,2,0,1};
            rvunram V[3,1]{1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=neg(o);
            rvhyperelliptic C[F]{z,t,t,t,t,o};
            rvjacbase L[C,W];
            rvjacbase Q[C,V];
            let zw=rvUnramElem[W]{};
            let ow=rvUnramElem[W]{1};
            let a=rvUnramElem[W]{0,1};
            let e1=rvJacPoint[L](zw,zw);
            let f1=rvJacPoint[L](ow,zw);
            let e2=add(add(e1,f1),rvJacPoint[L](a,zw));
            let f2=add(add(e1,f1),rvJacPoint[L](rvUnramFrob[1](a),zw));
            rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};
            rvenriques Y[C,Q,B,5];
            rvKummerEnriquesPoints[B,C,5,1];
            card (\y:RVEnriques(Y)=>S.hi);
            card (\y:RVEnriques(Y)=>rvEnriquesSplit(y));
            let x=rvKummerExceptional(zero[RVKummer(C)],z,z,o,z);
            let y=rvEnriquesOrbit[Y](x);
            assert rvEnriquesLift(y)==x;
            assert rvEnriquesRepresentative(y)==rvKummerK3Extend[Q](x);
            assert (\u:RVEnriques(Y)=>rvEnriquesOrbit[Y](rvEnriquesRepresentative(u)))==(\u:RVEnriques(Y)=>u);
        "
        )),
        concat!(
            "rvKummerEnriquesPoints[B, C, 5, 1] = points=13, normalized trace=1, split=12, nonsplit=1, twisted cover=2\n",
            "card = 13 / 13\n",
            "card = 12 / 13\n",
        )
    );
}

#[test]
fn rejects_nonsplit_lifts_and_nonstable_cover_points() {
    for (statement, expected) in [
        (
            "let x=rvKummerSmooth(rvKummer[L](ow,ow,rvUnramElem[W]{2,1},rvUnramElem[W]{1,2})); show rvEnriquesLift(rvEnriquesOrbit[Y](x));",
            "nonsplit Enriques orbit has no rational lift",
        ),
        (
            "let x=rvKummerTropeParam[B,0](ow,rvUnramElem[W]{1,1}); show rvEnriquesOrbit[Y](x);",
            "does not give a Frobenius-stable Enriques orbit",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{}\n{statement}", preamble())));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn checks_masks_models_and_nominal_quotient_types() {
    for (statement, expected) in [
        ("rvenriques Z[C,L,B,13];", "requires an even theta characteristic"),
        ("rvenriques Z[C,L,B,8];", "must be fixed by base Frobenius"),
        ("rvenriques Z[C,L,B,16];", "between 0 and 15"),
        ("rvenriques Z[C,C,B,2];", "quadratic extension of the base field"),
        ("rvhyperelliptic D[W]{zw,ow,zw,zw,zw,ow}; rvenriques Z[C,D,B,2];", "not the declared base extension"),
        ("rvenriques Y[C,L,B,2];", "quotient 'Y' is already defined"),
        ("card (\\y:RVEnriques(Missing)=>S.hi);", "unknown Enriques quotient 'Missing'"),
        ("show zero[RVEnriques(Y)];", "found RVEnriques(Y)"),
        ("show rvEnriquesSplit(S.hi);", "an RVEnriques point"),
        ("rvjacbase Q[C,W]; let p=rvKummerExceptional(zero[RVKummer(Q)],zw,zw,ow,zw); show rvEnriquesOrbit[Y](p);", "requires RVKummerK3(C) or RVKummerK3(L)"),
        ("rvenriques Z[C,L,B,2]; let x=rvKummerExceptional(zero[RVKummer(C)],z,z,o,z); assert rvEnriquesOrbit[Y](x)==rvEnriquesOrbit[Z](x);", "expected type RVEnriques(Y)"),
    ] {
        let error = diagnostic(interpret(&format!("{}\n{statement}", preamble())));
        assert!(error.contains(expected), "{error}");
    }
}
