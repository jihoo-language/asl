mod common;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,o,z,o,t,o};
    rvjacbase L[C,W];
";

#[test]
fn every_base_point_extends_and_lifts_and_projection_commutes_with_frobenius() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            let extend=\p:RVKummer(C)=>rvKummerExtend[L](p);
            assert (\p:RVKummer(C)=>rvKummerDescend[C](extend(p)))==(\p:RVKummer(C)=>p);
            assert (\p:RVKummer(C)=>rvJacKummerPoint(rvKummerLift(extend(p))))==extend;
            assert (\d:RVJac(C)=>extend(rvJacKummerPoint(d)))
                ==(\d:RVJac(C)=>rvJacKummerPoint(rvJacExtend[L](d)));
            assert (\d:RVJac(L)=>rvKummerFrob[C,1](rvJacKummerPoint(d)))
                ==(\d:RVJac(L)=>rvJacKummerPoint(rvJacFrob[C,1](d)));
            assert (\p:RVKummer(L)=>rvKummerFrob[C,2](p))==(\p:RVKummer(L)=>p);
            "
        ))),
        ""
    );
}

#[test]
fn normalized_nonprime_coefficients_use_relative_frobenius() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram U[3,1]{1,0,1};
            rvunram W[3,1]{2,1,0,0,1};
            let z=rvUnramElem[U]{};
            let o=rvUnramElem[U]{1};
            let a=rvUnramElem[U]{0,1};
            rvhyperelliptic C[U]{rvUnramElem[U]{2,1},o,z,z,z,o};
            rvjacbase L[C,W];
            let point=rvKummer[C](z,o,a,neg(o));
            let extended=rvKummerExtend[L](point);
            assert rvKummerCoord[2](extended)==rvUnramElem[W]{0,2,1,1};
            assert rvKummerFrob[C,1](extended)==extended;
            assert rvKummerDescend[C](extended)==point;
            "
        )),
        ""
    );
}

#[test]
fn descended_richelot_map_matches_extension_lifts_on_the_entire_quartic() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-base-descent")),
        concat!(
            "rvRichelotDescent[D] = relative degree=2, permutation=[1, 0, 2], rational kernel order=2\n",
            "RVKummer(D)[[1], [0], [2], [2]]\n",
        )
    );
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            rvrichelot DL[L]{{[0,1],[1],[1]}}{{[0,2],[1],[1]}}{{[],[1]}};
            rvrichelotdescend D[C,DL];
            let push=\p:RVKummer(C)=>rvKummerRichelotDesc[D](p);
            assert (\p:RVKummer(C)=>rvKummerExtend[DL](push(p)))
                ==(\p:RVKummer(C)=>rvKummerRichelot[DL](rvKummerExtend[L](p)));
            let lift=\p:RVKummer(C)=>rvKummerLift(rvKummerExtend[L](p));
            assert push==(\p:RVKummer(C)=>rvKummerDescend[D](rvJacKummerPoint(rvRichelotPush[DL](lift(p)))));
            assert (\d:RVJac(C)=>push(rvJacKummerPoint(d)))
                ==(\d:RVJac(C)=>rvJacKummerPoint(rvRichelotDescPush[D](d)));
            "
        ))),
        ""
    );
}

#[test]
fn descent_rejects_nonfixed_points_and_distinguishes_a_twist_only_lift() {
    let setup = format!(
        r"{F3}
        let zw=rvUnramElem[W]{{}};
        let ow=rvUnramElem[W]{{1}};
        let aw=rvUnramElem[W]{{0,1}};
        let k=rvKummer[L](zw,ow,aw,neg(ow));
        assert rvKummerFrob[C,1](k)==rvKummer[L](zw,ow,neg(aw),neg(ow));
        let twist=rvKummerExtend[L](rvKummer[C](z,o,o,o));
        let lift=rvKummerLift(twist);
        assert rvJacFrob[C,1](lift)==neg(lift);
        assert rvKummerFrob[C,1](twist)==twist;
        "
    );
    assert_eq!(stdout(interpret(&setup)), "");
    for (code, expected) in [
        (
            "show rvKummerDescend[C](k);",
            "Kummer coordinate lies outside the requested subfield",
        ),
        (
            "show rvJacDescend[C](lift);",
            "Mumford coefficient lies outside the requested subfield",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{setup}{code}")));
        assert!(problem.contains(expected), "{problem}");
    }
}

#[test]
fn names_and_curves_are_checked_for_each_transport() {
    let setup = format!(
        r"{F3}
        rvhyperelliptic Different[F]{{z,t,t,z,t,o}};
        rvhyperelliptic Twin[F]{{z,o,z,o,t,o}};
        rvrichelot DL[L]{{[0,1],[1],[1]}}{{[0,2],[1],[1]}}{{[],[1]}};
        rvrichelotdescend D[C,DL];
        assert rvKummerExtend[Twin](zero[RVKummer(C)])==zero[RVKummer(Twin)];
        "
    );
    for (code, expected) in [
        (
            "show rvKummerExtend[L](zero[RVKummer(Different)]);",
            "not the declared base extension",
        ),
        (
            "show rvKummerDescend[Different](zero[RVKummer(L)]);",
            "not the declared base extension",
        ),
        (
            "show rvKummerFrob[Different,0](zero[RVKummer(L)]);",
            "not the declared base extension",
        ),
        (
            "show rvKummerExtend[C](zero[RVKummer(L)]);",
            "requires an extension of the coefficient field",
        ),
        (
            "show rvKummerExtend[L](zero[RVJac(C)]);",
            "requires an RVKummer point",
        ),
        (
            "show rvKummerRichelotDesc[D](zero[RVKummer(Twin)]);",
            "expected type RVKummer(C)",
        ),
        (
            "show rvKummerRichelotDesc[DL](zero[RVKummer(L)]);",
            "has no stored Richelot descent",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{setup}{code}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
