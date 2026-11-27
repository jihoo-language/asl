mod common;

use common::{diagnostic, example, interpret, stdout};

const F11: &str = "rvunram F11[11,1]{0,1}; rvec E[F11](rvUnramElem[F11]{1},rvUnramElem[F11]{7}); rvellipticcm K[E];";

#[test]
fn recognizes_lang_scalings_and_the_same_map_over_a_larger_base_field() {
    let output = stdout(example("elliptic-cm-recognition"));
    assert_eq!(
        output.lines().collect::<Vec<_>>(),
        [
            "rvECCMCertify[K] = degree=15, trace=5, trace_product=10, [D]alpha=[A]+[B]pi with (A,B,D)=(1,-1,1)",
            "rvECCMCertify[K2] = degree=15, trace=5, trace_product=-85, [D]alpha=[A]+[B]pi with (A,B,D)=(14,1,3)",
            "rvECCMCertify[L] = degree=1, trace=-1, trace_product=-2, [D]alpha=[A]+[B]pi with (A,B,D)=(1,1,3)",
            "rvEllLocalOrder[C, 3] = r=1, bound=1, m=3, a=2, trace=-3, norm=3, discriminant=-3",
        ]
    );
}

#[test]
fn handles_zero_identity_doubling_opposites_and_inseparable_frobenius() {
    let output = stdout(interpret(&format!(
        r"
        {F11}
        let id = rvECChart[E,E](rvUnramElem[F11]{{1}});
        let zero_map = rvECSum(id,rvECNegate(id));
        let doubled = rvECSum(id,id);
        let pi = rvECFrobenius[E];
        assert rvECSum(zero_map,pi) == pi;
        assert rvECSum(pi,zero_map) == pi;
        assert rvECCMRecognize[K](zero_map) == rvQpQuad[K,0,0,1];
        assert rvECCMRecognize[K](doubled) == rvQpQuad[K,2,0,1];
        assert rvECCMRecognize[K](rvECNegate(id)) == rvQpQuad[K,-1,0,1];
        rvECCMCertify[K](pi);
    "
    )));
    assert!(
        output.contains("degree=11, trace=-3, trace_product=-13"),
        "{output}"
    );
}

#[test]
fn trace_polarization_matches_direct_composition_and_the_cm_quadratic_identity() {
    let output = stdout(interpret(
        r"
        rvunram F13[13,1]{0,1};
        rvec E[F13](rvUnramElem[F13]{0},rvUnramElem[F13]{2});
        rvellipticcm K[E];
        let alpha = rvECChart[E,E](rvUnramElem[F13]{3});
        let pi = rvECFrobenius[E];
        let ap = rvECMapCompose(alpha,pi);
        rvECCMCertify[K](alpha);
        rvECCMCertify[K](ap);
        let sq = rvECMapCompose(alpha,alpha);
        assert rvECSum(sq,rvECId[E]) == rvECNegate(alpha);
        assert rvECMapCompose(pi,alpha) == ap;
        assert rvECMapCompose(rvECSum(alpha,rvECId[E]),pi) == rvECSum(ap,pi);
    ",
    ));
    assert!(
        output.contains("degree=1, trace=-1, trace_product=-2"),
        "{output}"
    );
    assert!(
        output.contains("degree=13, trace=-2, trace_product=23"),
        "{output}"
    );
}

#[test]
fn checks_the_actual_embedding_of_nonconstant_curve_coefficients() {
    stdout(interpret(
        r"
        rvunram F25[5,1]{2,0,1};
        rvunram F625[5,1]{2,0,0,0,1};
        rvec E[F25](rvUnramElem[F25]{0,1},rvUnramElem[F25]{0,0});
        rvec E4[F625](rvUnramElem[F625]{0,0,1,0},rvUnramElem[F625]{0,0,0,0});
        let one = rvECChart[E,E](rvUnramElem[F25]{1,0});
        assert rvECExtend[E4,E4](one) == rvECChart[E4,E4](rvUnramElem[F625]{1,0,0,0});
    ",
    ));
}

#[test]
fn retains_nominal_endpoint_type_checks_and_rejects_incompatible_scalars_and_embeddings() {
    for (source, expected) in [
        (
            "show rvECChart[E,E](rvUnramElem[F11]{0});",
            "scale must be nonzero",
        ),
        (
            "show rvECChart[E,E](rvUnramElem[F11]{2});",
            "does not give the target equation",
        ),
        (
            "rvec E1[F11](rvUnramElem[F11]{1},rvUnramElem[F11]{7}); show rvECSum(rvECLang[E],rvECLang[E1]);",
            "rvECSum",
        ),
        (
            "rvec E1[F11](rvUnramElem[F11]{1},rvUnramElem[F11]{7}); show rvECCMRecognize[K](rvECChart[E,E1](rvUnramElem[F11]{1}));",
            "a rational endomorphism",
        ),
        (
            "rvec E1[F11](rvUnramElem[F11]{0},rvUnramElem[F11]{2}); show rvECExtend[E1,E](rvECLang[E]);",
            "does not give the declared endpoint",
        ),
        (
            "rvunram W[13,1]{0,1}; rvec E1[W](rvUnramElem[W]{0},rvUnramElem[W]{2}); show rvECCMRecognize[K](rvECLang[E1]);",
            "CM field does not match",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F11} {source}")));
        assert!(error.contains(expected), "{source}\n{error}");
    }
}
