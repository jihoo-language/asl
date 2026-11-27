mod common;

use common::{diagnostic, example, interpret, stdout};

const D8: &str = "rvcmclasspoly H8[-8]; rvunram F[11,1]{0,1};
    rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{8}); rvecendorder M[E,H8];";

#[test]
fn oriented_order_three_cycles_certify_the_chart_degree_and_point_action() {
    assert_eq!(
        stdout(example("elliptic-cm-principal-isogeny-cycles")),
        concat!(
            "rvECCMCycleData[K] = curve=E, ell=7, orientation=1, length=3, degree=343, trace=20, element=[10, 3], chart_scale=[8], chart_count=2, differential=[4], j_path=[[18], [21], [26], [18]]\n",
            "rvECCMCycleData[L] = curve=E, ell=7, orientation=0, length=3, degree=343, trace=20, element=[10, -3], chart_scale=[2], chart_count=2, differential=[16], j_path=[[18], [26], [21], [18]]\n",
            "rvECCMCycleStep[K,0] = source_a=[3], source_b=[15], target_a=[20], target_b=[9], kernel=[[28], [28], [15], [1]]\n",
            "rvECCMCycleStep[K,1] = source_a=[20], source_b=[9], target_a=[20], target_b=[6], kernel=[[22], [28], [6], [1]]\n",
            "rvECCMCycleStep[K,2] = source_a=[20], source_b=[6], target_a=[24], target_b=[29], kernel=[[1], [28], [3], [1]]\n",
            "rvECEndNorm = 343\n",
            "rvECCMCycleData[J] = curve=S, ell=3, orientation=1, length=1, degree=3, trace=2, element=[1, -1], chart_scale=[3], chart_count=2, differential=[4], j_path=[[3], [3]]\n",
            "rvECCMCycleVerify[J] = degree=3, rational_map_agrees=true\n",
        )
    );
}

#[test]
fn degree_three_opposite_eigenlines_are_distinguished_in_the_actual_maps() {
    let source = format!(
        "{D8}
         rveccmcycle A[M,3,0]; rveccmcycle B[M,3,1];
         assert rvECCMCycleElement[A]==rvECEndElem[M,1,1];
         assert rvECCMCycleElement[B]==rvECEndConj(rvECCMCycleElement[A]);
         assert rvECCMCycleApply[A](rvECInfinity[E])==rvECInfinity[E];
         let p=rvECPoint[E](rvUnramElem[F]{{9}},rvUnramElem[F]{{3}});
         assert rvECCMCycleApply[A](p)==neg(p);
         assert rvECCMCycleApply[B](p)==rvECInfinity[E];
         rvECCMCycleData[A]; rvECCMCycleVerify[A];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvECCMCycleData[A] = curve=E, ell=3, orientation=0, length=1, degree=3, trace=2, element=[1, 1], chart_scale=[5], chart_count=2, differential=[9], j_path=[[3], [3]]\n",
            "rvECCMCycleVerify[A] = degree=3, rational_map_agrees=true\n",
        )
    );
}

#[test]
fn gaussian_cm_selects_among_four_returning_charts() {
    let output = stdout(interpret(
        "rvcmclasspoly H4[-4]; rvunram F[13,1]{0,1};
         rvec E[F](rvUnramElem[F]{4},rvUnramElem[F]{0}); rvecendorder M[E,H4];
         rveccmcycle A[M,5,0]; rvECCMCycleData[A];
         assert rvECCMCycleElement[A]==rvECEndElem[M,2,1];
         rvECCMCycleVerify[A];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvECCMCycleData[A] = curve=E, ell=5, orientation=0, length=1, degree=5, trace=4, element=[2, 1], chart_scale=[2], chart_count=4, differential=[7], j_path=[[12], [12]]\n",
            "rvECCMCycleVerify[A] = degree=5, rational_map_agrees=true\n",
        )
    );
}

#[test]
fn nonconstant_extension_model_agrees_with_gaussian_cm_point_action() {
    let output = stdout(interpret(
        "rvcmclasspoly H4[-4]; rvunram L[13,1]{2,0,1};
         rvec E[L](rvUnramElem[L]{1,1},rvUnramElem[L]{0}); rvecendorder M[E,H4];
         rveccmcycle K[M,5,0]; rvECCMCycleData[K]; rvECCMCycleStep[K,0];
         let p=rvECPoint[E](rvUnramElem[L]{2,0},rvUnramElem[L]{2,7});
         let iP=rvECPoint[E](rvUnramElem[L]{11,0},rvUnramElem[L]{10,9});
         assert rvECCMCycleElement[K]==rvECEndElem[M,2,1];
         assert rvECCMCycleApply[K](p)==add(rvECScale[2](p),iP);",
    ));
    assert_eq!(
        output,
        concat!(
            "rvECCMCycleData[K] = curve=E, ell=5, orientation=0, length=1, degree=5, trace=4, element=[2, 1], chart_scale=[2, 0], chart_count=4, differential=[7, 0], j_path=[[12, 0], [12, 0]]\n",
            "rvECCMCycleStep[K,0] = source_a=[1, 1], source_b=[0, 0], target_a=[9, 9], target_b=[0, 0], kernel=[[10, 10], [0, 0], [1, 0]]\n",
        )
    );
}

#[test]
fn rejects_incompatible_curve_types_and_bad_cycle_addresses() {
    for (query, expected) in [
        ("rvECCMCycleData[K];", "unknown CM isogeny cycle 'K'"),
        ("rveccmcycle K[M,2,0];", "degree must be an odd prime"),
        (
            "rveccmcycle K[M,11,0];",
            "equals the residue characteristic",
        ),
        ("rveccmcycle K[M,5,0];", "prime must split"),
        ("rveccmcycle K[M,3,3];", "orientation must be 0 or 1"),
        (
            "rveccmcycle K[M,3,0]; rvECCMCycleStep[K,1];",
            "step index is out of range",
        ),
        (
            "rveccmcycle K[M,3,0]; rveccmcycle K[M,3,1];",
            "cycle is already defined",
        ),
        (
            "rvec S[F](rvUnramElem[F]{1},rvUnramElem[F]{8}); rveccmcycle K[M,3,0];
             show rvECCMCycleApply[K](rvECInfinity[S]);",
            "CM isogeny cycle argument",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{D8} {query}")));
        assert!(error.contains(expected), "{error}");
    }
    let error = diagnostic(interpret(
        "rvcmclasspoly H8[-8]; rvunram L[11,1]{1,0,1};
         rvec E[L](rvUnramElem[L]{1},rvUnramElem[L]{8}); rvecendorder M[E,H8];
         rveccmcycle K[M,3,0];",
    ));
    assert!(error.contains("has repeated eigenvalues"), "{error}");
}
