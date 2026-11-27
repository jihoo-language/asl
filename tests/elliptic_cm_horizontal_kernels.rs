mod common;

use common::{diagnostic, example, interpret, stdout};

const D48: &str = "rvcmclasspoly H3[-3]; rvcmclasspolylift H12[H3,2];
    rvcmclasspolylift H48[H12,2]; rvunram F[13,1]{0,1};
    rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{2}); rvecendorder M[E,H48];";

#[test]
fn degree_seven_distinct_kernels_have_the_same_other_cm_j() {
    assert_eq!(
        stdout(example("elliptic-cm-horizontal-kernels")),
        concat!(
            "rvECCMHorizontalData[M,7,0] = curve=E, eigenvalue=5, other_eigenvalue=4, x_candidate_degree=3, kernel=[[8], [7], [8], [1]]\n",
            "rvECCMHorizontalData[M,7,1] = curve=E, eigenvalue=4, other_eigenvalue=5, x_candidate_degree=3, kernel=[[12], [6], [12], [1]]\n",
            "rvECEndData[N] = curve=Q, q=13, trace=2, frobenius_discriminant=-48, discriminant=-48, j=[4], index=1, frobenius=[1, 1]\n",
        )
    );
}

#[test]
fn opposite_eigenvalues_share_the_x_equation_but_have_different_y_signs() {
    let output = stdout(interpret(
        "rvcmclasspoly H8[-8]; rvunram F[11,1]{0,1};
         rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{8});
         rvecendorder M[E,H8];
         let negative=rvECCMKernel[M,3,0];
         let positive=rvECCMKernel[M,3,1];
         assert negative==rveckernel[E,3]{[1],[1]};
         assert positive==rveckernel[E,3]{[2],[1]};
         let point=rvECPoint[E](rvUnramElem[F]{9},rvUnramElem[F]{3});
         assert rvECScale[3](point)==rvECInfinity[E];
         rvECCMHorizontalData[M,3,0];
         rvECCMHorizontalData[M,3,1];
         rveccmhorizontal I[M,Q,3,0];
         rvecisogeny J[E,R](positive);
         assert rvECA[Q]==rvUnramElem[F]{5};
         assert rvECB[Q]==rvUnramElem[F]{6};
         assert rvECA[R]==rvUnramElem[F]{3};
         assert rvECB[R]==rvUnramElem[F]{10};
         assert rvECEndJ[Q]==rvUnramElem[F]{3};
         assert rvECEndJ[R]==rvUnramElem[F]{3};
         assert rvECApply(rvECIso[J],point)==rvECInfinity[R];
         rvunram L[11,1]{5,0,1};
         rvec V[L](rvUnramElem[L]{1},rvUnramElem[L]{8});
         let y=rvUnramElem[L]{0,1};
         let p=rvECPoint[V](rvUnramElem[L]{10},y);
         assert rvECScale[3](p)==rvECInfinity[V];
         assert rvUnramFrob[1](y)==neg(y);
         rvecendorder N[R,H8];
         rvECEndData[N];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvECCMHorizontalData[M,3,0] = curve=E, eigenvalue=2, other_eigenvalue=1, x_candidate_degree=2, kernel=[[1], [1]]\n",
            "rvECCMHorizontalData[M,3,1] = curve=E, eigenvalue=1, other_eigenvalue=2, x_candidate_degree=2, kernel=[[2], [1]]\n",
            "rvECEndData[N] = curve=R, q=11, trace=6, frobenius_discriminant=-8, discriminant=-8, j=[3], index=1, frobenius=[3, 1]\n",
        )
    );
}

#[test]
fn order_three_artin_class_cycles_the_three_different_cm_targets() {
    let output = stdout(interpret(
        "rvcmclasspoly H3[-3]; rvcmclasspolylift H12[H3,2];
         rvcmclasspolylift H108[H12,3]; rvcmorder O[-108];
         rvunram F[31,1]{0,1}; rvec E[F](rvUnramElem[F]{3},rvUnramElem[F]{15});
         rvecendorder M[E,H108];
         rvCMClassOrder(rvCMArtinClass[O,7,1]);
         rveccmhorizontal I[M,First,7,1];
         rveccmhorizontal J[M,Inverse,7,0];
         assert rvECEndJ[E]==rvUnramElem[F]{18};
         assert rvECEndJ[First]==rvUnramElem[F]{21};
         assert rvECEndJ[Inverse]==rvUnramElem[F]{26};
         rvecendorder M1[First,H108]; rveccmhorizontal K[M1,Second,7,1];
         assert rvECEndJ[Second]==rvUnramElem[F]{26};
         rvecendorder M2[Second,H108]; rveccmhorizontal T[M2,Third,7,1];
         assert rvECEndJ[Third]==rvECEndJ[E];",
    ));
    assert_eq!(output, "rvCMClassOrder = 3\n");
}

#[test]
fn rejects_unsplit_bad_degree_and_repeated_eigenvalues() {
    for (query, expected) in [
        ("show rvECCMKernel[M,2,0];", "degree must be an odd prime"),
        ("show rvECCMKernel[M,9,0];", "9 is not prime"),
        (
            "show rvECCMKernel[M,13,0];",
            "equals the residue characteristic",
        ),
        (
            "show rvECCMKernel[M,3,0];",
            "must not divide the order discriminant",
        ),
        ("show rvECCMKernel[M,5,0];", "prime must split"),
        ("show rvECCMKernel[M,7,2];", "orientation must be 0 or 1"),
        (
            "rveccmhorizontal I[M,E,7,0];",
            "target name is already defined",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{D48} {query}")));
        assert!(error.contains(expected), "{error}");
    }
    let error = diagnostic(interpret(
        "rvcmclasspoly H8[-8]; rvunram L[11,1]{1,0,1};
         rvec E[L](rvUnramElem[L]{1},rvUnramElem[L]{8}); rvecendorder M[E,H8];
         show rvECCMKernel[M,3,0];",
    ));
    assert!(
        error.contains("divides the Frobenius order index and has repeated eigenvalues"),
        "{error}"
    );
}
