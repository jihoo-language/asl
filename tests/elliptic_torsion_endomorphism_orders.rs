mod common;

use common::{diagnostic, example, interpret, stdout};

const F13: &str = "rvunram F[13,1]{0,1};
    rvec E[F](rvUnramElem[F]{4},rvUnramElem[F]{0});
    rvec V[F](rvUnramElem[F]{2},rvUnramElem[F]{1});
    rvunram L[13,1]{11,0,1};
    rvec VL[L](rvUnramElem[L]{2},rvUnramElem[L]{1});
    let p=rvECPoint[E](rvUnramElem[F]{0},rvUnramElem[F]{0});
    let q=rvECPoint[E](rvUnramElem[F]{3},rvUnramElem[F]{0});
    let r=rvECPoint[VL](rvUnramElem[L]{2},rvUnramElem[L]{0});
    let s=rvECPoint[VL](rvUnramElem[L]{12,2},rvUnramElem[L]{0});
    rvecplainlevel A[E,E,2](p,q); rvecplainlevel B[V,VL,2](r,s);";

#[test]
fn same_trace_orders_and_composite_evidence_without_class_polynomials() {
    let output = stdout(example("elliptic-torsion-endomorphism-orders"));
    assert_eq!(
        output,
        concat!(
            "rvECEndTorsionCheck[E] = q=13, trace=6, frobenius_discriminant=-16, fundamental_discriminant=-4, frobenius_conductor=2, index_bounds=[2, 2], conductor_bounds=[1, 1], discriminant=-4, depths=[(2, 1, 1)], missing_levels=[]\n",
            "rvECEndTorsionCheck[V] = q=13, trace=6, frobenius_discriminant=-16, fundamental_discriminant=-4, frobenius_conductor=2, index_bounds=[1, 1], conductor_bounds=[2, 2], discriminant=-16, depths=[(2, 0, 0)], missing_levels=[]\n",
            "rvECEndData[M] = curve=E, q=13, trace=6, frobenius_discriminant=-16, discriminant=-4, j=[12], index=2, frobenius=[3, 2]\n",
            "rvECEndData[N] = curve=V, q=13, trace=6, frobenius_discriminant=-16, discriminant=-16, j=[1], index=1, frobenius=[3, 1]\n",
            "rvECEndTorsionCheck[U] = q=37, trace=2, frobenius_discriminant=-144, fundamental_discriminant=-4, frobenius_conductor=6, index_bounds=[2, 6], conductor_bounds=[1, 3], discriminant=none, depths=[(2, 1, 1), (3, 0, 1)], missing_levels=[3]\n",
            "rvECEndTorsionCheck[U] = q=37, trace=2, frobenius_discriminant=-144, fundamental_discriminant=-4, frobenius_conductor=6, index_bounds=[3, 6], conductor_bounds=[1, 2], discriminant=none, depths=[(2, 0, 1), (3, 1, 1)], missing_levels=[2]\n",
            "rvECEndTorsionCheck[U] = q=37, trace=2, frobenius_discriminant=-144, fundamental_discriminant=-4, frobenius_conductor=6, index_bounds=[6, 6], conductor_bounds=[1, 1], discriminant=-4, depths=[(2, 1, 1), (3, 1, 1)], missing_levels=[]\n",
            "rvECPlainLevelScalar[U6,2] = available_depth=1, scalar_depth=1, scalar=1, scalar_modulus=2, divided_level=3, divided_columns=[[0, 0], [0, 0]], divided_trace=0, divided_norm=9\n",
            "rvECEndData[R] = curve=U, q=37, trace=2, frobenius_discriminant=-144, discriminant=-4, j=[26], index=6, frobenius=[1, 6]\n",
        )
    );
}

#[test]
fn scalar_prefix_requires_more_precision_and_next_nonscalar_depth_is_exact() {
    let output = stdout(interpret(
        "rvunram F[17,1]{0,1}; rvunram L[17,1]{14,0,1};
         rvec E[F](rvUnramElem[F]{6},rvUnramElem[F]{3});
         rvec EL[L](rvUnramElem[L]{6},rvUnramElem[L]{3});
         rvecplainlevel B2[E,E,2](rvECPoint[E](rvUnramElem[F]{6},rvUnramElem[F]{0}),
                                  rvECPoint[E](rvUnramElem[F]{13},rvUnramElem[F]{0}));
         rvECEndTorsionCheck[E,B2];
         let p=rvECPoint[EL](rvUnramElem[L]{6,2},rvUnramElem[L]{1,7});
         let q=rvECPoint[EL](rvUnramElem[L]{13,1},rvUnramElem[L]{2,10});
         rvecplainlevel B4[E,EL,4](p,q);
         assert rvECLevelFrob[B4]==rvECLevelMat[B4,3,2,2,3];
         rvECPlainLevelScalar[B4,2]; rvECEndTorsionCheck[E,B2,B4,B2];
         rvecendtorsion M[E,B4,B2]; rvECEndData[M];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvECEndTorsionCheck[E] = q=17, trace=2, frobenius_discriminant=-64, fundamental_discriminant=-4, frobenius_conductor=4, index_bounds=[2, 4], conductor_bounds=[1, 2], discriminant=none, depths=[(2, 1, 2)], missing_levels=[4]\n",
            "rvECPlainLevelScalar[B4,2] = available_depth=2, scalar_depth=1, scalar=1, scalar_modulus=2, divided_level=2, divided_columns=[[1, 1], [1, 1]], divided_trace=0, divided_norm=4\n",
            "rvECEndTorsionCheck[E] = q=17, trace=2, frobenius_discriminant=-64, fundamental_discriminant=-4, frobenius_conductor=4, index_bounds=[2, 2], conductor_bounds=[2, 2], discriminant=-16, depths=[(2, 1, 1)], missing_levels=[]\n",
            "rvECEndData[M] = curve=E, q=17, trace=2, frobenius_discriminant=-64, discriminant=-16, j=[9], index=2, frobenius=[1, 2]\n",
        )
    );
}

#[test]
fn scalar_through_the_conductor_cap_needs_no_next_level() {
    let output = stdout(interpret(
        "rvunram F[17,1]{0,1};
         rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{0});
         rvecplainlevel B[E,E,4](rvECPoint[E](rvUnramElem[F]{1},rvUnramElem[F]{6}),
                                 rvECPoint[E](rvUnramElem[F]{3},rvUnramElem[F]{8}));
         assert rvECLevelFrob[B]==rvECLevelMat[B,1,0,0,1];
         rvECPlainLevelScalar[B,2]; rvECEndTorsionCheck[E,B];
         rvecendtorsion M[E,B];
         assert rvECEndFrobenius[M]==rvECEndElem[M,1,4];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvECPlainLevelScalar[B,2] = available_depth=2, scalar_depth=2, scalar=1, scalar_modulus=4\n",
            "rvECEndTorsionCheck[E] = q=17, trace=2, frobenius_discriminant=-64, fundamental_discriminant=-4, frobenius_conductor=4, index_bounds=[4, 4], conductor_bounds=[1, 1], discriminant=-4, depths=[(2, 2, 2)], missing_levels=[]\n",
        )
    );
}

#[test]
fn independent_certificates_agree_with_class_polynomials_and_order_arithmetic() {
    let output = stdout(interpret(&format!(
        "{F13}
         rvECLevelInvariants(rvECLevelFrob[A]); rvECLevelInvariants(rvECLevelFrob[B]);
         rvecendtorsion M[E,A]; rvecendtorsion N[V,B];
         assert rvECEndFrobenius[M]==rvECEndElem[M,3,2];
         assert rvECEndFrobenius[N]==rvECEndElem[N,3,1];
         assert rvECEndMul(rvECEndElem[N,0,1],rvECEndElem[N,0,1])==rvECEndElem[N,-4,0];
         rvcmclasspoly H4[-4]; rvcmclasspolylift H16[H4,2];
         rvecendorder D[E,H4]; rvecendorder T[V,H16];
         rvECEndData[M]; rvECEndData[D]; rvECEndData[N]; rvECEndData[T];"
    )));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        lines[0],
        "rvECLevelInvariants = modulus=2, trace=0, determinant=1"
    );
    assert_eq!(lines[0], lines[1]);
    assert_eq!(
        lines[2].split_once(" = ").unwrap().1,
        lines[3].split_once(" = ").unwrap().1
    );
    assert_eq!(
        lines[4].split_once(" = ").unwrap().1,
        lines[5].split_once(" = ").unwrap().1
    );
}

#[test]
fn odd_fundamental_and_frobenius_already_maximal_orders() {
    let output = stdout(interpret(
        "rvunram F[7,1]{0,1}; rvec E[F](rvUnramElem[F]{0},rvUnramElem[F]{2});
         rvecplainlevel B[E,E,3](rvECPoint[E](rvUnramElem[F]{0},rvUnramElem[F]{3}),
                                 rvECPoint[E](rvUnramElem[F]{3},rvUnramElem[F]{1}));
         rvecendtorsion M[E,B]; rvECEndData[M];
         rvunram G[11,1]{0,1}; rvec X[G](rvUnramElem[G]{1},rvUnramElem[G]{8});
         rvecendtorsion N[X]; rvECEndData[N];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvECEndData[M] = curve=E, q=7, trace=-1, frobenius_discriminant=-27, discriminant=-3, j=[0], index=3, frobenius=[-2, 3]\n",
            "rvECEndData[N] = curve=X, q=11, trace=6, frobenius_discriminant=-8, discriminant=-8, j=[3], index=1, frobenius=[3, 1]\n",
        )
    );
}

#[test]
fn rejects_unwitnessed_nominal_or_nonordinary_data() {
    for (query, expected) in [
        ("rvecendtorsion M[E];", "insufficient Frobenius scalar-depth evidence"),
        ("rvECEndTorsionCheck[E,B];", "same named base curve"),
        ("rvecplainlevel C[E,VL,2](r,s);", "not the declared torsion base extension"),
        ("rvecplainlevel C[V,VL,2](p,q);", "elliptic level generator"),
        ("rvecplainlevel C[E,E,2](p,p);", "dependent modulo the level"),
        ("rvECPlainLevelScalar[A,3];", "prime must divide the named level"),
        ("rvecendtorsion M[E,A]; rveccmaction D[M]; rveclevel C[D,2](p,q); rvECEndTorsionCheck[E,C];", "requires a plain level"),
        ("rvunram U[7,1]{0,1}; rvec X[U](rvUnramElem[U]{1},rvUnramElem[U]{0}); rvecendtorsion M[X];", "requires an ordinary elliptic curve"),
    ] {
        let error = diagnostic(interpret(&format!("{F13} {query}")));
        assert!(error.contains(expected), "{error}");
    }
}
