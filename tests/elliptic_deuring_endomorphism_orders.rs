mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn same_frobenius_polynomial_different_endomorphism_orders() {
    assert_eq!(
        stdout(example("elliptic-deuring-endomorphism-orders")),
        concat!(
            "rvECEndCheck[G,H16] = p=13, q=13, trace=6, frobenius_discriminant=-16, candidate_discriminant=-16, j=[12], residual=[11], ordinary=true, split=true, index=1, certified=false\n",
            "rvECEndCheck[N,H4] = p=13, q=13, trace=6, frobenius_discriminant=-16, candidate_discriminant=-4, j=[1], residual=[2], ordinary=true, split=true, index=2, certified=false\n",
            "rvECEndData[Maximal] = curve=G, q=13, trace=6, frobenius_discriminant=-16, discriminant=-4, j=[12], index=2, frobenius=[3, 2]\n",
            "rvECEndData[Nonmaximal] = curve=N, q=13, trace=6, frobenius_discriminant=-16, discriminant=-16, j=[1], index=1, frobenius=[3, 1]\n",
            "rvECEndTrace = 6\nrvECEndNorm = 13\n",
            "rvECEndCheck[Sup,H3] = p=5, q=5, trace=0, frobenius_discriminant=-20, candidate_discriminant=-3, j=[0], residual=[0], ordinary=false, split=false, index=none, certified=false\n",
        )
    );
}

#[test]
fn twists_and_constant_field_extension_preserve_the_end_order() {
    let output = stdout(interpret(
        "rvcmclasspoly H4[-4]; rvcmclasspolylift H16[H4,2];
         rvunram F[13,1]{0,1}; rvunram L[13,1]{11,0,1};
         rvec T[F](rvUnramElem[F]{3},rvUnramElem[F]{});
         rvec V[F](rvUnramElem[F]{8},rvUnramElem[F]{8});
         rvec G[L](rvUnramElem[L]{4},rvUnramElem[L]{});
         rvec N[L](rvUnramElem[L]{2},rvUnramElem[L]{1});
         rvecendorder R[T,H4]; rvecendorder W[V,H16];
         rvecendorder I[G,H4]; rvecendorder J[N,H16];
         assert rvECEndFrobenius[R]==rvECEndElem[R,-3,2];
         assert rvECEndFrobenius[W]==rvECEndElem[W,-3,1];
         assert rvECEndMul(rvECEndElem[I,3,2],rvECEndElem[I,3,2])==rvECEndFrobenius[I];
         assert rvECEndFrobenius[J]==rvECEndElem[J,5,6];
         assert rvECEndJ[G]==rvUnramElem[L]{12};
         rvECEndData[R]; rvECEndData[W]; rvECEndData[I]; rvECEndData[J];",
    ));
    assert_eq!(output, concat!(
        "rvECEndData[R] = curve=T, q=13, trace=-6, frobenius_discriminant=-16, discriminant=-4, j=[12], index=2, frobenius=[-3, 2]\n",
        "rvECEndData[W] = curve=V, q=13, trace=-6, frobenius_discriminant=-16, discriminant=-16, j=[1], index=1, frobenius=[-3, 1]\n",
        "rvECEndData[I] = curve=G, q=169, trace=10, frobenius_discriminant=-576, discriminant=-4, j=[12, 0], index=12, frobenius=[5, 12]\n",
        "rvECEndData[J] = curve=N, q=169, trace=10, frobenius_discriminant=-576, discriminant=-16, j=[1, 0], index=6, frobenius=[5, 6]\n",
    ));
}

#[test]
fn odd_discriminant_arithmetic_and_wide_norm() {
    let max = i128::MAX;
    let output = stdout(interpret(&format!(
        "rvcmclasspoly H[-3]; rvunram F[13,1]{{0,1}};
         rvec E[F](rvUnramElem[F]{{}},rvUnramElem[F]{{1}});
         rvecendorder O[E,H];
         let w=rvECEndElem[O,0,1];
         let x=rvECEndElem[O,7,2]; let y=rvECEndElem[O,-3,4];
         assert rvECEndMul(w,w)==add(w,rvECEndElem[O,-1,0]);
         assert rvECEndConj(w)==rvECEndElem[O,1,-1];
         assert rvECEndConj(rvECEndConj(x))==x;
         assert rvECEndMul(rvECEndMul(x,y),w)==rvECEndMul(x,rvECEndMul(y,w));
         assert rvECEndMul(x,rvECEndConj(x))==rvECEndElem[O,67,0];
         assert add(neg(x),x)==zero[RVECEnd(O)];
         assert rvECEndFrobenius[O]==rvECEndElem[O,-1,4];
         rvECEndTrace(x); rvECEndNorm(x); rvECEndNorm(rvECEndElem[O,{max},0]);"
    )));
    assert_eq!(output, concat!(
        "rvECEndTrace = 16\nrvECEndNorm = 67\n",
        "rvECEndNorm = 28948022309329048855892746252171976962977213799489202546401021394546514198529\n",
    ));
}

#[test]
fn genuinely_nonconstant_j_in_a_quadratic_field() {
    let output = stdout(interpret(
        "rvcmclasspoly H3[-3]; rvcmclasspolylift H12[H3,2]; rvcmclasspolylift H48[H12,2];
         rvunram L[7,1]{5,5,1};
         rvec E[L](rvUnramElem[L]{6,4},rvUnramElem[L]{4,5});
         let z=rvUnramElem[L]{0,1};
         assert rvECEndJ[E]==z;
         assert rvCMClassPolyEval[H48](z)==zero[RVUnram(L)];
         rvecendorder O[E,H48]; rvECEndData[O];",
    ));
    assert_eq!(
        output,
        "rvECEndData[O] = curve=E, q=49, trace=-2, frobenius_discriminant=-192, discriminant=-48, j=[0, 1], index=2, frobenius=[-1, 2]\n"
    );
}

#[test]
fn ordinary_root_in_the_characteristic_conductor_is_rejected() {
    let setup = "rvcmclasspoly H3[-3]; rvcmclasspolylift H147[H3,7];
         rvunram F[7,1]{0,1}; rvec E[F](rvUnramElem[F]{},rvUnramElem[F]{1});";
    let output = stdout(interpret(&format!("{setup} rvECEndCheck[E,H147];")));
    assert_eq!(output,
        "rvECEndCheck[E,H147] = p=7, q=7, trace=-4, frobenius_discriminant=-12, candidate_discriminant=-147, j=[0], residual=[0], ordinary=true, split=false, index=none, certified=false\n",
    );
    let error = diagnostic(interpret(&format!("{setup} rvecendorder Bad[E,H147];")));
    assert!(error.contains("must be a positive square"), "{error}");
}

#[test]
fn invalid_certificates_and_nominal_contexts_are_diagnosed() {
    let setup = "rvcmclasspoly H4[-4]; rvcmclasspolylift H16[H4,2]; rvcmclasspoly H3[-3];
        rvunram F[13,1]{0,1}; rvunram Twin[13,1]{0,1}; rvunram G[5,1]{0,1};
        rvec E[F](rvUnramElem[F]{4},rvUnramElem[F]{});
        rvec N[F](rvUnramElem[F]{2},rvUnramElem[F]{1});
        rvec Sup[G](rvUnramElem[G]{},rvUnramElem[G]{1});
        rvecendorder R[E,H4]; rvecendorder U[N,H16];";
    for (command, expected) in [
        ("rvecendorder Bad[Sup,H3];", "requires an ordinary curve"),
        ("rvecendorder Bad[E,H16];", "j-invariant is not a root"),
        ("rvecendorder R[E,H4];", "certificate is already defined"),
        ("show rvECEndElem[Missing,1,0];", "unknown elliptic endomorphism-order certificate"),
        ("show rvECEndMul(rvECEndFrobenius[R],rvECEndFrobenius[U]);", "expected type RVECEnd(R)"),
        ("assert rvECEndJ[E]==rvUnramElem[Twin]{12};", "expected type RVUnram(F)"),
        ("show rvECEndMul(rvECEndElem[R,170141183460469231731687303715884105727,0],rvECEndElem[R,2,0]);", "arithmetic exceeds signed i128 range"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}
