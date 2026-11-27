mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

const BASE: &str = "rvcmclasspoly H8[-8]; rvunram F[17,1]{0,1};
    rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{7}); rvecendorder M[E,H8];
    rveccmcycle K[M,3,1]; rveccmaction Frob[M]; rveccmaction Full[M,K,K];";

#[test]
fn explicit_full_order_actions_and_base_frobenius_match_the_example() {
    assert_eq!(
        stdout(example("elliptic-cm-order-actions")),
        concat!(
            "rvECCMActionData[FrobeniusOnly] = certificate=M, base=E, curve=E, index=2, discriminant=-32, full=false, generator=[0, 2], frobenius_weight=1, cycle_weights=[]\n",
            "rvECCMActionData[Full] = certificate=M, base=E, curve=E, index=1, discriminant=-8, full=true, generator=[0, 1], frobenius_weight=0, cycle_weights=[(\"K\", 1)]\n",
            "rvECMapDegree = 2\n",
            "rvECCMActionData[AL] = certificate=M, base=E, curve=EL, index=1, discriminant=-8, full=true, generator=[0, 1], frobenius_weight=0, cycle_weights=[(\"K\", 1)]\n",
        )
    );
}

#[test]
fn expanded_action_certifies_the_degree_two_map_exactly() {
    let output = stdout(interpret(&format!(
        "{BASE}
         let w=rvECCMActionGenerator[Full];
         let f=rvECCMActionMap[Full](w);
         rvECMapDegree(f); rvECMapCertify(f);
         assert rvECMapCompose(f,f)==rvECMultMap[E,-2];
         assert rvECMapXNum[0](f)==rvUnramElem[F]{{2}};
         assert rvECMapXNum[1](f)==rvUnramElem[F]{{1}};
         assert rvECMapXNum[2](f)==rvUnramElem[F]{{8}};
         assert rvECMapXDen[0](f)==rvUnramElem[F]{{15}};
         assert rvECMapYNum[0](f)==rvUnramElem[F]{{14}};
         assert rvECMapYNum[1](f)==rvUnramElem[F]{{10}};
         assert rvECMapYNum[2](f)==rvUnramElem[F]{{6}};
         assert rvECMapYDen[0](f)==rvUnramElem[F]{{4}};
         assert rvECMapYDen[1](f)==rvUnramElem[F]{{13}};
         assert rvECCMActionMap[Frob](rvECEndFrobenius[M])==rvECFrobenius[E];"
    )));
    assert!(output.starts_with("rvECMapDegree = 2\n"), "{output}");
}

#[test]
fn omega_agrees_with_an_independent_degree_two_isogeny_on_every_base_point() {
    let mut source = format!(
        "{BASE} let w=rvECCMActionGenerator[Full];
         assert rvECCMAct[Full](w,rvECInfinity[E])==rvECInfinity[E];
         assert rvECCMAct[Full](w,rvECPoint[E](rvUnramElem[F]{{2}},rvUnramElem[F]{{0}}))==rvECInfinity[E];"
    );
    // Quotient the point (2,0) by normalized degree-two Vélu and return
    // (A',B')=(4,12) to (1,7) with scale 5, differential 7=d(omega).
    for (x, y, u, v) in [
        (1, 3, 6, 5),
        (1, 14, 6, 12),
        (5, 1, 1, 3),
        (5, 16, 1, 14),
        (6, 5, 6, 12),
        (6, 12, 6, 5),
        (7, 0, 2, 0),
        (8, 0, 2, 0),
        (12, 8, 1, 3),
        (12, 9, 1, 14),
    ] {
        write!(
            source,
            "assert rvECCMAct[Full](w,rvECPoint[E](rvUnramElem[F]{{{x}}},rvUnramElem[F]{{{y}}}))
                    ==rvECPoint[E](rvUnramElem[F]{{{u}}},rvUnramElem[F]{{{v}}});"
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn extension_actions_preserve_addition_multiplication_and_frobenius() {
    assert_eq!(
        stdout(interpret(&format!(
            "{BASE}
             rvunram L[17,1]{{14,0,1}}; rvec EL[L](rvUnramElem[L]{{1}},rvUnramElem[L]{{7}});
             rveccmactionbase AL[Full,EL]; rveccmactionbase OL[Frob,EL];
             let p=rvECPoint[EL](rvUnramElem[L]{{0,1}},rvUnramElem[L]{{2,1}});
             let t=rvECPoint[EL](rvUnramElem[L]{{1}},rvUnramElem[L]{{3}});
             let w=rvECCMActionGenerator[AL]; let pi=rvECEndFrobenius[M];
             let z=rvECEndElem[M,2,-3];
             let fp=rvECPoint[EL](rvUnramElem[L]{{0,16}},rvUnramElem[L]{{2,16}});
             let wp=rvECPoint[EL](rvUnramElem[L]{{13,6}},rvUnramElem[L]{{7,9}});
             assert rvECCMAct[AL](w,add(p,t))==add(wp,rvECCMAct[AL](w,t));
             assert rvECCMAct[AL](add(w,z),p)==add(wp,rvECCMAct[AL](z,p));
             assert rvECCMAct[AL](rvECEndMul(w,z),p)==rvECCMAct[AL](w,rvECCMAct[AL](z,p));
             assert rvECCMAct[AL](w,fp)==rvECCMAct[OL](pi,wp);
             let beta=rvECCMActionGenerator[OL];
             assert rvECCMAct[OL](beta,p)==add(fp,rvECScale[-3](p));
             assert rvECCMAct[AL](beta,p)==rvECCMAct[OL](beta,p);
             assert rvECMapApply(rvECCMActionMap[AL](w),p)==wp;
             assert rvECMapApply(rvECCMActionMap[OL](pi),p)==fp;
             let m=rvECEndElem[M,0,-170141183460469231731687303715884105728];
             assert rvECCMAct[AL](m,p)==rvECScale[-170141183460469231731687303715884105728](wp);"
        ))),
        ""
    );
}

#[test]
fn two_proper_suborders_saturate_by_a_signed_bezout_combination() {
    let output = stdout(interpret(
        "rvcmclasspoly H8[-8]; rvunram F[73,1]{0,1};
         rvec E[F](rvUnramElem[F]{63},rvUnramElem[F]{64}); rvecendorder M[E,H8];
         rveccmcycle K[M,17,0]; rveccmcycle L[M,19,0];
         assert rvECCMCycleElement[K]==rvECEndElem[M,-3,2];
         assert rvECCMCycleElement[L]==rvECEndElem[M,1,-3];
         rveccmaction A[M,K]; rveccmaction B[M,L]; rveccmaction C[M,K,L];
         assert rvECCMActionGenerator[A]==rvECEndElem[M,0,2];
         assert rvECCMActionGenerator[B]==rvECEndElem[M,0,3];
         let w=rvECCMActionGenerator[C];
         assert w==rvECEndElem[M,0,1];
         let p=rvECPoint[E](rvUnramElem[F]{0},rvUnramElem[F]{8});
         let wp=rvECPoint[E](rvUnramElem[F]{39},rvUnramElem[F]{70});
         assert rvECCMAct[C](w,p)==wp;
         rvECCMActionData[C];",
    ));
    // Independent degree-two Vélu oracle: kernel (45,0), quotient
    // (33,72), and return scale 67 for d(omega) = -1/6 in F73.
    assert_eq!(
        output,
        "rvECCMActionData[C] = certificate=M, base=E, curve=E, index=1, discriminant=-8, full=true, generator=[0, 1], frobenius_weight=0, cycle_weights=[(\"K\", -1), (\"L\", -1)]\n"
    );
}

#[test]
fn nonprime_coefficients_and_kernels_follow_the_actual_residue_embedding() {
    let output = stdout(interpret(
        "rvcmclasspoly H4[-4]; rvunram F[13,1]{2,0,1}; rvunram L[13,1]{5,0,1};
         rvec E[F](rvUnramElem[F]{1,1},rvUnramElem[F]{0}); rvecendorder M[E,H4];
         rvec H[L](rvUnramElem[L]{1,4},rvUnramElem[L]{0});
         rveccmcycle C[M,5,0]; rveccmaction A[M,C]; rveccmactionbase B[A,H];
         let p=rvECPoint[H](rvUnramElem[L]{2,0},rvUnramElem[L]{2,2});
         let ip=rvECPoint[H](rvUnramElem[L]{11,0},rvUnramElem[L]{10,10});
         assert rvECCMAct[B](rvECCMActionGenerator[A],p)==ip;
         assert rvECCMAct[B](rvECCMCycleElement[C],p)==add(rvECScale[2](p),ip);
         rvECCMActionData[B];",
    ));
    assert_eq!(output,
        "rvECCMActionData[B] = certificate=M, base=E, curve=H, index=1, discriminant=-4, full=true, generator=[0, 1], frobenius_weight=0, cycle_weights=[(\"C\", 1)]\n"
    );
}

#[test]
fn frobenius_over_a_nonprime_base_uses_its_full_residue_degree() {
    assert_eq!(
        stdout(interpret(
            "rvcmclasspoly H8[-8]; rvunram F[17,1]{14,0,1}; rvunram L[17,1]{14,0,0,0,1};
             rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{7}); rvecendorder M[E,H8];
             assert rvECEndFrobenius[M]==rvECEndElem[M,1,12];
             rvec EL[L](rvUnramElem[L]{1},rvUnramElem[L]{7});
             rveccmaction A[M]; rveccmactionbase B[A,EL];
             let p=rvECPoint[EL](rvUnramElem[L]{0,1},rvUnramElem[L]{1,0,6,9});
             let fp=rvECPoint[EL](rvUnramElem[L]{0,16},rvUnramElem[L]{1,0,6,8});
             assert rvECCMAct[B](rvECEndFrobenius[M],p)==fp;",
        )),
        ""
    );
}

#[test]
fn rejects_unwitnessed_elements_nominal_mismatches_and_wrong_extensions() {
    for (query, expected) in [
        ("show rvECCMAct[Frob](rvECEndElem[M,0,1],rvECInfinity[E]);", "element lies outside the realized suborder"),
        ("rvecendorder N[E,H8]; show rvECCMAct[Full](rvECEndElem[N,0,1],rvECInfinity[E]);", "CM action endomorphism"),
        ("rvecendorder N[E,H8]; rveccmcycle C[N,3,0]; rveccmaction A[M,K,C];", "same named endomorphism-order certificate"),
        ("rvec G[F](rvUnramElem[F]{1},rvUnramElem[F]{7}); show rvECCMAct[Full](rvECCMActionGenerator[Full],rvECInfinity[G]);", "CM action point"),
        ("rvec G[F](rvUnramElem[F]{1},rvUnramElem[F]{6}); rveccmactionbase B[Full,G];", "not the declared CM base extension"),
        ("rvunram G[13,1]{0,1}; rvec B[G](rvUnramElem[G]{1},rvUnramElem[G]{7}); rveccmactionbase BB[Full,B];", "residue embedding requires an extension"),
        ("show rvECCMActionMap[Full](rvECEndElem[M,0,-170141183460469231731687303715884105728]);", "expanded CM action degree exceeds the supported integer range"),
        ("rveccmaction Full[M];", "CM action is already defined"),
        ("rveccmaction A[M,K,Missing];", "unknown CM isogeny cycle 'Missing'"),
    ] {
        let error = diagnostic(interpret(&format!("{BASE} {query}")));
        assert!(error.contains(expected), "{error}");
    }
}
