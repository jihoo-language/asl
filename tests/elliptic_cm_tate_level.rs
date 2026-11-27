mod common;

use common::{diagnostic, example, interpret, stdout};

const SETUP: &str = "rvcmclasspoly H8[-8]; rvunram F[17,1]{0,1};
    rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{7}); rvecendorder M[E,H8];
    rveccmcycle K[M,3,1]; rveccmaction A[M,K];
    rvunram L[17,1]{14,0,1}; rvec EL[L](rvUnramElem[L]{1},rvUnramElem[L]{7});
    rveccmactionbase AL[A,EL];
    let p4=rvECPoint[EL](rvUnramElem[L]{7,2},rvUnramElem[L]{3,8});
    let q4=rvECPoint[EL](rvUnramElem[L]{8,6},rvUnramElem[L]{6,6});
    let p12=rvECPoint[EL](rvUnramElem[L]{0,0},rvUnramElem[L]{0,5});
    let q12=rvECPoint[EL](rvUnramElem[L]{0,1},rvUnramElem[L]{2,1});";

#[test]
fn independent_level_four_and_composite_twelve_fixtures() {
    assert_eq!(
        stdout(example("elliptic-cm-tate-level")),
        concat!(
            "rvECLevelData[T4] = curve=EL, level=4, basis_pairing=[4, 0], frobenius_columns=[[1, 2], [2, 1]]\n",
            "rvECLevelData[T12] = curve=EL, level=12, basis_pairing=[2, 8], frobenius_columns=[[11, 0], [6, 7]]\n",
            "rvECLevelInvariants = modulus=12, trace=0, determinant=2\n",
        )
    );
}

#[test]
fn torsion_coordinates_ring_operations_and_projections_commute_with_maps() {
    let source = format!(
        "{SETUP}
         rveclevel B[AL,4](p4,q4); rveclevel C[AL,12](p12,q12);
         let w=rvECEndElem[M,0,1]; let pi=rvECEndFrobenius[M]; let z=rvECEndElem[M,1,-3];
         let wc=rvECLevelMatrix[C](w); let wb=rvECLevelMatrix[B](w);
         let pc=rvECLevelFrob[C]; let pb=rvECLevelFrob[B];
         assert rvECLevelMatrix[C](pi)==pc;
         assert rvECLevelMatrix[C](add(w,z))==add(wc,rvECLevelMatrix[C](z));
         assert rvECLevelMatrix[C](rvECEndMul(w,z))==rvECLevelMul(wc,rvECLevelMatrix[C](z));
         assert rvECLevelMul(wc,wc)==rvECLevelMat[C,-2,0,0,-2];
         assert rvECLevelMul(wc,pc)==rvECLevelMul(pc,wc);
         assert add(wc,neg(wc))==zero[RVECLevelMat(C)];
         assert add(rvECLevelVec[C,-1,2],rvECLevelVec[C,1,-2])==zero[RVECLevelVec(C)];
         assert (\\v:RVECLevelVec(C)=>rvECLevelCoord[C](rvECLevelPoint(v)))==(\\v:RVECLevelVec(C)=>v);
         let v=rvECLevelVec[C,5,2]; let b=rvECLevelProject[B](v);
         assert rvECLevelProject[B](rvECLevelApply(wc,v))==rvECLevelApply(wb,b);
         assert rvECLevelProject[B](rvECLevelApply(pc,v))==rvECLevelApply(pb,b);
         assert rvECLevelWeil[B](rvECScale[3](p12),rvECScale[3](q12))
             ==rvUnramPow[3](rvECLevelBasisPairing[C]);
         assert rvECLevelMatrix[C](add(w,rvECEndElem[M,0,12]))==wc;
         let huge=rvECEndElem[M,0,-170141183460469231731687303715884105728];
         assert rvECLevelMatrix[C](huge)==rvECLevelMat[C,4,0,0,8];
         let p=rvECLevelPoint(v); let fp=rvECLevelPoint(rvECLevelApply(pc,v));
         assert rvECX(fp)==rvUnramFrob[1](rvECX(p));
         assert rvECY(fp)==rvUnramFrob[1](rvECY(p));"
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn true_odd_pairing_and_ramified_two_torsion_with_common_support() {
    let source = format!(
        "{SETUP}
         rveclevel B[AL,2](rvECScale[6](p12),rvECScale[6](q12));
         let p=rvECLevelPoint(rvECLevelVec[B,1,0]); let q=rvECLevelPoint(rvECLevelVec[B,0,1]);
         assert rvECLevelWeil[B](p,q)==rvUnramElem[L]{{16}};
         assert rvECLevelWeil[B](p,add(p,q))==rvUnramElem[L]{{16}};
         assert rvECLevelWeil[B](p,p)==rvUnramElem[L]{{1}};
         assert rvECLevelWeil[B](p,rvECInfinity[EL])==rvUnramElem[L]{{1}};
         rvcmclasspoly H3[-3]; rvunram T[7,1]{{0,1}};
         rvec X[T](rvUnramElem[T]{{0}},rvUnramElem[T]{{2}}); rvecendorder N[X,H3];
         rveccmaction J[N];
         let u=rvECPoint[X](rvUnramElem[T]{{0}},rvUnramElem[T]{{3}});
         let v=rvECPoint[X](rvUnramElem[T]{{3}},rvUnramElem[T]{{1}});
         rveclevel C[J,3](u,v);
         assert rvECLevelBasisPairing[C]==rvUnramElem[T]{{4}};
         assert rvECLevelFrob[C]==rvECLevelMat[C,1,0,0,1];
         assert rvECLevelMatrix[C](rvECEndFrobenius[N])==rvECLevelFrob[C];
         assert rvECLevelWeil[C](v,u)==rvUnramElem[T]{{2}};
         assert rvECLevelWeil[C](u,add(u,v))==rvUnramElem[T]{{4}};
         assert rvECLevelWeil[C](u,neg(u))==rvUnramElem[T]{{1}};"
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn rejects_bad_levels_dependent_or_nonnominal_points_and_outside_suborders() {
    for (query, expected) in [
        ("rveclevel B[AL,1](p4,q4);", "at least two and prime to the characteristic"),
        ("rveclevel B[AL,17](p4,q4);", "at least two and prime to the characteristic"),
        ("rveclevel B[AL,2000000000](p4,q4);", "elliptic level basis table exceeds the supported size"),
        ("rveclevel B[AL,3](p4,q4);", "not killed by the level"),
        ("rveclevel B[AL,4](p4,rvECScale[2](p4));", "dependent modulo the level"),
        ("rveclevel B[AL,12](p12,rvECScale[2](q12));", "dependent modulo the level"),
        ("rveclevel B[AL,4](rvECInfinity[E],q4);", "elliptic level generator"),
        ("rveclevel B[AL,4](p4,q4); show rvECLevelWeil[B](p4,p12);", "outside the certified torsion level"),
        ("rveclevel B[AL,4](p4,q4); rveclevel C[AL,12](p12,q12); show rvECLevelProject[C](rvECLevelVec[B,1,0]);", "requires the same named curve and a dividing level"),
        ("rveclevel B[AL,4](p4,q4); rveclevel C[AL,4](p4,q4); show rvECLevelMul(rvECLevelFrob[B],rvECLevelFrob[C]);", "elliptic level matrix factor"),
        ("rveccmaction A0[M]; rveccmactionbase L0[A0,EL]; rveclevel B[L0,4](p4,q4); show rvECLevelMatrix[B](rvECEndElem[M,0,1]);", "element lies outside the realized suborder"),
    ] {
        let error = diagnostic(interpret(&format!("{SETUP} {query}")));
        assert!(error.contains(expected), "{error}");
    }
}
