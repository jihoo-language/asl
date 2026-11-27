mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn runnable_example_verifies_constant_and_nonconstant_lifts() {
    assert_eq!(
        stdout(example("elliptic-serre-tate-canonical-lifts")),
        concat!(
            "RVUnram(W)[1886]\n",
            "RVUnram(W)[1722]\n",
            "RVUnram(W)[481]\n",
            "RVUnram(U)[287, 71]\n",
            "RVUnram(U)[275, 272]\n",
        )
    );
}

#[test]
fn truncation_and_both_certified_frobenius_conjugates() {
    let output = stdout(interpret(
        "rvcmclasspoly H3[-3]; rvcmclasspolylift H12[H3,2]; rvcmclasspolylift H48[H12,2];
         rvunram F[7,1]{5,5,1}; rvunram W2[7,2]{5,5,1}; rvunram W3[7,3]{5,5,1};
         rvunram Alt[7,2]{12,5,1};
         let a=rvUnramElem[F]{6,4}; let b=rvUnramElem[F]{4,5};
         rvec E[F](a,b);
         rvec S[F](rvUnramFrob[1](a),rvUnramFrob[1](b));
         rvecendorder C[E,H48]; rvecendorder D[S,H48];
         let j=rvECEndCanonicalJ[C,W3];
         let conjugate=rvECEndCanonicalJ[D,W3];
         assert rvECEndCanonicalJ[C,F]==rvUnramElem[F]{0,1};
         assert rvECEndCanonicalJ[C,W2]==rvUnramElem[W2]{42,22};
         assert rvECEndCanonicalJ[C,Alt]==rvUnramElem[Alt]{0,15};
         assert j==rvUnramElem[W3]{287,71};
         assert rvUnramFrob[1](j)==conjugate;
         assert rvUnramFrob[1](conjugate)==j;
         assert conjugate==rvUnramElem[W3]{275,272};
         let c=rvUnramElem[W3]{-2835810000};
         let d=rvUnramElem[W3]{6549518250000};
         assert add(add(rvUnramMul(j,j),rvUnramMul(c,j)),d)==zero[RVUnram(W3)];
         assert add(add(rvUnramMul(conjugate,conjugate),rvUnramMul(c,conjugate)),d)==zero[RVUnram(W3)];
         assert rvECEndCanonicalPhi[C,W3]==zero[RVUnram(W3)];
         assert rvECEndCanonicalPhi[D,W3]==zero[RVUnram(W3)];
         assert rvECEndUnitRoot[C,F]==rvUnramElem[F]{5};
         assert rvECEndNonunitRoot[C,F]==zero[RVUnram(F)];
         assert rvECEndUnitRoot[C,W2]==rvUnramElem[W2]{47};
         assert rvECEndNonunitRoot[C,W2]==zero[RVUnram(W2)];
         assert rvECEndUnitRoot[C,W3]==rvUnramElem[W3]{194};
         assert rvECEndNonunitRoot[C,W3]==rvUnramElem[W3]{147};",
    ));
    assert!(output.is_empty(), "{output}");
}

#[test]
fn embedding_is_multiplicative_and_conjugation_gives_trace_and_norm() {
    let output = stdout(interpret(
        "rvcmclasspoly H4[-4]; rvcmclasspolylift H16[H4,2];
         rvunram F[13,1]{0,1}; rvunram W2[13,2]{0,1}; rvunram W[13,3]{0,1};
         rvec E[F](rvUnramElem[F]{2},rvUnramElem[F]{1});
         rvecendorder O[E,H16];
         assert rvECEndCanonicalJ[O,F]==rvUnramElem[F]{1};
         assert rvECEndCanonicalJ[O,W2]==rvUnramElem[W2]{27};
         assert rvECEndCanonicalJ[O,W]==rvUnramElem[W]{1886};
         let x=rvECEndElem[O,7,2];
         let y=rvECEndElem[O,1,-3];
         let a=rvECEndEmbed[O,W](x);
         let b=rvECEndEmbedConj[O,W](x);
         assert a==rvUnramElem[W]{1248};
         assert b==rvUnramElem[W]{963};
         assert b==rvECEndEmbed[O,W](rvECEndConj(x));
         assert add(a,b)==rvUnramElem[W]{14};
         assert rvUnramMul(a,b)==rvUnramElem[W]{65};
         assert rvECEndEmbed[O,W](rvECEndMul(x,y))==rvUnramMul(a,rvECEndEmbed[O,W](y));
         assert rvECEndEmbed[O,W](add(x,y))==add(a,rvECEndEmbed[O,W](y));",
    ));
    assert!(output.is_empty(), "{output}");
}

#[test]
fn integral_modular_evaluation_uses_more_than_the_residue_congruence() {
    let output = stdout(interpret(
        "rvunram A[2,5]{0,1}; rvunram B[3,3]{0,1};
         assert rvModularPolyWittEval[2](rvUnramElem[A]{},rvUnramElem[A]{1})==rvUnramElem[A]{17};
         assert rvModularPolyWittEval[3](rvUnramElem[B]{},rvUnramElem[B]{1})==rvUnramElem[B]{10};
         rvcmclasspoly H4[-4]; rvunram F[5,1]{0,1}; rvunram W[5,3]{0,1};
         rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{});
         rvecendorder O[E,H4];
         assert rvECEndCanonicalJ[O,W]==rvUnramElem[W]{1728};
         assert rvECEndCanonicalPhi[O,W]==zero[RVUnram(W)];",
    ));
    assert!(output.is_empty(), "{output}");
}

#[test]
fn presentation_checks_and_nominal_embedding_errors() {
    let setup = "rvcmclasspoly H3[-3]; rvcmclasspolylift H12[H3,2]; rvcmclasspolylift H48[H12,2];
         rvunram F[7,1]{5,5,1}; rvunram W[7,3]{5,5,1};
         rvunram Wrong[7,3]{1,0,1}; rvunram WrongPrime[11,3]{1,0,1};
         rvec E[F](rvUnramElem[F]{6,4},rvUnramElem[F]{4,5});
         rvecendorder C[E,H48]; rvecendorder D[E,H48];";
    for (query, expected) in [
        (
            "show rvECEndCanonicalJ[C,Wrong];",
            "curve's residue presentation",
        ),
        (
            "show rvECEndUnitRoot[C,WrongPrime];",
            "curve's residue presentation",
        ),
        (
            "show rvECEndCanonicalJ[Missing,W];",
            "unknown elliptic endomorphism-order certificate",
        ),
        (
            "show rvECEndEmbed[C,W](rvECEndElem[D,0,1]);",
            "expected type RVECEnd(C)",
        ),
        (
            "show rvModularPolyEval[7](rvUnramElem[W]{},rvUnramElem[W]{});",
            "requires a precision-one RVUnram field",
        ),
        (
            "show rvModularPolyWittEval[7](rvUnramElem[W]{},rvUnramElem[Wrong]{});",
            "expected type RVUnram(W)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {query}")));
        assert!(error.contains(expected), "{error}");
    }
}
