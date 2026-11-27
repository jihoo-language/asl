mod common;

use common::{diagnostic, example, interpret, stdout};

const DEMO: &str = "fricke-brandt-hecke-orders";
const F23: &str = "rvunram F[23,1]{-5,0,1}; rvfrickess B[F];";
const ORDER: &str = "rvfrickeheckeorder A[B,0]{[-4,1],[0,1],[-5,0,1]};";

#[test]
fn exact_crt_smith_and_local_congruences_agree_with_independent_f23_results() {
    let out = stdout(example(DEMO));
    assert!(out.contains("cokernel Smith = [1, 1, 1, 55]; index = 55"));
    assert!(out.contains(
        "factor (ascending) = [-4, 1]; certified irreducible; rank = 1; projector denominator = 11"
    ));
    assert!(out.contains(
        "factor (ascending) = [0, 1]; certified irreducible; rank = 1; projector denominator = 5"
    ));
    assert!(out.contains("factor (ascending) = [-5, 0, 1]; certified irreducible; rank = 2; projector denominator = 55"));
    assert!(out.contains("B5 on block 2 = [-1, -1] / 1"));
    assert!(out.contains("prime 11; T eigenvalue = 4; residue root = 4; exact joint congruence depth = 1; residual valuations: minimal polynomial=1, B3=∞, B5=1, Frob=∞"));
    assert!(out.contains("scalar block ambient eigenvector = [1, 3, 6, 1]\nalgebraic block ambient eigenvector = [1, 3, 6, 1]"));
    assert!(out
        .contains("prime 5; T eigenvalue = 4; residue root = 4; exact joint congruence depth = 0"));
    assert!(out.contains("RVFrickeSSQDiv(B23){[4, 1, 2, 4] / 11}"));
}

#[test]
fn f41_joint_congruence_depth_distinguishes_hecke_from_the_generator_alone() {
    let out = stdout(example(DEMO));
    assert!(out.contains("cokernel Smith = [1, 1, 1, 1, 10]; index = 10"));
    assert!(out.contains(
        "factor (ascending) = [-4, 1]; certified irreducible; rank = 1; projector denominator = 5"
    ));
    assert!(out.contains(
        "factor (ascending) = [2, 1]; certified irreducible; rank = 1; projector denominator = 2"
    ));
    assert!(out.contains("factor (ascending) = [2, -4, 0, 1]; certified irreducible; rank = 3; projector denominator = 10"));
    assert!(out.contains("B5 on block 2 = [2, 0, -1] / 1"));
    assert!(out.contains("Frob on block 1 = [-1] / 1"));
    assert!(out.contains("Frob on block 2 = [1, 0, 0] / 1"));
    assert!(out.contains("exact joint congruence depth = 1; residual valuations: minimal polynomial=2, B3=∞, B5=1, Frob=∞"));
    assert!(out.contains("scalar block ambient eigenvector = [1, 1, 1, 1, 1]\nalgebraic block ambient eigenvector = [1, 1, 1, 1, 1]"));
    assert!(out.contains("prime 2; T eigenvalue = -2; residue root = 0; exact joint congruence depth = 1; residual valuations: minimal polynomial=1, B3=∞, B5=∞, Frob=1"));
    assert!(out.contains("RVFrickeSSQDiv(B41){[0, 1, -1, 0, 0] / 2}"));
}

#[test]
fn mixed_operators_and_grouped_or_repeated_factors_keep_their_true_block_structure() {
    let out = stdout(example(DEMO));
    assert!(out.contains("T = B3 + (2) B5; saturated-block cokernel Smith = [1, 1, 1, 55]"));
    assert!(out.contains("B3 on block 2 = [-2, -1] / 1"));
    assert!(out.contains("B5 on block 2 = [1, 1] / 1"));
    assert!(out.contains(
        "factor (ascending) = [1, 1]; certified irreducible; rank = 2; projector denominator = 55"
    ));
    assert!(out.contains("block 2 has multiplicity; a single polynomial image is uncertified"));
    let source = format!(
        "{F23} {ORDER}
         rvfrickeheckeorder C[B,0]{{[-5,0,1],[0,-4,1]}};
         rvFrickeOrderData[C];
         let v=rvFrickeSSBasis[B,1];
         assert rvFrickeOrderProject[C,0](v)==rvFrickeOrderProject[A,2](v);"
    );
    let out = stdout(interpret(&source));
    assert!(out.contains("factor (ascending) = [0, -4, 1]; certified reducible; rank = 2"));
}

#[test]
fn projections_cancel_in_wide_arithmetic_and_extract_typed_numerators() {
    let max = i128::MAX;
    let source = format!(
        "{F23} {ORDER}
         let e=rvFrickeSSBasis[B,0];
         assert rvFrickeOrderNumerator(rvFrickeOrderProject[A,0](e))==rvFrickeSSDiv[B]{{4,1,2,4}};
         let v=rvFrickeSSDiv[B]{{{max},-{max},0,0}};
         assert rvFrickeOrderIntegral(rvFrickeOrderProject[A,0](v))==zero[RVFrickeSSDiv(B)];"
    );
    stdout(interpret(&source));
    let failure = format!(
        "{F23} {ORDER} show rvFrickeOrderIntegral(rvFrickeOrderProject[A,0](rvFrickeSSBasis[B,0]));"
    );
    assert!(
        diagnostic(interpret(&failure)).contains("Brandt projection has nontrivial denominator 11")
    );
}

#[test]
fn bad_factor_nominal_and_bad_level_diagnostics() {
    for (body, message) in [
        (
            "rvfrickeheckeorder A[B,0]{[-4,1]};",
            "do not annihilate the chosen operator",
        ),
        (
            "rvfrickeheckeorder A[B,1]{[-10,1],[-4,1],[1,1],[1,1]};",
            "not the minimal polynomial",
        ),
        (
            "rvfrickeheckeorder A[B,0]{[-8,2],[0,1],[-5,0,1]};",
            "positive-degree monic integral polynomials",
        ),
        (
            "rvfrickeheckeorder A[B,0]{[-4,1],[0,1],[-5,0,1]}; rvfrickess X[F]; show rvFrickeOrderProject[A,0](rvFrickeSSBasis[X,0]);",
            "Brandt eigenprojection",
        ),
        (
            "rvfrickeheckeorder A[B,0]{[-4,1],[0,1],[-5,0,1]}; rvFrickeOrderCongruence[A,0,1,5];",
            "requires a rational eigenline and a certified higher-degree eigenfield",
        ),
        (
            "rvfrickeheckeorder A[B,0]{[-4,1],[0,1],[-5,0,1]}; rvFrickeOrderCongruence[A,0,2,9];",
            "9 is not prime",
        ),
    ] {
        let result = diagnostic(interpret(&format!("{F23} {body}")));
        assert!(result.contains(message), "{result}");
    }
    let good = stdout(interpret(
        "rvunram F[5,1]{-2,0,1}; rvfrickess B[F]; rvfrickeheckeorder A[B,0]{[-4,1]}; rvFrickeOrderData[A];",
    ));
    assert!(good.contains("T = B3; saturated-block cokernel Smith = [1]; index = 1"));
    assert!(!good.contains("B5"));
    let bad = diagnostic(interpret(
        "rvunram F[5,1]{-2,0,1}; rvfrickess B[F]; rvfrickeheckeorder A[B,1]{[-4,1]};",
    ));
    assert!(bad.contains("c=0 in characteristic five"));
}

#[test]
fn higher_rank_fraction_free_elimination_keeps_exact_smith_invariants() {
    for (source, result) in [
        (
            "rvunram F[83,1]{-2,0,1}; rvfrickess B[F];
             rvfrickeheckeorder A[B,0]{[-4,1],[1,1],[-4,2,1],[-25,-4,30,5,-10,-1,1]};
             rvFrickeOrderData[A];",
            "cokernel Smith = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 10742]; index = 21484",
        ),
        (
            "rvunram F[89,1]{-3,0,1}; rvfrickess B[F];
             rvfrickeheckeorder A[B,0]{[-4,1],[1,1],[-2,1],[-1,2,1],[-1,-9,-16,-4,3,1]};
             rvFrickeOrderData[A];",
            "cokernel Smith = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1540]; index = 6160",
        ),
    ] {
        let out = stdout(interpret(source));
        assert!(out.contains(result), "{out}");
    }
}
