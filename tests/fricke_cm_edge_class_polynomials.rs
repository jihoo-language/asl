mod common;

use common::{diagnostic, example, interpret, stdout};

const DEMO: &str = "fricke-cm-edge-class-polynomials";

#[test]
fn exact_norms_retain_automorphism_multiplicities_and_match_both_resultants() {
    let out = stdout(example(DEMO));
    assert!(out.contains("rvFrickeCMPolyData[N3] = degree=3, coefficients=[15252992, 184512, 744, 1], order=-3, directed multiplicity strata=[(3, [248, 1])]"));
    assert!(out.contains("order=-4, directed multiplicity strata=[(1, [104, 1]), (2, [-544, 1])]"));
    assert!(out.contains("order=-7, directed multiplicity strata=[(1, [4073, 1]), (2, [23, 1])]"));
    assert_eq!(
        out.matches("quadratic norm and Bareiss resultant agree")
            .count(),
        4
    );
}

#[test]
fn exact_interorder_edges_keep_both_orientations() {
    let out = stdout(example(DEMO));
    assert!(out.contains("degree=1, coefficients=[248, 1], orders=(-3, -12), directed degrees=[3, 1], directed multiplicity strata=[(3, 1, [248, 1])]"));
    assert!(out.contains("coefficients=[-544, 1], orders=(-4, -16), directed degrees=[2, 1]"));
    assert!(out.contains("coefficients=[4073, 1], orders=(-7, -28), directed degrees=[1, 1]"));
    assert!(out.contains("degree=2, coefficients=[13399504, -53504, 1], orders=(-12, -48), directed degrees=[2, 2], directed multiplicity strata=[(1, 1, [13399504, -53504, 1])]"));
    assert!(out.contains("rvFrickeCMPolyCoeff[E1248, 0] = 13399504"));
}

#[test]
fn horizontal_distinct_and_both_kinds_of_loop_are_certified() {
    let out = stdout(example(DEMO));
    assert!(out.contains("orders=(-4, -4), directed degrees=[1, 1], directed multiplicity strata=[(1, 1, [104, 1])], distinct endpoints=[1], Fricke-fixed loops=[104, 1], nonfixed paired loops=[1]"));
    assert!(out.contains("orders=(-7, -7), directed degrees=[2, 2], directed multiplicity strata=[(2, 2, [23, 1])], distinct endpoints=[1], Fricke-fixed loops=[1], nonfixed paired loops=[23, 1]"));
    assert!(out.contains("Fricke-fixed loops=[-152, 1], nonfixed paired loops=[1]"));
    assert!(out.contains("orders=(-63, -63), directed degrees=[8, 8]"));
    assert!(out.contains("distinct endpoints=[17529534886321441, -49488637933732, 67515454374, -733, 1], Fricke-fixed loops=[1], nonfixed paired loops=[1]"));
}

#[test]
fn finite_endpoint_inseparability_labels_and_witt_arithmetic_are_exact() {
    let out = stdout(example(DEMO));
    assert!(out.contains("order-compatible factors (-12, -48)=[[[2], [1]], [[7], [1]]]"));
    assert!(out.contains("repeated endpoint, order-compatible factors (-12, -48)=[[[10], [1]], [[10], [1]]], residue order labels overlap"));
    assert!(out.contains("endpoint polynomial=[[1], [4], [1]], factors (degree, multiplicity, coefficients)=[(2, 1, [[1], [4], [1]])], nonsplit, Frobenius exchanges the roots"));
    assert!(out.contains("rvFrickeCMPolyFactors[E1248] = factors (degree, multiplicity, coefficients)=[(1, 2, [[6], [1]])]"));
    assert!(out.contains("rvFrickeCMPolyFactors[N7] = factors (degree, multiplicity, coefficients)=[(1, 3, [[3], [1]])]"));
    let binary = stdout(interpret(
        "rvunram F[2,1]{0,1}; rvFrickeCMEndpoints(rvUnramElem[F]{});",
    ));
    assert!(binary.contains("repeated endpoint, characteristic-two polynomial reduction"));
}

#[test]
fn reduction_collisions_and_wrong_domains_do_not_create_an_edge() {
    let out = stdout(example(DEMO));
    assert!(out.contains("rvFrickeCMPolyData[Empty] = degree=0, coefficients=[1], orders=(-3, -7), directed degrees=[0, 0]"));
    for (source, fragment) in [
        (
            "rvcmclasspoly A[-3]; rvcmclasspoly B[-7]; rvfrickecmedge E[A,B]; rvunram F[5,1]{0,1}; rvFrickeCMEdgeEndpoints[E](rvUnramElem[F]{2});",
            "outside the reduced characteristic-zero Fricke edge support",
        ),
        (
            "rvcmclasspoly A[-3]; rvfrickecmnorm R[A]; rvunram F[5,1]{0,1}; rvFrickeCMEdgeEndpoints[R](rvUnramElem[F]{2});",
            "requires a declared edge polynomial",
        ),
        (
            "rvunram F[5,2]{0,1}; rvFrickeCMEndpoints(rvUnramElem[F]{});",
            "requires a residue field at precision one",
        ),
        (
            "rvfrickecmnorm R[Unknown];",
            "unknown CM class polynomial",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(fragment), "{error}");
    }
}
