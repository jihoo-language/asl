mod common;

use common::{diagnostic, example, interpret, stdout};

const DEMO: &str = "fricke-semistable-monodromy";

#[test]
fn actual_nodes_have_the_verified_widths_gram_and_smith_groups() {
    let out = stdout(example(DEMO));
    assert!(out.contains("X0(14) at bad p=7; two P1 components; 2 geometric nodes; toric rank 1"));
    assert!(out.contains("edge widths = [1, 2]; exceptional edge indices after resolution = [1]"));
    assert!(out.contains("Gram = [[3]]; Smith = [3]; component order = 3"));
    assert!(out.contains("X0(46) at bad p=23; two P1 components; 6 geometric nodes; toric rank 5"));
    assert!(out.contains("Smith = [1, 1, 1, 1, 11]; component order = 11"));
    assert!(out.contains("X0(82) at bad p=41; two P1 components; 10 geometric nodes; toric rank 9"));
    assert!(out.contains("edge widths = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; exceptional edge indices after resolution = []"));
    assert!(out.contains("Smith = [1, 1, 1, 1, 1, 1, 1, 1, 10]; component order = 10"));
    assert!(out.contains("X0(10) at bad p=5; two P1 components; 1 geometric nodes; toric rank 0"));
    assert!(out.contains("Gram = []; Smith = []; component order = 1"));
}

#[test]
fn monodromy_and_closed_node_frobenius_preserve_the_expected_linking_forms() {
    let out = stdout(example(DEMO));
    assert!(out.contains("N = [[0, 3], [0, 0]], Φ = [[1, 0], [0, 7]]; N²=0; NΦ=pΦN"));
    assert!(out.contains("bad Euler factor (ascending) = [1, -1]; conductor exponent = 1"));
    assert!(out.contains(
        "bad Euler factor (ascending) = [1, -5, 10, -10, 5, -1]; conductor exponent = 5"
    ));
    assert!(out.contains(
        "bad Euler factor (ascending) = [1, -3, 0, 8, -6, -6, 8, 0, -3, 1]; conductor exponent = 9"
    ));
    assert!(out.contains("bad Euler factor (ascending) = [1]; conductor exponent = 0"));
    assert!(out.contains("component linking (mod Z) = 1/3"));
    assert!(out.contains("component linking (mod Z) = 5/11"));
    assert!(out.contains("component linking (mod Z) = 9/10"));
}

#[test]
fn ordered_hecke_spectra_and_integral_component_bridges_match_independent_computations() {
    let out = stdout(example(DEMO));
    assert!(out.contains("ordered Γ0(2) Hecke coefficient [3,1,1] = -16777216"));
    assert!(out.contains("ordered Γ0(2) Hecke coefficient [5,4,4] = 905720"));
    assert!(
        out.contains("oriented Γ0(2) T3 on cycles = [[-2]]; characteristic (descending) = [1, 2]")
    );
    assert!(out.contains("integral component bridge = [[-2]]"));
    for (coefficients, action) in [
        ("[1, 0, -10, 0, 25, 0]", "4 mod 11"),
        ("[1, 0, -20, 0, 80, -64]", "6 mod 11"),
        ("[1, 2, -10, -16, 40, 40, -76, -24, 56, -16]", "4 mod 10"),
        (
            "[1, 6, -4, -80, -80, 288, 464, -224, -640, -256]",
            "6 mod 10",
        ),
    ] {
        assert!(out.contains(&format!(
            "characteristic (descending) = {coefficients}; full Manin J+ and Fricke quotient actions certified; component action = {action}"
        )));
    }
}

#[test]
fn nominal_cycle_dual_and_component_diagnostics_are_checked() {
    let prefix = "rvmanin M[7]; rvunram F[7,1]{-3,0,1}; rvfrickess B[F]; rvfrickebad S[M,B];";
    for (body, expected) in [
        ("show rvFrickeCycle[S]{1,2};", "requires 1 coordinates"),
        (
            "show rvFrickeMonodromy(rvFrickeDual[S]{1});",
            "has the wrong input lattice",
        ),
        (
            "show rvFrickeComponent[S,3];",
            "Smith residue is out of range",
        ),
        ("show rvFrickeBadNode[S,2];", "node index is out of range"),
        (
            "show rvFrickeBadAtkin[3](rvFrickeCycle[S]{1});",
            "index must be 2 or p",
        ),
        (
            "show rvFrickeBadHecke[7](rvFrickeCycle[S]{1});",
            "prime to the bad characteristic",
        ),
        (
            "rvfrickebad R[M,B]; rvFrickeLink(rvFrickeComponent[S,1],rvFrickeComponent[R,1]);",
            "rvFrickeLink",
        ),
        (
            "rvmanin M2[11]; rvfrickebad R[M2,B];",
            "level must be twice the Brandt characteristic",
        ),
    ] {
        let result = diagnostic(interpret(&format!("{prefix} {body}")));
        assert!(result.contains(expected), "{result}");
    }
}
