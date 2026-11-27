mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_frobenius_solutions_and_mirror_inverse_match_reference() {
    let output = stdout(example("dwork-quartic-mirror-periods"));
    assert!(output.contains(
        "f=[1, 24, 2520, 369600, 63063000, 11732745024], g=[0, 104, 12276, 5632160/3, 327270650, 308431885104/5], h=[0, 0, 5408, 1146912, 688332760/3, 46687269200]"
    ));
    assert!(output.contains("u=[1, 12, 1188, 170544, 28779300, 5318414640]"));
    assert!(output.contains("W=[1, 128, 24576, 5242880,"));
    assert!(
        output.contains("q(z)/z=[1, 104, 15188, 2585184, 480222434, 94395247376, 19304070595976]")
    );
    assert!(output.contains("z(q)/q=[1, -104, 6444, -311744, 13018830, -493025760, 17411253944]"));
    assert!(output.contains("certified both compositions"));
    assert!(output.contains("certified Ore equations and symmetric-square identities"));
    let endpoints = stdout(interpret("rvDworkMirrorMap[0]; rvDworkMirrorMap[12];"));
    assert!(endpoints.contains("rvDworkMirrorMap[0] = log(q/z)=[0], q(z)/z=[1], z(q)/q=[1]"));
    assert!(endpoints
        .lines()
        .last()
        .unwrap()
        .ends_with("certified both compositions"));
}

#[test]
fn actual_ore_elements_apply_through_the_logarithmic_jet() {
    let output = stdout(example("dwork-quartic-mirror-periods"));
    assert!(output.contains("rvOreCoeff = -176"));
    for (rank, sol) in [(3, 0), (3, 1), (3, 2), (2, 0), (2, 1)] {
        assert!(output.contains(&format!("rvDworkMirrorApply[{rank}, {sol}, 5] = 0")));
    }
    assert!(output.contains("rvDworkMirrorApply[2, 1, 2] = coefficients of 1, log z, (log z)^2/2 = [[1, 116, 23244], [0, 12, 2376], [0, 0, 0]]"));
    let output = stdout(interpret(
        "let L=rvDworkPicardFuchs[];
         let z=rvOre[1,0];
         let theta=rvOre[0,1];
         let bracket=add(rvOreMul(theta,z),neg(rvOreMul(z,theta)));
         assert bracket == z;
         rvDworkMirrorApply[3,2,4](rvOreMul(theta,L));
         rvDworkMirrorApply[3,0,2](bracket);
         rvDworkMirrorApply[3,2,0](rvOre[0,2]);",
    ));
    assert!(output.contains("rvDworkMirrorApply[3, 2, 4] = 0"));
    assert!(output.contains("rvDworkMirrorApply[3, 0, 2] = coefficients of 1, log z, (log z)^2/2 = [[0, 1, 24], [0, 0, 0], [0, 0, 0]]"));
    assert!(output.contains(
        "rvDworkMirrorApply[3, 2, 0] = coefficients of 1, log z, (log z)^2/2 = [[1], [0], [0]]"
    ));
}

#[test]
fn log_shifts_and_exact_horizontal_form_are_checked() {
    let output = stdout(example("dwork-quartic-mirror-periods"));
    assert!(output.contains("rvDworkMirrorMonodromy[5, 1] = [[1, 1, 1/2], [0, 1, 1], [0, 0, 1]]"));
    assert!(output.contains("A=D^-1*[[0, 1-256z, 0], [0, 0, 1-256z], [24z, 176z, 384z]]"));
    assert!(output
        .contains("B=[[-48z, -128z, 1-256z], [-128z, -1+256z, 0], [1-256z, 0, 0]], det B=D^3"));
    assert!(output.contains("certified symmetric square and horizontal pairing"));
    let output = stdout(interpret(
        "rvDworkMirrorMonodromy[8,-3]; rvDworkMirrorMonodromy[8,0];",
    ));
    assert!(output.contains("[[1, -3, 9/2], [0, 1, -3], [0, 0, 1]]"));
    assert!(output.contains("[[1, 0, 0], [0, 1, 0], [0, 0, 1]]"));
}

#[test]
fn brauer_connection_certifies_mod_p_and_reports_actual_higher_difference() {
    let output = stdout(example("dwork-quartic-mirror-periods"));
    assert!(output.contains("mod 7, z=4, terms=2, b_7=6, hypergeometric sum=6, scaled sum=6, difference=0, agreeing digits=1"));
    assert!(output.contains("mod 49, z=11, terms=13, b_49=22, hypergeometric sum=22, scaled sum=43, difference=28, agreeing digits=1"));
    assert!(output.contains("b_49=22, b_7=20, unit=6"));
    let other = stdout(interpret(
        "rvDworkMirrorCongruence[13,2,3]; rvDworkMirrorCongruence[13,4,2];",
    ));
    assert!(other.contains("mod 2197, z=1821, terms=550, b_2197=1292, hypergeometric sum=1201"));
    assert!(other.contains("difference=1196, agreeing digits=1"));
    assert!(other
        .contains("b_169=0, hypergeometric sum=0, scaled sum=0, difference=0, agreeing digits=2"));
}

#[test]
fn domain_and_operator_work_limits_are_enforced() {
    for (source, message) in [
        ("rvDworkMirrorMap[13];", "degree must be at most 12"),
        (
            "rvDworkMirrorApply[4,0,1](rvOre[0,0]);",
            "rank must be 2 or 3",
        ),
        (
            "rvDworkMirrorApply[2,2,1](rvOre[0,0]);",
            "solution index is outside",
        ),
        (
            "rvDworkMirrorApply[3,2,1](rvOre[0,100000000]);",
            "exceeds the computation limit",
        ),
        (
            "rvunram F[5,1]{0,1}; rvDworkMirrorApply[3,2,1](rvUnramElem[F]{1});",
            "requires an Ore operator",
        ),
        ("rvDworkMirrorCongruence[7,0,2];", "psi nonzero modulo p"),
        ("rvDworkMirrorCongruence[7,1,2];", "singular Dwork quartic"),
        (
            "rvDworkMirrorCongruence[7,3,0];",
            "precision must be positive",
        ),
    ] {
        let failure = diagnostic(interpret(source));
        assert!(failure.contains(message), "{failure}");
    }
}
