mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F7: &str = include_str!("../examples/hyperelliptic-k3-pencil-bsd-factors.asl");
const F23: &str = include_str!("../examples/hyperelliptic-k3-pencil-selmer.asl");

fn f7() -> &'static str {
    F7.split_once("rvK3PencilBSDFactors[P]").unwrap().0
}

fn f23() -> &'static str {
    F23.split_once("rvK3PencilSelmer[P]").unwrap().0
}

#[test]
fn rank_one_separates_factors_with_or_without_a_marked_generator() {
    let output = stdout(example("hyperelliptic-k3-pencil-bsd-factors"));
    let provenance = "rank=1, torsion=(Z/2)^2, prime-to-characteristic certificate=I0* at infinity, characteristic certificate=I2 at t=[6], height denominator=4, dim Sha[2]=0, Sha*Reg=1";
    assert_eq!(
        output,
        format!(
            "rvK3PencilBSDFactors[P] = {provenance}, Sha order divides=4, Sha order=1, regulator=1\n\
             rvK3MordellBSDFactors[M, 5, 1] = {provenance}, zero trope=4, section height=1, generator index=1, regulator=1, Sha order=1\n\
             rvK3MordellBSDFactors[M, 5, -3] = {provenance}, zero trope=4, section height=9, generator index=3, regulator=1, Sha order=1\n"
        ),
    );

    let mut source = format!("{} rvK3PencilLFunction[P]; rvK3MordellLattice[M];", f7());
    for mask in [5, 6, 12, 14] {
        writeln!(source, "rvK3MordellBSDFactors[M,{mask},2];").unwrap();
    }
    source.push_str("rvk3mordell Q[P,5]; rvK3MordellBSDFactors[Q,4,1];");
    let output = stdout(interpret(&source));
    assert!(output.contains("trivial rank=17, Mordell-Weil rank=1"));
    assert!(output.contains("trivial rank=16, trivial |det|=1024, marked rank=17"));
    assert_eq!(
        output
            .matches("section height=4, generator index=2, regulator=1, Sha order=1")
            .count(),
        4
    );
    assert!(output
        .contains("zero trope=5, section height=1, generator index=1, regulator=1, Sha order=1"));
}

#[test]
fn exact_signed_multiples_include_the_smallest_signed_integer() {
    let output = stdout(interpret(&format!(
        "{} rvK3MordellBSDFactors[M,5,{}];",
        f7(),
        i128::MIN,
    )));
    assert!(output.contains("section height=28948022309329048855892746252171976963317496166410141009864396001978282409984"));
    assert!(output.contains(
        "generator index=170141183460469231731687303715884105728, regulator=1, Sha order=1"
    ));
}

#[test]
fn rank_two_keeps_nonunique_factors_and_reports_exact_arithmetic_constraints() {
    let f13 = include_str!("../examples/hyperelliptic-k3-pencil-two-descent.asl")
        .split_once("rvK3PencilTorsion2[P]")
        .unwrap()
        .0;
    for (setup, expected, certificate) in [
        (
            f23(),
            "height denominator=4, dim Sha[2]=2, Sha*Reg=4, Sha order divides=64, Sha order candidates=[4, 16, 64]",
            "prime-to-characteristic certificate=I0* at t=[0], characteristic certificate=I4 at t=[21]",
        ),
        (
            f13,
            "height denominator=2, dim Sha[2]=0, Sha*Reg=9/4, Sha order divides=9, Sha order candidates=[1, 9]",
            "prime-to-characteristic certificate=I0* at infinity, characteristic certificate=I2 at t=[12]",
        ),
    ] {
        let output = stdout(interpret(&format!("{setup} rvK3PencilBSDFactors[P];")));
        assert!(output.contains("rank=2, torsion=(Z/2)^2"));
        assert!(output.contains(certificate), "{output}");
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn rank_two_can_separate_the_factors_when_only_one_square_order_survives() {
    let setup = include_str!("../examples/hyperelliptic-k3-pencil-two-descent.asl")
        .split_once("rvK3PencilTorsion2[P]")
        .unwrap()
        .0
        .replace("F[13,1]", "F[5,1]")
        .replace("z,eleven,two,nine,three,o", "z,neg(o),z,z,z,o");
    let output = stdout(interpret(&format!("{setup} rvK3PencilBSDFactors[P];")));
    assert!(output.contains("rank=2, torsion=(Z/2)^2"));
    assert!(output.contains("Sha*Reg=5/4, Sha order divides=20, Sha order=1, regulator=5/4"));
}

#[test]
fn wrong_rank_zero_height_and_nominal_inputs_have_specific_diagnostics() {
    for (setup, operation, expected) in [
        (
            f23(),
            "rvk3mordell M[P,0]; rvK3MordellBSDFactors[M,3,1];",
            "requires arithmetic Mordell-Weil rank one",
        ),
        (
            f7(),
            "rvK3MordellBSDFactors[M,7,1];",
            "section must have positive height",
        ),
        (
            f7(),
            "rvK3MordellBSDFactors[M,5,0];",
            "section multiple must be nonzero",
        ),
        (
            f7(),
            "rvK3MordellBSDFactors[M,0,1];",
            "section trope must contain exactly one",
        ),
        (f7(), "rvK3MordellBSDFactors[P,5,1];", "unknown K3 Mordell"),
        (f7(), "rvK3PencilBSDFactors[M];", "unknown K3 pencil"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{operation}")));
        assert!(error.contains(expected), "{error}");
    }
}
