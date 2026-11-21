mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn ic_trace_basis_and_rank_three_adjoint_tensor_square() {
    assert_eq!(
        stdout(example("rv-satake-ic")),
        concat!(
            "rvSphCoeff = [1, 1]\n",
            "rvSatakeCoeff = {0: 2}\n",
            "rvSatakeCoeff = {4: 2}\n",
            "rvSatakeAt = 20070400\n",
        ),
    );
}

#[test]
fn satake_inverse_reconstructs_signed_coefficients_and_det_shifts() {
    let program = r"
        let first = rvSphBasis[[-4,-1,2]];
        let second = rvSphBasis[[-5,0,1]];
        let third = rvSphBasis[[-2,-2,-2]];
        let q = rvSphQ[3];
        let a = add(rvSphMul(q,first),add(neg(second),third));
        assert rvSatakeInverse(rvSatake(a)) == a;
        assert rvSatakeInverse(zero[RVSatake(3)]) == zero[RVSphHecke(3)];
        assert rvSatakeInverse(rvSatV[0,2]) == rvSphQ[0];
        assert rvSatSchur[[]] == rvSatV[0,0];
        assert rvSphIC[[]] == rvSphId[0];
        let negative = rvSatSchur[[-3,-1]];
        rvSatakeAt[2,[2,3]](negative);
        assert rvSatake(rvSphIC[[-3,-1]]) == rvSatakeMul(rvSatV[2,2],negative);
        assert rvSphIC[[-3,-1]] == rvSphMul(rvSphBasis[[-3,-3]],rvSphIC[[0,2]]);
    ";
    assert_eq!(stdout(interpret(program)), "rvSatakeAt = 19/216\n");
}

#[test]
fn inverse_does_not_need_to_negate_the_lowest_signed_leading_coefficient() {
    let mut program = "let c0 = neg(rvSatV[0,0]);\n".to_owned();
    for number in 1..=127 {
        let previous = number - 1;
        program.push_str(&format!("let c{number} = add(c{previous},c{previous});\n"));
    }
    program.push_str("rvSphCoeff[[]](rvSatakeInverse(c127));");
    assert_eq!(
        stdout(interpret(&program)),
        format!("rvSphCoeff = [{}]\n", i128::MIN),
    );
}

#[test]
fn partial_inverse_checks_integral_q_image_and_valid_coweights() {
    for (program, expected) in [
        (
            "show rvSatakeInverse(rvSphId[2]);",
            "rvSatakeInverse requires a Satake polynomial",
        ),
        (
            "show rvSatakeInverse(rvSatV[2,-2]);",
            "requires nonnegative v powers",
        ),
        (
            "show rvSatakeInverse(rvSatV[2,1]);",
            "requires even v powers",
        ),
        (
            "show rvSatakeInverse(rvSatSchur[[0,1]]);",
            "after dividing by v^1 (found v^0)",
        ),
        (
            "show rvSatakeInverse(rvSatV[1,170141183460469231731687303715884105726]);",
            "q degree is out of range",
        ),
        (
            "show rvSatSchur[[1,0]];",
            "coweight must be weakly increasing",
        ),
        (
            "show rvSphIC[[0,1000000]];",
            "Satake Laurent calculation exceeds the storage limit",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
