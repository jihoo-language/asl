mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn kloosterman_factor_and_wild_two_pole_newton_polygon() {
    assert_eq!(
        stdout(example("rank-as-two-pole")),
        concat!(
            "rvASTwoGenus[Kloosterman] = 4\n",
            "rvASTwoPrank[Kloosterman] = 4\n",
            "rvASTwoCharacterL[Kloosterman, 1] = [[1, 0, 0, 0], [2, 0, 1, 1], [5, 0, 0, 0]]\n",
            "rvASTwoLPolynomial[Kloosterman] = [1, 6, 31, 96, 261, 480, 775, 750, 625]\n",
            "rvASTwoNewton[Kloosterman] = [0: 4, 1: 4]\n",
            "rvASTwoPoints[Kloosterman, 1] = 12\n",
            "rvASTwoPoints[Kloosterman, 3] = 72\n",
            "rvASTwoGenus[Wild] = 12\n",
            "rvASTwoPrank[Wild] = 4\n",
            "rvASTwoCharacterL[Wild, 1] = [[1, 0, 0, 0], [2, 0, 2, 0], [-1, -6, -1, -2], [-4, -8, -2, -6], [-25, -30, -20, -25], [-50, 0, 0, -50], [-125, 0, 0, 0]]\n",
            "rvASTwoCharacterNewton[Wild, 1] = [0: 1, 1/4: 1, 1/2: 2, 3/4: 1, 1: 1]\n",
            "rvASTwoNewton[Wild] = [0: 4, 1/4: 4, 1/2: 8, 3/4: 4, 1: 4]\n",
            "rvASTwoPoints[Wild, 6] = 16172\n",
            "rvASTwoGenus[Pure] = 14\n",
            "rvASTwoGenus[Cubic] = 14\n",
            "rvASTwoGenus[Quartic] = 14\n",
            "rvASTwoPrank[Pure] = 4\n",
            "rvASTwoPrank[Cubic] = 4\n",
            "rvASTwoPrank[Quartic] = 4\n",
            "rvASTwoCharacterNewton[Pure, 1] = [0: 1, 1/2: 5, 1: 1]\n",
            "rvASTwoCharacterNewton[Cubic, 1] = [0: 1, 3/8: 2, 1/2: 1, 5/8: 2, 1: 1]\n",
            "rvASTwoCharacterNewton[Quartic, 1] = [0: 1, 1/4: 1, 1/2: 3, 3/4: 1, 1: 1]\n",
            "rvASTwoPoints[Pure, 7] = 77497\n",
            "rvASTwoPoints[Cubic, 7] = 78407\n",
            "rvASTwoPoints[Quartic, 7] = 77457\n",

        )
    );
}

#[test]
fn pole_reversal_scaling_reduction_and_constant_twist_at_all_coefficients() {
    let mut source = String::from(
        r"
        rvastwo A[5] {(4,1),(-2,1)};
        rvastwo Flip[5] {(-4,1),(2,1)};
        rvastwo Scale[5] {(4,1),(-2,-1)};
        rvastwo Reduce[5] {(20,1),(-50,1),(-3,2),(-15,-2),(0,-5)};
        rvastwo Twist[5] {(4,1),(-2,1),(0,-1)};
        rvASTwoReduced[Reduce];
        rvASTwoReduced[Twist];
        ",
    );
    for k in 0..=6 {
        for j in 1..5 {
            for label in ["Flip", "Scale", "Reduce"] {
                writeln!(
                    &mut source,
                    "assert rvASTwoCoeff[A,{j},{k}] == rvASTwoCoeff[{label},{j},{k}];"
                )
                .unwrap();
            }
            writeln!(&mut source, "assert rvASTwoCoeff[Twist,{j},{k}] == rvCycloMul(rvCyclo[5,{}],rvASTwoCoeff[A,{j},{k}]);",-j*k).unwrap();
        }
    }
    source.push_str(
        r"
        assert rvASTwoCoeff[A,1,7] == zero[RVCyclo(5)];
        rvastwo K[5] {(1,1),(-1,1)};
        rvastwo KT[5] {(1,1),(-1,1),(0,1)};
        rvASTwoPoints[K,1];
        rvASTwoPoints[KT,1];
    ",
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvASTwoReduced[Reduce] = [(-2, 1), (4, 1)]\nrvASTwoReduced[Twist] = [(-2, 1), (0, 4), (4, 1)]\nrvASTwoPoints[K, 1] = 12\nrvASTwoPoints[KT, 1] = 2\n"
    );
}

#[test]
fn two_pole_validation_and_budget_errors() {
    for (statement, expected) in [
        ("rvastwo Bad[2] {(-1,1),(1,1)};", "requires an odd prime"),
        ("rvastwo Bad[9] {(-1,1),(1,1)};", "9 is not prime"),
        (
            "rvastwo Bad[3] {(-3,1),(-1,-1),(1,1)};",
            "both zero and infinity",
        ),
        ("rvASTwoPoints[C,0];", "extension degree must be positive"),
        (
            "rvASTwoCharacterNewton[C,3];",
            "character index must lie in 1..p-1",
        ),
        (
            "show rvASTwoCoeff[D,1,0];",
            "unknown two-pole Artin–Schreier cover",
        ),
        ("rvastwo C[3] {(-1,1),(1,1)};", "already defined"),
        (
            "rvastwo Large[5] {(-1,1),(21,1)}; rvASTwoCharacterL[Large,1];",
            "computation limit",
        ),
    ] {
        let source = format!("rvastwo C[3] {{(-1,1),(1,1)}};\n{statement}");
        let error = diagnostic(interpret(&source));
        assert!(error.contains(expected), "{statement}\n{error}");
    }
}
