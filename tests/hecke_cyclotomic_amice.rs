mod common;

use common::{diagnostic, example, interpret, stdout};

const DELTA: &str = "
    rvheckepacket Delta(rvHeckeT[12,2]);
    rvheckeeigenperiod P[Delta,plus,2];
    rvmazurordinary O[P,11,4];
";

#[test]
fn eleven_adic_delta_coefficients_and_primitive_character() {
    assert_eq!(
        stdout(example("hecke-cyclotomic-amice")),
        concat!(
            "rvMazurAmiceCoeff[0] = precision=4, modulus=14641, residue=1584\n",
            "rvMazurAmiceCoeff[1] = precision=1, modulus=11, residue=2\n",
            "rvMazurAmiceCoeff[1] = precision=2, modulus=121, residue=79\n",
            "rvMazurAmiceIwasawa[4] = mu=0, lambda=1\n",
            "rvMazurAmiceIwasawa[4] = mu=0, lambda=0\n",
            "rvMazurAmiceCharacter[1] = order=11, modulus=14641, polynomial=[14618, 12285, 4190, 3860, 7995, 7368, 3300, 8129, 2870, 7259]\n",
        )
    );

    let program = format!(
        "{DELTA}
         let a = rvMazurAmice[O,2,0,0];
         assert add(a,neg(a)) == zero[RVMazurAmice(O,2,0,0)];
         assert rvMazurAmiceGroupCoeff[11](a) == rvMazurAmiceGroupCoeff[0](a);
         show rvMazurAmiceGroupCoeff[0](a);
         rvMazurAmiceCoeff[2](a);
         rvMazurAmiceCoeff[3](a);
         rvMazurAmiceCoeff[4](a);
         rvMazurAmiceCoeff[5](a);
         rvMazurAmiceCoeff[0](rvMazurAmice[O,2,0,2]);
         let weighted = rvMazurAmice[O,2,2,8];
         rvMazurAmiceCoeff[0](weighted);
         rvMazurAmiceIwasawa[1](weighted);
         rvMazurAmiceCharacter[0](a);
         rvMazurAmiceIwasawa[0](a);
         rvMazurAmiceIwasawa[25](rvMazurAmice[O,2,0,1]);"
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "RVMazurOrdinary(O){13552}\n",
            "rvMazurAmiceCoeff[2] = precision=1, modulus=11, residue=10\n",
            "rvMazurAmiceCoeff[3] = precision=1, modulus=11, residue=7\n",
            "rvMazurAmiceCoeff[4] = precision=1, modulus=11, residue=1\n",
            "rvMazurAmiceCoeff[5] = precision=1, modulus=11, residue=3\n",
            "rvMazurAmiceCoeff[0] = precision=4, modulus=14641, residue=13342\n",
            "rvMazurAmiceCoeff[0] = precision=4, modulus=14641, residue=8162\n",
            "rvMazurAmiceIwasawa[1] = mu=0, lambda=1\n",
            "rvMazurAmiceCharacter[0] = order=1, modulus=14641, polynomial=[1584]\n",
            "rvMazurAmiceIwasawa[0] = mu=unknown, lambda=unknown, first_unit_degree_at_least=1\n",
            "rvMazurAmiceIwasawa[25] = mu=unknown, lambda=unknown, first_unit_degree_at_least=11\n",
        )
    );
}

#[test]
fn levels_prove_only_the_available_t_adic_digits() {
    let source = format!(
        "{DELTA}
         let high = rvMazurAmice[O,3,0,0];
         let low = rvMazurAmiceProject(high);
         rvMazurAmiceCoeff[11](high);
         rvMazurAmiceCoeff[11](low);
         rvMazurAmiceCoeff[120](high);
         rvMazurAmiceCoeff[121](high);
         rvMazurAmiceCoeff[1](rvMazurAmiceProject(low));
         rvMazurAmiceCharacter[1](high);
         rvMazurAmiceCharacter[1](low);"
    );
    let report = stdout(interpret(&source));
    let lines: Vec<_> = report.lines().collect();
    assert_eq!(
        &lines[..5],
        [
            "rvMazurAmiceCoeff[11] = precision=1, modulus=11, residue=4",
            "rvMazurAmiceCoeff[11] = precision=0, unknown",
            "rvMazurAmiceCoeff[120] = precision=1, modulus=11, residue=1",
            "rvMazurAmiceCoeff[121] = precision=0, unknown",
            "rvMazurAmiceCoeff[1] = precision=0, unknown",
        ]
    );
    assert_eq!(lines[5], lines[6]);
}

#[test]
fn an_odd_period_chart_change_and_split_quadratic_embeddings() {
    let source = format!(
        "{DELTA}
         rvheckeeigenperiod Q[Delta,minus,1];
         rvheckeeigenperiod R[Delta,minus,5];
         rvmazurordinary A[Q,11,4]; rvmazurordinary B[R,11,4];
         let odd = rvMazurAmice[A,2,1,0]; let other = rvMazurAmice[B,2,1,0];
         rvMazurAmiceCoeff[0](odd); rvMazurAmiceCoeff[0](other);
         rvMazurAmiceIwasawa[4](odd);
         assert rvMazurAmiceProject(odd) == rvMazurAmice[A,1,1,0];
         let t = rvHeckeT[24,2];
         rvheckepacket E(t);
         rvheckepacket F(add(rvHeckeScale[5](t),rvHeckeScale[7](rvHeckeUnit[24])));
         rvheckeeigenperiod U[E,plus,2]; rvheckeeigenperiod V[F,plus,2];
         rvheckeresidual S[E,13,2]; rvheckeresidual W[F,13,4];
         rvmazurordinary C[U,S,5]; rvmazurordinary D[V,W,5];
         let c = rvMazurAmice[C,2,2,4]; let d = rvMazurAmice[D,2,2,4];
         assert rvMazurAmiceProject(c) == rvMazurAmice[C,1,2,4];
         rvMazurAmiceCharacter[1](c); rvMazurAmiceCharacter[1](d);"
    );
    let report = stdout(interpret(&source));
    let lines: Vec<_> = report.lines().collect();
    assert_eq!(
        &lines[..3],
        [
            "rvMazurAmiceCoeff[0] = precision=4, modulus=14641, residue=2309",
            "rvMazurAmiceCoeff[0] = precision=4, modulus=14641, residue=7889",
            "rvMazurAmiceIwasawa[4] = mu=0, lambda=0",
        ]
    );
    assert_eq!(lines[3], lines[4]);
    assert!(!lines[3].ends_with("polynomial=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]"));
}

#[test]
fn invalid_levels_characters_and_nominal_grades() {
    for (source, expected) in [
        ("show rvMazurAmice[O,0,0,0];", "level must be at least one"),
        (
            "show rvMazurAmice[O,2,0,10];",
            "character index must be between 0 and 9",
        ),
        ("show rvMazurAmice[O,2,11,0];", "between 0 and k−2"),
        (
            "show rvMazurAmice[O,19,0,0];",
            "exceeds the unsigned 64-bit range",
        ),
        (
            "show rvMazurAmice[O,7,0,0];",
            "exceeds the finite table limit",
        ),
        (
            "show rvMazurAmiceProject(rvMazurAmice[O,1,0,0]);",
            "cannot project below",
        ),
        (
            "rvMazurAmiceCharacter[2](rvMazurAmice[O,2,0,0]);",
            "conductor exceeds",
        ),
        (
            "rvMazurAmiceCoeff[0](rvMazurTateMeasure[O,1,0,1]);",
            "RVMazurAmice(O,n,r,i)",
        ),
        (
            "rvmazurordinary Other[P,11,4];
             show add(rvMazurAmice[O,1,0,0],rvMazurAmice[Other,1,0,0]);",
            "RVMazurAmice(O,1,0,0)",
        ),
    ] {
        let result = diagnostic(interpret(&format!("{DELTA}{source}")));
        assert!(result.contains(expected), "{result}");
    }
}
