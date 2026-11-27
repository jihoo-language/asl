mod common;

use common::{diagnostic, example, interpret, stdout};

const DELTA: &str = "
    rvheckepacket Delta(rvHeckeT[12,2]);
    rvheckeeigenperiod P[Delta,plus,2];
    rvmazurordinary O[P,11,4];
";

#[test]
fn delta_closed_disc_counts_and_hensel_enclosures() {
    let output = stdout(example("hecke-weight-disc-zeros"));
    assert_eq!(
        output,
        concat!(
            "rvMazurDiscZeros = precision=2, mu=1, lambda=1, reduction=[1, 2], residual_roots=[(5, 1)], roots_in_nonrational_residues=0\n",
            "RVMazurDiscRoot(O,2,0,0,0,1){x=5, x_precision=1, x_modulus=11, t=55, t_precision=2, t_modulus=121}\n",
            "rvMazurDiscZeros = precision=3, mu=1, lambda=1, reduction=[1, 2], residual_roots=[(5, 1)], roots_in_nonrational_residues=0\n",
            "RVMazurDiscRoot(O,3,0,0,0,1){x=115, x_precision=2, x_modulus=121, t=1265, t_precision=3, t_modulus=1331}\n",
            "rvMazurDiscZeros = precision=3, mu=2, lambda=1, reduction=[2, 2], residual_roots=[(10, 1)], roots_in_nonrational_residues=0\n",
            "RVMazurDiscRoot(O,3,0,0,55,2){x=10, x_precision=1, x_modulus=11, t=1265, t_precision=3, t_modulus=1331}\n",
            "rvMazurDiscZeros = precision=3, mu=1, lambda=0, reduction=[9], residual_roots=[], roots_in_nonrational_residues=0\n",
            "rvMazurDiscZeros = precision=1, unresolved\n",
        )
    );
}

fn group_coefficients(line: &str) -> Vec<u64> {
    line.split_once("{[")
        .unwrap()
        .1
        .strip_suffix("]}")
        .unwrap()
        .split(", ")
        .map(|value| value.parse().unwrap())
        .collect()
}

fn evaluate_group_polynomial(coefficients: &[u64], t: u64, modulus: u64) -> u64 {
    let mut power = 1;
    let mut sum = 0;
    for &coefficient in coefficients {
        sum = (sum + coefficient * power) % modulus;
        power = (power * (1 + t)) % modulus;
    }
    sum
}

#[test]
fn independent_finite_evaluation_and_precision_at_the_tail() {
    let program = format!(
        "{DELTA}
        let a = rvMazurAmice[O,3,0,0];
        show a; show rvMazurAmiceProject(a);
        let d = rvMazurDisc[0,1](a);
        rvMazurDiscCoeff[0](d); rvMazurDiscCoeff[1](d);
        rvMazurDiscCoeff[2](d); rvMazurDiscCoeff[3](d);
        rvMazurDiscCoeff[9000](d);
        let shifted = rvMazurDisc[-66,2](a);
        show rvMazurDiscRoot[0](shifted);
        assert rvMazurDiscRestrict[15467380314588111975607936701444009612,1](rvMazurDisc[-11,1](a)) ==
               rvMazurDisc[170141183460469231731687303715884105721,2](a);
        rvMazurDiscCoeff[0](rvMazurDiscRestrict[0,1](d));
        rvMazurDiscCoeff[1](rvMazurDiscRestrict[0,1](d));
        rvMazurDiscZeros(rvMazurDisc[0,1](rvMazurAmice[O,3,0,2]));
        rvmazurordinary L[P,11,1];
        rvMazurDiscZeros(rvMazurDisc[0,1](rvMazurAmice[L,3,0,0]));
        card (\\r : RVMazurDiscRoot(O,1,0,0,0,1) => S.hi);"
    );
    let output = stdout(interpret(&program));
    let lines: Vec<_> = output.lines().collect();
    let high = group_coefficients(lines[0]);
    let low = group_coefficients(lines[1]);
    for (row, modulus, expected_root) in [(high, 1331, 1265), (low, 121, 55)] {
        let zeros: Vec<_> = (0..modulus)
            .step_by(11)
            .filter(|&t| evaluate_group_polynomial(&row, t, modulus) == 0)
            .collect();
        assert_eq!(zeros, [expected_root]);
    }
    assert_eq!(
        &lines[2..],
        [
            "rvMazurDiscCoeff[0] = precision=3, modulus=1331, residue=253",
            "rvMazurDiscCoeff[1] = precision=3, modulus=1331, residue=869",
            "rvMazurDiscCoeff[2] = precision=3, modulus=1331, residue=1210",
            "rvMazurDiscCoeff[3] = precision=3, modulus=1331, residue=0",
            "rvMazurDiscCoeff[9000] = precision=3, modulus=1331, residue=0",
            "RVMazurDiscRoot(O,3,0,0,-66,2){x=0, x_precision=1, x_modulus=11, t=1265, t_precision=3, t_modulus=1331}",
            "rvMazurDiscCoeff[0] = precision=4, modulus=14641, residue=1584",
            "rvMazurDiscCoeff[1] = precision=4, modulus=14641, residue=9559",
            "rvMazurDiscZeros = precision=3, mu=0, lambda=0, reduction=[10], residual_roots=[], roots_in_nonrational_residues=0",
            "rvMazurDiscZeros = precision=1, unresolved",
            "card = 0 / 0",
        ]
    );
}

#[test]
fn unit_period_chart_rescaling_and_split_quadratic_packets() {
    let program = format!(
        "{DELTA}
        rvheckeeigenperiod Q[Delta,plus,4];
        rvmazurordinary B[Q,11,4];
        let chart = rvMazurDisc[0,1](rvMazurAmice[B,3,0,0]);
        rvMazurDiscZeros(chart);
        show rvMazurDiscRoot[5](chart);

        let t = rvHeckeT[24,2];
        rvheckepacket E(t);
        rvheckepacket F(add(rvHeckeScale[5](t),rvHeckeScale[7](rvHeckeUnit[24])));
        rvheckeeigenperiod U[E,plus,2]; rvheckeeigenperiod V[F,plus,2];
        rvheckeresidual S[E,31,28]; rvheckeresidual W[F,31,23];
        rvmazurordinary C[U,S,4]; rvmazurordinary D[V,W,4];
        let a = rvMazurAmice[C,3,0,26];
        let c = rvMazurDisc[0,1](a);
        let d = rvMazurDisc[0,1](rvMazurAmice[D,3,0,26]);
        rvMazurDiscZeros(c); rvMazurDiscZeros(d);
        rvMazurDiscCoeff[0](c); rvMazurDiscCoeff[0](d);
        show rvMazurDiscRoot[11](c); show rvMazurDiscRoot[11](d);
        show a;"
    );
    let output = stdout(interpret(&program));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        &lines[..2],
        [
            "rvMazurDiscZeros = precision=3, mu=1, lambda=1, reduction=[7, 3], residual_roots=[(5, 1)], roots_in_nonrational_residues=0",
            "RVMazurDiscRoot(B,3,0,0,0,1){x=115, x_precision=2, x_modulus=121, t=1265, t_precision=3, t_modulus=1331}",
        ]
    );
    assert_eq!(lines[2], "rvMazurDiscZeros = precision=3, mu=1, lambda=1, reduction=[4, 25], residual_roots=[(11, 1)], roots_in_nonrational_residues=0");
    assert_eq!(lines[2], lines[3]);
    assert_eq!(
        lines[4],
        "rvMazurDiscCoeff[0] = precision=3, modulus=29791, residue=10695"
    );
    assert_eq!(lines[4], lines[5]);
    for (line, ordinary) in [(lines[6], "C"), (lines[7], "D")] {
        assert_eq!(line, format!("RVMazurDiscRoot({ordinary},3,0,26,0,1){{x=755, x_precision=2, x_modulus=961, t=23405, t_precision=3, t_modulus=29791}}"));
    }
    let row = group_coefficients(lines[8]);
    let modulus = 29_791;
    let zeros: Vec<_> = (0..modulus)
        .step_by(31)
        .filter(|&t| evaluate_group_polynomial(&row, t, modulus) == 0)
        .collect();
    assert_eq!(zeros, [23_405]);
}

#[test]
fn invalid_charts_roots_and_unresolved_discs_have_distinct_diagnostics() {
    for (source, expected) in [
        ("show rvMazurDisc[1,1](a);", "center must belong to pZ_p"),
        ("show rvMazurDisc[0,0](a);", "depth must be at least one"),
        (
            "show rvMazurDisc[0,20](a);",
            "radius exceeds the unsigned 64-bit",
        ),
        (
            "show rvMazurDiscRestrict[0,0](d);",
            "refinement must be at least one",
        ),
        (
            "show rvMazurDiscRestrict[170141183460469231731687303715884105727,1](d);",
            "center exceeds the signed 128-bit",
        ),
        (
            "show rvMazurDiscRoot[11](d);",
            "root residue must be between 0 and 10",
        ),
        ("show rvMazurDiscRoot[4](d);", "residue 4 has no root"),
        (
            "show rvMazurDiscRoot[5](rvMazurDisc[0,1](rvMazurAmiceProject(a)));",
            "unresolved at the available precision",
        ),
        ("rvMazurDiscZeros(a);", "RVMazurDisc(O,n,r,i,a,e)"),
        ("show rvMazurDisc[0,1](d);", "RVMazurAmice(O,n,r,i)"),
    ] {
        let program = format!(
            "{DELTA}
            let a = rvMazurAmice[O,2,0,0];
            let d = rvMazurDisc[0,1](a);
            {source}"
        );
        let issue = diagnostic(interpret(&program));
        assert!(issue.contains(expected), "{issue}");
    }
}
