mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const ODD_LEVEL: &str = include_str!("../examples/hyperelliptic-generalized-kummer-cohomology.asl");
const EVEN_LEVEL: &str = include_str!("../examples/hyperelliptic-generalized-kummer.asl");

#[test]
fn example_reports_odd_cohomology_and_true_negative_traces() {
    assert_eq!(
        stdout(example("hyperelliptic-generalized-kummer-cohomology")),
        concat!(
            "rvGenKummerBetti[3, 2] = 7\n",
            "rvGenKummerBetti[3, 3] = 8\n",
            "rvGenKummerBetti[3, 4] = 108\n",
            "rvJacHilbBetti[2, 3] = 32\n",
            "rvGenKummerTrace[B, C, 2, 2] = 175\n",
            "rvGenKummerTrace[B, C, 2, 3] = -1000\n",
            "rvGenKummerTrace[B, C, 2, 4] = 67500\n",
            "rvJacHilbTrace[C, 2, 2, 3] = -4000\n",
            "rvGenKummerLevel[B, C, 2] = 593676, smooth generalized Kummer fiber, length=3, rational torsion order=81\n",
        )
    );
}

#[test]
fn geometric_betti_numbers_match_k3_fourfold_and_higher_euler_characteristics() {
    let mut source = String::new();
    for n in 1..=5 {
        for i in 0..=4 * (n - 1) {
            writeln!(source, "rvGenKummerBetti[{n},{i}];").unwrap();
        }
    }
    for i in 0..=8 {
        writeln!(source, "rvJacHilbBetti[2,{i}];").unwrap();
    }
    let values = values(&stdout(interpret(&source)));
    assert_eq!(values[0], 1);
    assert_eq!(&values[1..6], &[1, 0, 22, 0, 1]);
    assert_eq!(&values[6..15], &[1, 0, 7, 8, 108, 8, 7, 0, 1]);
    assert_eq!(alternating(&values[15..28]), 448);
    assert_eq!(alternating(&values[28..45]), 750);
    assert_eq!(&values[45..], &[1, 4, 13, 32, 44, 32, 13, 4, 1]);
}

#[test]
fn prime_level_traces_sum_to_count_and_display_the_middle_torsion_contribution() {
    let mut source = ODD_LEVEL.to_string();
    for r in 1..=2 {
        for i in 0..=8 {
            writeln!(source, "rvGenKummerTrace[B,C,{r},{i}];").unwrap();
        }
        writeln!(source, "rvGenKummerLevel[B,C,{r}];").unwrap();
        for i in 0..=8 {
            writeln!(source, "rvJacHilbTrace[C,2,{r},{i}];").unwrap();
        }
        writeln!(source, "rvJacHilbPoints[C,2,{r}];").unwrap();
    }
    let output = stdout(interpret(&source));
    let numbers = values(&output.lines().skip(9).collect::<Vec<_>>().join("\n"));
    for group in numbers.chunks_exact(20) {
        assert_eq!(alternating(&group[..9]), group[9]);
        assert_eq!(alternating(&group[10..19]), group[19]);
    }
    assert_eq!(numbers[2], 10 + 5);
    assert_eq!(numbers[3], 2 * 5 * 0);
    assert_eq!(
        &numbers[20..29],
        &[1, 0, 175, -1000, 67500, -25000, 109375, 0, 390625]
    );
    // At scalar Frobenius the untwisted middle contribution has dimension
    // 28. Every one of the 80 nonzero 3-torsion points adds 25^2.
    assert_eq!(numbers[24] - 28 * 25 * 25, 80 * 25 * 25);
}

#[test]
fn level_two_matches_the_resolved_k3_traces_over_both_fields() {
    let mut source = EVEN_LEVEL.to_string();
    for r in 1..=2 {
        for i in 0..=4 {
            writeln!(
                source,
                "rvGenKummerTrace[B,C,{r},{i}]; rvKummerHilbTrace[C,1,{r},{i}];"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let numbers = values(&output.lines().skip(12).collect::<Vec<_>>().join("\n"));
    assert_eq!(numbers.len(), 20);
    for pair in numbers.chunks_exact(2) {
        assert_eq!(pair[0], pair[1]);
    }
    assert_eq!((numbers[4], numbers[14]), (14, 166));
}

#[test]
fn handles_range_stability_and_overflow_without_huge_allocations() {
    let output = stdout(interpret(&format!(
        r"{ODD_LEVEL}
        rvGenKummerBetti[1,1]; rvGenKummerBetti[340282366920938463463374607431768211455,2];
        rvGenKummerBetti[2,340282366920938463463374607431768211455];
        rvJacHilbBetti[0,0]; rvJacHilbBetti[0,1];
        rvGenKummerTrace[B,C,340282366920938463463374607431768211455,9];
        rvJacHilbTrace[C,0,340282366920938463463374607431768211455,0];
        "
    )));
    assert_eq!(values(&output)[9..], [0, 7, 0, 1, 0, 0, 1]);
    for (statement, expected) in [
        (
            "rvGenKummerBetti[0,0];",
            "generalized Kummer length must be positive",
        ),
        (
            "rvGenKummerBetti[10000,10000];",
            "coefficient storage limit",
        ),
        (
            "rvGenKummerTrace[B,C,0,2];",
            "extension degree must be positive",
        ),
        (
            "rvGenKummerTrace[B,C,30,8];",
            "Frobenius trace exceeds signed 128-bit range",
        ),
        (
            "rvGenKummerTrace[B,C,1000,8];",
            "exact intermediate exceeds 512 bits",
        ),
        (
            "rvhyperelliptic E[F]{oF,oF,zF,zF,zF,oF}; rvGenKummerTrace[B,E,1,9];",
            "not the declared base extension",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{ODD_LEVEL}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn values(output: &str) -> Vec<i128> {
    output
        .lines()
        .map(|line| {
            line.split_once(" = ")
                .unwrap()
                .1
                .split(',')
                .next()
                .unwrap()
                .parse()
                .unwrap()
        })
        .collect()
}

fn alternating(values: &[i128]) -> i128 {
    values
        .iter()
        .enumerate()
        .map(|(i, &value)| if i.is_multiple_of(2) { value } else { -value })
        .sum()
}
