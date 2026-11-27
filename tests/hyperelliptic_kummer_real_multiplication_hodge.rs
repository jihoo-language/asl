mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn number(line: &str) -> u128 {
    line.rsplit_once(" = ").unwrap().1.parse().unwrap()
}

#[test]
fn example_detects_extra_real_multiplication_and_orientation_classes() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-real-multiplication-hodge")),
        concat!(
            "rvHilbHodgeInvariant[2, 4, 14] = 106\n",
            "rvHilbRealMultiplicationHodge[2, 4, 14, 3] = 108\n",
            "rvHilbRealMultiplicationPrimitiveHodge[2, 4, 14, 3] = 94\n",
            "rvHilbRealMultiplicationHodge[3, 4, 14, 3] = 122\n",
            "rvHilbSOOscillator[3, 6] = 8\n",
        )
    );
}

fn torus_counts(cap: usize, weights: &[i32]) -> Vec<Vec<[u128; 2]>> {
    let width = 2 * cap + 1;
    let mut counts = vec![vec![[0; 2]; width]; cap + 1];
    counts[0][cap][0] = 1;
    for part in 1..=cap {
        for &charge in weights {
            for n in part..=cap {
                for b in 0..width {
                    let previous = b as i32 - charge;
                    if (0..width as i32).contains(&previous) {
                        let [even, odd] = counts[n - part][previous as usize];
                        counts[n][b][0] += odd;
                        counts[n][b][1] += even;
                    }
                }
            }
        }
    }
    counts
}

fn so3_invariants(cap: usize) -> Vec<[u128; 2]> {
    // Every nontrivial integer-spin SO(3) irrep has one weight 0 and
    // one weight 1. Subtract, keeping even/odd tensor parity separate.
    torus_counts(cap, &[0, -1, 1])
        .into_iter()
        .map(|row| [row[cap][0] - row[cap + 1][0], row[cap][1] - row[cap + 1][1]])
        .collect()
}

#[test]
fn abstract_oscillators_match_rank_one_characters_past_the_first_gram_relation() {
    let cap = 24;
    let so3 = so3_invariants(cap);
    let so2 = torus_counts(cap, &[-1, 1]);
    let mut program = String::new();
    let mut expected = Vec::new();
    for n in 0..=cap {
        writeln!(program, "rvHilbSOOscillator[2,{n}];").unwrap();
        expected.push(so2[n][cap].iter().sum());
        writeln!(program, "rvHilbSOOscillator[3,{n}];").unwrap();
        expected.push(so3[n].iter().sum());
    }
    let got: Vec<_> = stdout(interpret(&program)).lines().map(number).collect();
    assert_eq!(got, expected);

    let mut free_gram = vec![0; cap + 1];
    free_gram[0] = 1;
    for i in 1..=cap {
        for j in i..=cap - i {
            for n in i + j..=cap {
                free_gram[n] += free_gram[n - i - j];
            }
        }
    }
    // O(3) fixes the even part. The first 4-by-4 Gram minor is in
    // degree 2(1+2+3+4)=20. SO(3) also fixes determinant orientations.
    assert!(so3[..6].iter().all(|pair| pair[1] == 0));
    assert_eq!(so3[6][1], 1);
    assert_eq!(
        so3[..20].iter().map(|pair| pair[0]).collect::<Vec<_>>(),
        free_gram[..20]
    );
    assert_eq!(so3[20][0] + 1, free_gram[20]);
}

#[test]
fn multiple_real_places_match_an_independent_graded_character_product() {
    let cap = 12;
    let b: Vec<u128> = so3_invariants(cap)
        .iter()
        .map(|pair| pair.iter().sum())
        .collect();
    for e in [2, 4, 7] {
        let rho = 23 - 3 * e;
        let mut weights = vec![0; rho - 1];
        weights.extend([-1, 1]);
        let mut charged = torus_counts(cap, &weights);
        for _ in 0..e {
            charged = (0..=cap)
                .map(|n| {
                    (0..=2 * cap)
                        .map(|z| {
                            [
                                (0..=n)
                                    .map(|j| b[j] * (charged[n - j][z][0] + charged[n - j][z][1]))
                                    .sum(),
                                0,
                            ]
                        })
                        .collect()
                })
                .collect();
        }
        let mut program = String::new();
        let mut expected = Vec::new();
        for &n in &[2, 3, 6, 9, 12] {
            for j in 0..=2 * n {
                let z = cap - n + j;
                let total = charged[n][z][0];
                let primitive = if j == 0 {
                    total
                } else if j <= n {
                    total - charged[n][z - 1][0]
                } else {
                    0
                };
                writeln!(
                    program,
                    "rvHilbRealMultiplicationHodge[{n},{},{rho},{e}];",
                    2 * j
                )
                .unwrap();
                writeln!(
                    program,
                    "rvHilbRealMultiplicationPrimitiveHodge[{n},{},{rho},{e}];",
                    2 * j
                )
                .unwrap();
                expected.extend([total, primitive]);
            }
        }
        let got: Vec<_> = stdout(interpret(&program)).lines().map(number).collect();
        assert_eq!(got, expected, "rho={rho}, e={e}");
    }
}

#[test]
fn field_degree_one_agrees_with_orthogonal_branching_in_each_grade() {
    for rho in [1, 2, 7, 15, 20] {
        let mut program = String::new();
        for n in 2..=9 {
            for j in 0..=2 * n {
                let d = 2 * j;
                writeln!(program, "rvHilbRealMultiplicationHodge[{n},{d},{rho},1];").unwrap();
                writeln!(program, "rvHilbHodgeInvariant[{n},{d},{rho}];").unwrap();
                writeln!(
                    program,
                    "rvHilbRealMultiplicationPrimitiveHodge[{n},{d},{rho},1];"
                )
                .unwrap();
                writeln!(program, "rvHilbPrimitiveHodgeInvariant[{n},{d},{rho}];").unwrap();
            }
        }
        let report = stdout(interpret(&program));
        let values: Vec<_> = report.lines().map(number).collect();
        assert!(!values.is_empty());
        for pair in values.chunks_exact(2) {
            assert_eq!(pair[0], pair[1], "rho={rho}");
        }
    }
}

#[test]
fn hilbert_square_and_cube_have_one_quadratic_form_per_real_place() {
    let mut program = String::new();
    let mut expected = Vec::new();
    for rho in 1u128..=20 {
        for e in 1..=(23 - rho) / 3 {
            if !(23 - rho).is_multiple_of(e) {
                continue;
            }
            let middle = rho * (rho + 1) / 2 + e;
            for (n, d, all, primitive) in [
                (2, 2, rho, rho - 1),
                (2, 4, middle, middle - rho),
                (3, 4, middle + rho, middle),
            ] {
                writeln!(program, "rvHilbRealMultiplicationHodge[{n},{d},{rho},{e}];").unwrap();
                writeln!(
                    program,
                    "rvHilbRealMultiplicationPrimitiveHodge[{n},{d},{rho},{e}];"
                )
                .unwrap();
                expected.extend([all, primitive]);
            }
        }
    }
    let got: Vec<_> = stdout(interpret(&program)).lines().map(number).collect();
    assert_eq!(got, expected);
}

#[test]
fn bounds_use_the_excess_grade_and_validate_the_given_group() {
    let report = stdout(interpret(
        r"
        rvHilbRealMultiplicationHodge[0,0,2,7];
        rvHilbRealMultiplicationPrimitiveHodge[0,0,2,7];
        rvHilbRealMultiplicationHodge[0,2,2,7];
        rvHilbRealMultiplicationHodge[2,3,2,7];
        rvHilbRealMultiplicationHodge[2,10,2,7];
        rvHilbRealMultiplicationPrimitiveHodge[2,6,2,7];
        rvHilbRealMultiplicationHodge[1000000000000000000,2,2,7];
        rvHilbRealMultiplicationHodge[1000000000000000000,3999999999999999998,2,7];
        rvHilbRealMultiplicationPrimitiveHodge[1000000000000000000,4,2,7];
        rvHilbSOOscillator[1000000000000000000,3];
        ",
    ));
    assert_eq!(
        report.lines().map(number).collect::<Vec<_>>(),
        [1, 1, 0, 0, 0, 0, 2, 2, 11, 1]
    );
    for (input, expected) in [
        (
            "rvHilbRealMultiplicationHodge[2,4,0,1];",
            "rank must be in 1..=21",
        ),
        (
            "rvHilbRealMultiplicationHodge[2,4,2,0];",
            "positive divisor",
        ),
        (
            "rvHilbRealMultiplicationHodge[2,4,2,2];",
            "positive divisor",
        ),
        (
            "rvHilbRealMultiplicationHodge[2,4,21,1];",
            "at least three transcendental directions",
        ),
        (
            "rvHilbRealMultiplicationHodge[1,2,2,7];",
            "whose LLV type is D12",
        ),
        ("rvHilbSOOscillator[1,1];", "dimension must be at least two"),
        (
            "rvHilbSOOscillator[3,1000000000000000000];",
            "polynomial computation limit",
        ),
        (
            "rvHilbRealMultiplicationHodge[1000000000000000000,1000000000000000000,2,7];",
            "polynomial computation limit",
        ),
    ] {
        assert!(diagnostic(interpret(input)).contains(expected), "{input}");
    }
}
