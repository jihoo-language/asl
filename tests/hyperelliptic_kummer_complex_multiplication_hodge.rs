mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn number(line: &str) -> u128 {
    line.rsplit_once(" = ").unwrap().1.parse().unwrap()
}

#[test]
fn the_middle_hilbert_cube_detects_new_unitary_skew_forms() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-complex-multiplication-hodge")),
        concat!(
            "rvHilbRealMultiplicationHodge[2, 4, 15, 2] = 122\n",
            "rvHilbComplexMultiplicationHodge[2, 4, 15, 2] = 122\n",
            "rvHilbRealMultiplicationHodge[3, 6, 15, 2] = 816\n",
            "rvHilbComplexMultiplicationHodge[3, 6, 15, 2] = 818\n",
            "rvHilbComplexMultiplicationPrimitiveHodge[3, 6, 15, 2] = 681\n",
            "rvHilbGLOscillator[2, 3, 0] = 2\n",
            "rvHilbGLOscillator[2, 3, 1] = 1\n",
        )
    );
    let report = stdout(interpret(
        r"
        rvHilbComplexMultiplicationHodge[3,2,15,2];
        rvHilbComplexMultiplicationHodge[3,4,15,2];
        rvHilbRealMultiplicationPrimitiveHodge[3,6,15,2];
        ",
    ));
    assert_eq!(
        report.lines().map(number).collect::<Vec<_>>(),
        [15, 137, 679]
    );
}

fn torus(cap: usize, weights: &[(i32, i32)]) -> Vec<Vec<Vec<u128>>> {
    let width = 2 * cap + 1;
    let mut counts = vec![vec![vec![0; width]; width]; cap + 1];
    counts[0][cap][cap] = 1;
    for mode in 1..=cap {
        for &(a, b) in weights {
            for n in mode..=cap {
                for i in 0..width {
                    let x = i as i32 - a;
                    if !(0..width as i32).contains(&x) {
                        continue;
                    }
                    for j in 0..width {
                        let y = j as i32 - b;
                        if (0..width as i32).contains(&y) {
                            let previous = counts[n - mode][x as usize][y as usize];
                            counts[n][i][j] += previous;
                        }
                    }
                }
            }
        }
    }
    counts
}

#[test]
fn determinant_sectors_match_direct_gl1_and_gl2_weyl_extraction() {
    let cap = 18;
    let rank_one = torus(cap, &[(1, 0), (-1, 0)]);
    let rank_two = torus(cap, &[(1, 0), (-1, 0), (0, 1), (0, -1)]);
    let mut expected = Vec::new();
    let mut program = String::new();
    for k in 0..=cap {
        for s in -4..=4 {
            let x = (cap as i32 + s) as usize;
            writeln!(program, "rvHilbGLOscillator[1,{k},{s}];").unwrap();
            expected.push(rank_one[k][x][cap]);
            // Weyl's numerator subtracts the adjacent root shift.
            writeln!(program, "rvHilbGLOscillator[2,{k},{s}];").unwrap();
            expected.push(rank_two[k][x][x] - rank_two[k][x + 1][x - 1]);
        }
    }
    let got: Vec<_> = stdout(interpret(&program)).lines().map(number).collect();
    assert_eq!(got, expected);

    let mut free = vec![0; cap + 1];
    free[0] = 1;
    // GL pairings are ordered between a W mode and a W-dual mode.
    for a in 1..=cap {
        for b in 1..=cap - a {
            for k in a + b..=cap {
                free[k] += free[k - a - b];
            }
        }
    }
    let invariant = |k: usize| rank_two[k][cap][cap] - rank_two[k][cap + 1][cap - 1];
    let first_relation = 2 * (1 + 2 + 3);
    assert_eq!(
        (0..first_relation).map(invariant).collect::<Vec<_>>(),
        free[..first_relation]
    );
    assert_eq!(invariant(first_relation) + 1, free[first_relation]);
}

#[test]
fn rank_two_transcendental_cm_matches_the_full_hodge_diagonal() {
    for n in 2..=9 {
        let mut program = format!("rvHilbLLVHodgeTable[{n}];\n");
        for p in 0..=2 * n {
            writeln!(
                program,
                "rvHilbComplexMultiplicationHodge[{n},{},21,1];",
                2 * p
            )
            .unwrap();
            writeln!(
                program,
                "rvHilbComplexMultiplicationPrimitiveHodge[{n},{},21,1];",
                2 * p
            )
            .unwrap();
        }
        let report = stdout(interpret(&program));
        let mut lines = report.lines();
        let table = lines
            .next()
            .unwrap()
            .split_once("hodge=[[")
            .unwrap()
            .1
            .strip_suffix("]]")
            .unwrap();
        let mut previous = 0;
        for (p, row) in table.split("], [").enumerate() {
            let diagonal = row.split(", ").nth(p).unwrap().parse::<u128>().unwrap();
            assert_eq!(number(lines.next().unwrap()), diagonal, "n={n}, p={p}");
            let primitive = if p <= n { diagonal - previous } else { 0 };
            assert_eq!(number(lines.next().unwrap()), primitive, "n={n}, p={p}");
            previous = diagonal;
        }
        assert_eq!(lines.next(), None);
    }
}

#[test]
fn multiple_unitary_factors_match_the_direct_torus_product() {
    let cap = 10;
    let chars = torus(cap, &[(1, 0), (-1, 0), (0, 1), (0, -1)]);
    let b: Vec<_> = chars
        .into_iter()
        .map(|row| row[cap][cap] - row[cap + 1][cap - 1])
        .collect();
    let mut program = String::new();
    let mut expected = Vec::new();
    for e in [2, 5] {
        let rho = 23 - 4 * e;
        let mut fixed = vec![(0, 0); rho - 1];
        fixed.extend([(1, 0), (-1, 0)]);
        let mut graded: Vec<Vec<_>> = torus(cap, &fixed)
            .into_iter()
            .map(|plane| plane.into_iter().map(|row| row[cap]).collect())
            .collect();
        for _ in 0..e {
            graded = (0..=cap)
                .map(|n| {
                    (0..=2 * cap)
                        .map(|z| (0..=n).map(|j| b[j] * graded[n - j][z]).sum())
                        .collect()
                })
                .collect();
        }
        for (n, row) in graded.iter().enumerate().skip(2) {
            for p in 0..=2 * n {
                let z = cap - n + p;
                let total = row[z];
                let primitive = if p == 0 {
                    total
                } else if p <= n {
                    total - row[z - 1]
                } else {
                    0
                };
                writeln!(
                    program,
                    "rvHilbComplexMultiplicationHodge[{n},{},{rho},{e}];",
                    2 * p
                )
                .unwrap();
                writeln!(
                    program,
                    "rvHilbComplexMultiplicationPrimitiveHodge[{n},{},{rho},{e}];",
                    2 * p
                )
                .unwrap();
                expected.extend([total, primitive]);
            }
        }
    }
    let got: Vec<_> = stdout(interpret(&program)).lines().map(number).collect();
    assert_eq!(got, expected);
}

#[test]
fn determinant_and_degree_boundaries_preserve_sparse_large_inputs() {
    let report = stdout(interpret(
        r"
        rvHilbComplexMultiplicationHodge[0,0,21,1];
        rvHilbComplexMultiplicationHodge[0,2,21,1];
        rvHilbComplexMultiplicationHodge[2,3,21,1];
        rvHilbComplexMultiplicationPrimitiveHodge[2,6,21,1];
        rvHilbComplexMultiplicationHodge[1000000000000000000,2,1,11];
        rvHilbGLOscillator[1000000000,500000000499999999,1];
        rvHilbGLOscillator[1000000000,500000000500000000,1];
        rvHilbGLOscillator[1000000000,500000000500000001,-1];
        rvHilbGLOscillator[1000000000,500000000500000002,1];
        rvHilbGLOscillator[1,170141183460469231731687303715884105729,-170141183460469231731687303715884105728];
        rvHilbGLOscillator[340282366920938463463374607431768211455,340282366920938463463374607431768211455,1];
        ",
    ));
    assert_eq!(
        report.lines().map(number).collect::<Vec<_>>(),
        [1, 0, 0, 0, 1, 0, 1, 1, 3, 1, 0]
    );

    for (input, expected) in [
        (
            "rvHilbComplexMultiplicationHodge[2,4,22,1];",
            "rank must be in 1..=21",
        ),
        (
            "rvHilbComplexMultiplicationHodge[2,4,2,1];",
            "CM field degree 2e",
        ),
        (
            "rvHilbComplexMultiplicationHodge[2,4,21,0];",
            "CM field degree 2e",
        ),
        (
            "rvHilbComplexMultiplicationHodge[2,4,21,340282366920938463463374607431768211455];",
            "CM field degree 2e",
        ),
        (
            "rvHilbComplexMultiplicationHodge[1,2,21,1];",
            "whose LLV type is D12",
        ),
        ("rvHilbGLOscillator[0,0,0];", "dimension must be positive"),
        (
            "rvHilbGLOscillator[2,1000000000000000000,0];",
            "polynomial computation limit",
        ),
    ] {
        assert!(diagnostic(interpret(input)).contains(expected), "{input}");
    }
}
