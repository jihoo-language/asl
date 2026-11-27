mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/hyperelliptic-witt-norm.asl");
const HUGE: u128 = u128::MAX;

type Matrix = Vec<Vec<i128>>;

fn geometry() -> &'static str {
    let start = EXAMPLE.find("rvunram F").unwrap();
    let end = EXAMPLE.find("let extensionSpectrum").unwrap();
    &EXAMPLE[start..end]
}

fn numbers(output: &str) -> Vec<i128> {
    output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn example_matches_the_trace_of_the_actual_weil_restriction() {
    assert_eq!(
        stdout(example("hyperelliptic-witt-norm")),
        concat!(
            "RVWittSpectrum{[1, 1], [1]}\n",
            "rvWittSpectrumNormGhost[2, 1] = -1\n",
            "rvWittSpectrumNormGhost[2, 2] = 1\n",
            "rvWittSpectrumNormGhost[2, 1] = 248\n",
            "rvKummerHilbWeilRestrictPoints[C, 1, 2, 1] = 248\n",
            "rvWittSpectrumNormGhost[2, 2] = 61504\n",
            "rvKummerHilbWeilRestrictPoints[C, 1, 2, 2] = 61504\n",
        )
    );
}

fn cyclic_tensor_matrices(f: &[[i128; 3]; 3], odd: [bool; 3], n: u32) -> [Matrix; 2] {
    let words: Vec<Vec<usize>> = (0..3usize.pow(n))
        .map(|mut index| {
            (0..n)
                .map(|_| {
                    let digit = index % 3;
                    index /= 3;
                    digit
                })
                .collect()
        })
        .collect();
    let mut ranks = [0; 2];
    let positions: Vec<_> = words
        .iter()
        .map(|word| {
            let parity = word.iter().filter(|&&i| odd[i]).count() % 2;
            ranks[parity] += 1;
            (parity, ranks[parity] - 1)
        })
        .collect();
    let mut result = ranks.map(|rank| vec![vec![0; rank]; rank]);
    for ((parity, col), word) in positions.iter().copied().zip(words) {
        let first = word[0];
        let other_odd = word[1..].iter().filter(|&&i| odd[i]).count();
        let sign = if odd[first] && !other_odd.is_multiple_of(2) {
            -1
        } else {
            1
        };
        for (last, row) in f.iter().enumerate() {
            let coefficient = row[first];
            if coefficient == 0 {
                continue;
            }
            let index = word[1..]
                .iter()
                .copied()
                .chain([last])
                .rev()
                .fold(0, |number, letter| 3 * number + letter);
            let (target_parity, row) = positions[index];
            assert_eq!(target_parity, parity);
            result[parity][row][col] += sign * coefficient;
        }
    }
    result
}

fn traces(matrix: &Matrix, through: usize) -> Vec<i128> {
    let n = matrix.len();
    let mut power: Matrix = (0..n)
        .map(|i| (0..n).map(|j| i128::from(i == j)).collect())
        .collect();
    (0..through)
        .map(|_| {
            power = power
                .iter()
                .map(|row| {
                    (0..n)
                        .map(|j| row.iter().zip(matrix).map(|(&a, b)| a * b[j]).sum())
                        .collect()
                })
                .collect();
            power.iter().enumerate().map(|(i, row)| row[i]).sum()
        })
        .collect()
}

fn determinant(power_traces: &[i128], rank: usize) -> Vec<i128> {
    let mut coefficients = vec![1];
    for k in 1..=rank {
        let numerator: i128 = coefficients
            .iter()
            .rev()
            .zip(power_traces)
            .map(|(a, b)| a * b)
            .sum();
        assert_eq!(numerator % k as i128, 0);
        coefficients.push(-numerator / k as i128);
    }
    coefficients
}

#[test]
fn full_norm_and_ghosts_match_signed_cyclic_tensor_matrices() {
    for (input, f, odd, n) in [
        (
            "rvWittSpectrum{[1,-1,1],[1,1]}",
            [[0, -1, 0], [1, 1, 0], [0, 0, -1]],
            [false, false, true],
            2,
        ),
        (
            "rvWittSpectrum{[1,-2],[1,-1,1]}",
            [[2, 0, 0], [0, 0, -1], [0, 1, 1]],
            [false, true, true],
            3,
        ),
        (
            "rvWittSpectrum{[1],[1,-3,3,-2]}",
            [[2, 0, 0], [0, 0, -1], [0, 1, 1]],
            [true, true, true],
            2,
        ),
        (
            "rvWittSpectrum{[1],[1,-3,3,-2]}",
            [[2, 0, 0], [0, 0, -1], [0, 1, 1]],
            [true, true, true],
            3,
        ),
    ] {
        let [a, b] = cyclic_tensor_matrices(&f, odd, n);
        let last = 9;
        let pos = traces(&a, last.max(a.len()));
        let neg = traces(&b, last.max(b.len()));
        let p = determinant(&pos, a.len());
        let q = determinant(&neg, b.len());
        let mut source = format!("let x = {input}; let nx = rvWittSpectrumNorm[{n}](x); assert nx == rvWittSpectrum{{{p:?},{q:?}}};");
        for k in 1..=last {
            writeln!(
                source,
                "rvWittSpectrumGhost[{k}](nx); rvWittSpectrumNormGhost[{n},{k}](x);"
            )
            .unwrap();
        }
        let expected: Vec<_> = pos
            .iter()
            .zip(&neg)
            .take(last)
            .flat_map(|(a, b)| [a - b, a - b])
            .collect();
        assert_eq!(numbers(&stdout(interpret(&source))), expected);
    }
}

#[test]
fn transitivity_multiplicativity_and_quadratic_distributivity() {
    let output = stdout(interpret(
        r"
        let unit = rvWittSpectrumLine[1];
        let x = add(rvWittSpectrumLine[2],neg(unit));
        let y = add(unit,neg(rvWittSpectrumLine[-1]));
        let n2x = rvWittSpectrumNorm[2](x);
        let n3x = rvWittSpectrumNorm[3](x);
        assert rvWittSpectrumNorm[1](x) == x;
        assert rvWittSpectrumNorm[3](n2x) == rvWittSpectrumNorm[6](x);
        assert rvWittSpectrumNorm[2](n3x) == rvWittSpectrumNorm[6](x);
        assert rvWittSpectrumNorm[2](rvWittSpectrumTensor(x,y)) == rvWittSpectrumTensor(n2x,rvWittSpectrumNorm[2](y));
        assert rvWittSpectrumAdams[3](n3x) == rvWittSpectrumTensor(rvWittSpectrumTensor(x,x),x);
        assert rvWittSpectrumNorm[2](add(x,y)) == add(add(n2x,rvWittSpectrumNorm[2](y)),rvWittSpectrumVerschiebung[2](rvWittSpectrumTensor(x,y)));
        assert rvWittSpectrumNorm[2](neg(rvWittSpectrumLine[-3])) == rvWittSpectrumLine[3];
        assert rvWittSpectrumNorm[3](neg(rvWittSpectrumLine[-3])) == neg(rvWittSpectrumLine[-3]);
    ",
    ));
    assert!(output.is_empty());
}

#[test]
fn trace_only_queries_do_not_build_large_determinants() {
    let output = stdout(interpret(&format!(
        r"
        let empty = zero[RVWittSpectrum];
        let unit = rvWittSpectrumLine[1];
        let difference = add(rvWittSpectrumLine[2],neg(unit));
        assert rvWittSpectrumNorm[{HUGE}](empty) == empty;
        assert rvWittSpectrumNorm[{HUGE}](unit) == unit;
        assert rvWittSpectrumNorm[{HUGE}](neg(unit)) == neg(unit);
        rvWittSpectrumNormGhost[{HUGE},{HUGE}](difference);
        rvWittSpectrumNormGhost[{HUGE},5](difference);
        rvWittSpectrumNormGhost[{HUGE},{HUGE}](empty);
        rvWittSpectrumNormGhost[100,7](rvWittSpectrumLine[2]);
        rvWittSpectrumNormGhost[512,512](neg(unit));
    "
    )));
    assert_eq!(numbers(&output), [1, 1, 0, 128, 1]);
    let source = format!(
        "{} show rvWittSpectrumNorm[2](rvK3FrobHilbSpectrum[Idle,C,1,2]);",
        geometry()
    );
    assert!(diagnostic(interpret(&source)).contains("determinant degree limit"));
}

#[test]
fn weil_restriction_uses_the_lcm_field_and_the_gcd_number_of_factors() {
    let output = stdout(interpret(&format!(
        r"
        {}
        rvunram E[3,1]{{1,2,0,1}};
        rvjacbase U[C,E];
        card (\x:RVKummerK3(U)=>S.hi);
        rvKummerHilbWeilRestrictPoints[C,1,3,1];
        rvKummerHilbWeilRestrictPoints[C,1,3,3];
        rvKummerHilbPoints[C,2,6];
        rvKummerHilbWeilRestrictPoints[C,2,2,3];
        rvKummerHilbWeilRestrictPoints[C,2,3,2];
        rvKummerHilbPoints[C,2,2];
        rvKummerHilbWeilRestrictPoints[C,2,2,2];
        rvKummerHilbWeilRestrictPoints[C,0,{HUGE},340282366920938463463374607431768211454];
    ",
        geometry()
    )));
    let (points, counts) = output.split_once('\n').unwrap();
    assert_eq!(points, "card = 992 / 992");
    let counts = numbers(counts);
    assert_eq!(&counts[..2], [992, 992i128.pow(3)]);
    assert_eq!(counts[2..5], [counts[2]; 3]);
    assert_eq!(counts[6], counts[5].pow(2));
    assert_eq!(counts[7], 1);
}

#[test]
fn invalid_degrees_and_arithmetic_limits_have_direct_diagnostics() {
    for (input, expected) in [
        (
            "show rvWittSpectrumNorm[0](rvWittSpectrumLine[1]);",
            "norm degree must be positive",
        ),
        ("show rvWittSpectrumNorm[2](S.hi);", "rvWittSpectrumNorm"),
        (
            "rvWittSpectrumNormGhost[0,1](rvWittSpectrumLine[1]);",
            "must be positive",
        ),
        (
            "rvWittSpectrumNormGhost[1,0](rvWittSpectrumLine[1]);",
            "must be positive",
        ),
        ("rvWittSpectrumNormGhost[2,4](S.hi);", "RVWittSpectrum"),
        (
            "show rvWittSpectrumNorm[9](rvWittSpectrum{[1,-2,1]});",
            "determinant degree limit",
        ),
        (
            "rvWittSpectrumNormGhost[127,127](rvWittSpectrum{[1,-2,1]});",
            "signed 128-bit range",
        ),
        (
            "rvWittSpectrumNormGhost[512,512](rvWittSpectrum{[1,-2,1]});",
            "512 bits",
        ),
    ] {
        let error = diagnostic(interpret(input));
        assert!(error.contains(expected), "{error}");
    }
    for (input, expected) in [
        (
            "rvKummerHilbWeilRestrictPoints[C,0,0,1];",
            "must be positive",
        ),
        (
            "rvKummerHilbWeilRestrictPoints[C,0,1,0];",
            "must be positive",
        ),
        ("rvKummerHilbWeilRestrictPoints[Unknown,0,1,1];", "Unknown"),
        (
            "rvKummerHilbWeilRestrictPoints[C,1,340282366920938463463374607431768211455,2];",
            "512 bits",
        ),
        (
            "rvKummerHilbWeilRestrictPoints[C,1,4,24];",
            "unsigned 128-bit range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{} {input}", geometry())));
        assert!(error.contains(expected), "{error}");
    }
}
