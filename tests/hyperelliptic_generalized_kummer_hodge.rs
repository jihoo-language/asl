mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn example_reports_hodge_numbers_and_primitive_eigenvalue_multiplicities() {
    assert_eq!(
        stdout(example("hyperelliptic-generalized-kummer-hodge")),
        concat!(
            "rvK3HilbHodge[2, 2, 2] = 232\n",
            "rvJacHilbHodge[2, 1, 1] = 9\n",
            "rvGenKummerHodge[3, 2, 1] = 4\n",
            "rvGenKummerHodge[3, 2, 2] = 96\n",
            "rvGenKummerTranslateHodge[2, 2, 1, 1] = 4\n",
            "rvGenKummerHodgeSpectrum[2, 2, 1, 1] = [(1, 12), (2, 8)]\n",
            "rvGenKummerTranslateHodge[3, 3, 2, 2] = 15\n",
            "rvGenKummerHodgeSpectrum[3, 3, 2, 2] = [(1, 42), (3, 27)]\n",
        ),
    );
}

#[test]
fn full_hodge_diamonds_refine_all_graded_dimensions_and_characters() {
    for (hodge, graded) in [
        ("rvK3HilbHodge", "rvKummerHilbBetti"),
        ("rvJacHilbHodge", "rvJacHilbBetti"),
    ] {
        for n in 0..=3 {
            check_diamond(&format!("{hodge}[{n}"), &format!("{graded}[{n}"), 2 * n);
        }
    }
    for n in 1..=5 {
        let dimension = 2 * (n - 1);
        check_diamond(
            &format!("rvGenKummerHodge[{n}"),
            &format!("rvGenKummerBetti[{n}"),
            dimension,
        );
        for h in (1..=n).filter(|h| n % h == 0) {
            check_diamond(
                &format!("rvGenKummerTranslateHodge[{n},{h}"),
                &format!("rvGenKummerTranslateChar[{n},{h}"),
                dimension,
            );
        }
    }
}

#[test]
fn spectra_reconstruct_the_character_of_every_power_in_low_lengths() {
    for (n, h) in [(2, 2), (3, 3), (4, 2), (4, 4), (6, 6)] {
        let dimension = 2 * (n - 1);
        let mut source = String::new();
        for p in 0..=dimension {
            for q in 0..=dimension {
                writeln!(source, "rvGenKummerHodgeSpectrum[{n},{h},{p},{q}];").unwrap();
                for k in 1..=h {
                    let order = h / gcd(h, k);
                    writeln!(source, "rvGenKummerTranslateHodge[{n},{order},{p},{q}];").unwrap();
                }
                writeln!(source, "rvGenKummerHodge[{n},{p},{q}];").unwrap();
            }
        }
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        for p in 0..=dimension {
            for q in 0..=dimension {
                let row = lines.next().unwrap();
                let numbers: Vec<i128> = rhs(row)
                    .split(|c: char| !c.is_ascii_digit())
                    .filter(|component| !component.is_empty())
                    .map(|component| component.parse().unwrap())
                    .collect();
                assert!(numbers.len().is_multiple_of(2));
                let spectrum: Vec<_> = numbers
                    .chunks_exact(2)
                    .map(|pair| (pair[0] as usize, pair[1]))
                    .collect();
                for k in 1..=h {
                    let trace = spectrum
                        .iter()
                        .map(|&(d, count)| count * small_root_trace(d, k))
                        .sum::<i128>();
                    assert_eq!(trace, number(lines.next().unwrap()), "{row}, power {k}");
                }
                let rank = spectrum
                    .iter()
                    .map(|&(d, count)| count * small_root_trace(d, h))
                    .sum::<i128>();
                assert_eq!(rank, number(lines.next().unwrap()), "{row}");
                if p == 3 && q == 3 && h == 4 {
                    assert_eq!(spectrum, [(1, 148), (2, 96), (4, 64)]);
                }
            }
        }
        assert!(lines.next().is_none());
    }
}

#[test]
fn signed_hodge_specialization_matches_the_constant_elliptic_genus() {
    for n in 0..=3 {
        let mut source = diamond_source(&format!("rvK3HilbHodge[{n}"), 2 * n);
        for p in 0..=2 * n {
            writeln!(
                source,
                "rvK3EllipticGenus[{n},0,{}];",
                p as isize - n as isize
            )
            .unwrap();
        }
        let output = stdout(interpret(&source));
        let (diamond, rest) = read_diamond(&output, 2 * n);
        assert_eq!(signed_specialization(&diamond), rest);
    }
    for (n, expected) in [
        (1, vec![1]),
        (2, vec![2, 20, 2]),
        (3, vec![3, 6, 90, 6, 3]),
        (4, vec![4, 8, 44, 336, 44, 8, 4]),
    ] {
        let dimension = 2 * (n - 1);
        let source = diamond_source(&format!("rvGenKummerHodge[{n}"), dimension);
        let output = stdout(interpret(&source));
        let (diamond, _) = read_diamond(&output, dimension);
        assert_eq!(signed_specialization(&diamond), expected);
        for q in 0..=dimension {
            assert_eq!(diamond[0][q], i128::from(q.is_multiple_of(2)));
        }
    }
    for n in 1..=3 {
        let output = stdout(interpret(&diamond_source(
            &format!("rvJacHilbHodge[{n}"),
            2 * n,
        )));
        let (diamond, _) = read_diamond(&output, 2 * n);
        assert_eq!(signed_specialization(&diamond), vec![0; 2 * n + 1]);
    }
}

#[test]
fn stable_and_dual_bidegrees_handle_wide_inputs_with_checked_limits() {
    let max = u128::MAX;
    let length = max / 2 + 2;
    let output = stdout(interpret(&format!(
        "rvK3HilbHodge[0,0,0];
         rvK3HilbHodge[0,0,1];
         rvK3HilbHodge[{max},0,{max}];
         rvJacHilbHodge[{max},0,{max}];
         rvGenKummerHodge[{length},{max},{max}];
         rvGenKummerHodge[{}, {max},1];
         rvGenKummerHodgeSpectrum[{max},{max},1,1];
         rvGenKummerHodgeSpectrum[{max},{max},1,0];
         rvGenKummerHodgeSpectrum[1,1,0,0];",
        length - 1,
    )));
    let values: Vec<_> = output.lines().map(rhs).collect();
    assert_eq!(
        values,
        ["1", "0", "0", "2", "5", "0", "[(1, 5)]", "[]", "[(1, 1)]"],
    );
    for (statement, expected) in [
        ("rvGenKummerHodge[0,0,0];", "length must be positive"),
        ("rvGenKummerTranslateHodge[4,0,1,1];", "positive divisor"),
        ("rvGenKummerHodgeSpectrum[4,3,10,10];", "positive divisor"),
        (
            "rvJacHilbHodge[2,-1,1];",
            "expected a nonnegative decimal integer",
        ),
        (
            "rvK3HilbHodge[10000,10000,10000];",
            "coefficient table exceeds the storage limit",
        ),
        (
            "rvGenKummerHodgeSpectrum[10000,10000,9999,9999];",
            "coefficient table exceeds the storage limit",
        ),
    ] {
        let message = diagnostic(interpret(statement));
        assert!(message.contains(expected), "{message}");
    }
}

fn check_diamond(hodge: &str, graded: &str, dimension: usize) {
    let mut source = diamond_source(hodge, dimension);
    for i in 0..=2 * dimension {
        writeln!(source, "{graded},{i}];").unwrap();
    }
    let output = stdout(interpret(&source));
    let (diamond, expected) = read_diamond(&output, dimension);
    let mut degrees = vec![0; 2 * dimension + 1];
    for p in 0..=dimension {
        for q in 0..=dimension {
            let value = diamond[p][q];
            assert_eq!(value, diamond[q][p]);
            assert_eq!(value, diamond[dimension - p][dimension - q]);
            degrees[p + q] += value;
        }
    }
    assert_eq!(degrees, expected, "{hodge}");
}

fn diamond_source(hodge: &str, dimension: usize) -> String {
    let mut source = String::new();
    for p in 0..=dimension {
        for q in 0..=dimension {
            writeln!(source, "{hodge},{p},{q}];").unwrap();
        }
    }
    source
}

fn read_diamond(output: &str, dimension: usize) -> (Vec<Vec<i128>>, Vec<i128>) {
    let mut lines = output.lines().map(number);
    let matrix = (0..=dimension)
        .map(|_| lines.by_ref().take(dimension + 1).collect())
        .collect();
    (matrix, lines.collect())
}

fn signed_specialization(diamond: &[Vec<i128>]) -> Vec<i128> {
    diamond
        .iter()
        .enumerate()
        .map(|(p, row)| {
            row.iter()
                .enumerate()
                .map(|(q, &value)| if (p + q) % 2 == 0 { value } else { -value })
                .sum()
        })
        .collect()
}

fn rhs(line: &str) -> &str {
    line.split_once(" = ").unwrap().1
}

fn number(line: &str) -> i128 {
    rhs(line).parse().unwrap()
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn small_root_trace(d: usize, k: usize) -> i128 {
    match (d, k % d) {
        (1, _) | (2, 0) => 1,
        (2, _) | (3, 1 | 2) => -1,
        (3 | 4 | 6, 0) => 2,
        (4, 2) | (6, 3) => -2,
        (4, _) => 0,
        (6, 1 | 5) => 1,
        (6, 2 | 4) => -1,
        _ => panic!("unexpected root order {d}"),
    }
}
