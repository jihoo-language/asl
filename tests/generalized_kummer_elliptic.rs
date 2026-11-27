mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

type Row = BTreeMap<i32, i128>;

#[test]
fn worked_example_matches_the_k3_and_composite_lengths() {
    assert_eq!(
        stdout(example("generalized-kummer-elliptic")),
        concat!(
            "rvGenKummerElliptic[1, 0, 0] = 1\n",
            "rvGenKummerElliptic[2, 2, 0] = 1616\n",
            "rvK3EllipticGenus[1, 2, 0] = 1616\n",
            "rvGenKummerElliptic[3, 0, 0] = 90\n",
            "rvGenKummerElliptic[3, 0, 1] = 6\n",
            "rvGenKummerElliptic[3, 0, 2] = 3\n",
            "rvGenKummerElliptic[3, 1, 0] = 4848\n",
            "rvGenKummerElliptic[3, 2, 4] = 90\n",
            "rvGenKummerElliptic[4, 0, 0] = 336\n",
            "rvGenKummerElliptic[4, 2, 5] = 44\n",
        )
    );
}

#[test]
fn entire_length_two_rows_agree_with_the_k3_dmv_v_product() {
    let mut source = String::new();
    for d in 0..=2 {
        for r in -5..=5 {
            writeln!(
                source,
                "rvGenKummerElliptic[2,{d},{r}]; rvK3EllipticGenus[1,{d},{r}];"
            )
            .unwrap();
        }
    }
    let values = coefficients(&stdout(interpret(&source)));
    assert_eq!(values.len(), 66);
    for pair in values.chunks_exact(2) {
        assert_eq!(pair[0], pair[1]);
    }
}

#[test]
fn constant_rows_and_euler_numbers_match_the_divisor_formula() {
    let mut source = String::new();
    for n in [1i32, 3, 4, 6, 12] {
        for r in -n..=n {
            writeln!(source, "rvGenKummerElliptic[{n},0,{r}];").unwrap();
        }
    }
    let values = coefficients(&stdout(interpret(&source)));
    let mut entries = values.into_iter();
    for n in [1i32, 3, 4, 6, 12] {
        let divisors: Vec<_> = (1..=n).filter(|d| n % d == 0).collect();
        let mut euler = 0;
        for r in -n..=n {
            let expected: i128 = divisors
                .iter()
                .map(|&a| {
                    // Coefficient of y^(1-a)(1+y+...+y^(a-1))^2.
                    i128::from((a - r.abs()).max(0)) * i128::from(n) * i128::from(n / a).pow(3)
                })
                .sum();
            let reported = entries.next().unwrap();
            assert_eq!(reported, expected, "length {n}, exponent {r}");
            euler += reported;
        }
        assert_eq!(
            euler,
            i128::from(n).pow(3) * i128::from(divisors.into_iter().sum::<i32>())
        );
    }
    assert!(entries.next().is_none());
}

#[test]
fn mixed_rows_multiply_back_to_the_independent_integral_hecke_numerator() {
    let mut source = String::new();
    for n in [3, 4] {
        for d in 0..=3 {
            for r in -8..=8 {
                writeln!(source, "rvGenKummerElliptic[{n},{d},{r}];").unwrap();
            }
        }
    }
    let values = coefficients(&stdout(interpret(&source)));
    let seed = literal_seed(12);
    for (n, reported) in [3, 4].into_iter().zip(values.chunks_exact(4 * 17)) {
        let genus: Vec<Row> = reported
            .chunks_exact(17)
            .map(|row| {
                (-8..=8)
                    .zip(row.iter().copied())
                    .filter(|(_, value)| *value != 0)
                    .collect()
            })
            .collect();
        for d in 0..=3 {
            let mut convolution = Row::new();
            for j in 0..=d {
                for (&r, &a) in &seed[j] {
                    for (&s, &e) in &genus[d - j] {
                        *convolution.entry(r + s).or_default() += a * e;
                    }
                }
            }
            let mut hecke = Row::new();
            for a in 1..=n {
                if n % a == 0 && d % a == 0 {
                    for (&r, &value) in &seed[n * d / (a * a)] {
                        *hecke.entry(r * a as i32).or_default() +=
                            (n * (n / a).pow(3)) as i128 * value;
                    }
                }
            }
            convolution.retain(|_, value| *value != 0);
            hecke.retain(|_, value| *value != 0);
            assert_eq!(convolution, hecke, "length {n}, q degree {d}");
            for (&r, &value) in &genus[d] {
                assert_eq!(genus[d].get(&(-r)), Some(&value));
            }
            if d != 0 {
                assert_eq!(genus[d].values().sum::<i128>(), 0);
            }
        }
        // Elliptic translation of the q^0 y^-1 coefficient by one period.
        let index = n - 1;
        assert_eq!(
            genus[index - 1].get(&(2 * index as i32 - 1)),
            genus[0].get(&-1)
        );
    }
}

fn literal_seed(top: usize) -> Vec<Row> {
    let mut rows = vec![Row::new(); top + 1];
    rows[0] = Row::from([(-1, 1), (0, -2), (1, 1)]);
    for j in 1..=top {
        for y in [-1, 1] {
            let old = rows.clone();
            for d in j..=top {
                for (&r, &value) in &old[d - j] {
                    *rows[d].entry(r + y).or_default() -= 2 * value;
                }
                if d >= 2 * j {
                    for (&r, &value) in &old[d - 2 * j] {
                        *rows[d].entry(r + 2 * y).or_default() += value;
                    }
                }
            }
        }
        // One unbounded factor per copy gives (1-q^j)^-4.
        for _ in 0..4 {
            for d in j..=top {
                let previous = rows[d - j].clone();
                for (r, value) in previous {
                    *rows[d].entry(r).or_default() += value;
                }
            }
        }
    }
    rows
}

#[test]
fn handles_zero_length_unbounded_input_and_outside_support() {
    for (source, expected) in [
        (
            "rvGenKummerElliptic[0,0,0];",
            "generalized Kummer length must be positive",
        ),
        (
            "rvGenKummerElliptic[340282366920938463463374607431768211455,2,0];",
            "coefficient storage limit",
        ),
        (
            "rvGenKummerElliptic[4,340282366920938463463374607431768211455,0];",
            "coefficient storage limit",
        ),
        (
            "rvGenKummerElliptic[4,5000,0];",
            "coefficient storage limit",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{error}");
    }
    assert_eq!(
        coefficients(&stdout(interpret(
            "rvGenKummerElliptic[1,340282366920938463463374607431768211455,0]; rvGenKummerElliptic[4,2,-1000]; rvGenKummerElliptic[4,0,-4];"
        ))),
        vec![0, 0, 0],
    );
}

fn coefficients(output: &str) -> Vec<i128> {
    output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}
