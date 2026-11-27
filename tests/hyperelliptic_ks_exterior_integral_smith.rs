mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-ks-exterior-integral-smith.asl")
        .split_once("// Integral divided")
        .unwrap()
        .0
}

fn pairs(text: &str) -> Vec<(u128, u128)> {
    let list = text.strip_prefix('[').unwrap().strip_suffix(']').unwrap();
    if list.is_empty() {
        return Vec::new();
    }
    list.split("), (")
        .map(|pair| {
            let (value, count) = pair.trim_matches(['(', ')']).split_once(", ").unwrap();
            (value.parse().unwrap(), count.parse().unwrap())
        })
        .collect()
}

fn report(line: &str) -> (Vec<u128>, Vec<(u128, u128)>) {
    let (_, body) = line.split_once(" = dimension=").unwrap();
    let (heading, factors) = body.split_once(", factors=").unwrap();
    let fields: Vec<_> = heading.split(", ").collect();
    let numeric = fields
        .iter()
        .map(|part| part.rsplit('=').next().unwrap().parse::<u128>().unwrap())
        .collect();
    (numeric, pairs(factors))
}

#[test]
fn example_exercises_wide_scale_middle_torsion_and_nonuniform_factors() {
    let output = stdout(example("hyperelliptic-ks-exterior-integral-smith"));
    let mut lines = output.lines();
    assert_eq!(
        report(lines.next().unwrap()),
        (vec![2016, 120, 1896, 1896, 3, 2], vec![(1, 120)])
    );
    let (middle, factors) = report(lines.next().unwrap());
    assert_eq!(
        middle,
        [
            1832624140942590534,
            608990588904281634,
            1223633552038308900,
            1223633552038308900,
            1,
            4
        ]
    );
    assert_eq!(factors.last(), Some(&(3063060, 1026)));
    assert_eq!(factors.iter().map(|p| p.1).sum::<u128>(), middle[1]);
    assert!(factors
        .windows(2)
        .all(|pair| pair[1].0.is_multiple_of(pair[0].0)));
    assert_eq!(
        report(lines.next().unwrap()).1,
        [(1, 120), (2, 256), (4, 120)]
    );
    let (huge, factors) = report(lines.next().unwrap());
    assert_eq!(&huge[4..], [1u128 << 127, 16]);
    assert_eq!(factors, [(1, 601080390)]);
    assert!(lines.next().is_none());
}

#[test]
fn every_degree_and_order_matches_the_rational_jordan_ranks() {
    let mut source = setup().to_owned();
    for name in ["A", "B3"] {
        for degree in 0..=64 {
            writeln!(source, "rvKSExteriorData[{name},{degree}];").unwrap();
            for order in 0..=33 {
                writeln!(source, "rvKSExteriorSmith[{name},{degree},{order}];").unwrap();
            }
        }
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    for _ in 0..130 {
        let data = lines.next().unwrap();
        let dimension: u128 = data
            .split_once("dimension=")
            .unwrap()
            .1
            .split(',')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        let jordan = pairs(data.split_once("jordan=").unwrap().1);
        for order in 0u128..=33 {
            let rank: u128 = jordan
                .iter()
                .map(|&(size, count)| size.saturating_sub(order) * count)
                .sum();
            let (header, factors) = report(lines.next().unwrap());
            assert_eq!(
                header[..4],
                [dimension, rank, dimension - rank, dimension - rank]
            );
            assert_eq!(factors.iter().map(|p| p.1).sum::<u128>(), rank);
            assert!(factors.windows(2).all(|p| p[1].0.is_multiple_of(p[0].0)));
        }
    }
    assert!(lines.next().is_none());
}

type Column = BTreeMap<u64, i64>;

fn sparse_column(line: &str) -> Column {
    let (_, terms) = line.split_once("{[").unwrap();
    let terms = terms.strip_suffix("] / 1}").unwrap();
    if terms.is_empty() {
        return Column::new();
    }
    terms
        .split("), (")
        .map(|term| {
            let (mask, coefficient) = term.trim_matches(['(', ')']).split_once(", ").unwrap();
            (mask.parse().unwrap(), coefficient.parse().unwrap())
        })
        .collect()
}

fn modular_rank(output: &str, prime: i64) -> usize {
    let mut pivots: BTreeMap<u64, Column> = BTreeMap::new();
    for mut column in output.lines().map(sparse_column) {
        while let Some((&position, &coefficient)) = column.first_key_value() {
            let coefficient = coefficient.rem_euclid(prime);
            if coefficient == 0 {
                column.remove(&position);
            } else if let Some(prior) = pivots.get(&position) {
                let pivot = prior[&position];
                for value in column.values_mut() {
                    *value = (*value * pivot).rem_euclid(prime);
                }
                for (&row, &entry) in prior {
                    let value = column.entry(row).or_default();
                    *value = (*value - coefficient * entry).rem_euclid(prime);
                }
                column.retain(|_, value| *value != 0);
            } else {
                pivots.insert(position, column);
                break;
            }
        }
    }
    pivots.len()
}

#[test]
fn nonuniform_degree_two_spectra_match_full_pbw_matrices_modulo_primes() {
    let output = stdout(interpret(&format!(
        r"{}
        rvKSExteriorSmith[Mixed,1,1]; rvKSExteriorSmith[Mixed,2,1];
        rvKSExteriorSmith[Mixed,2,2]; rvKSExteriorSmith[Mixed,62,2];
        rvKSExteriorSmith[Mixed,63,1];
        ",
        setup()
    )));
    let reports: Vec<_> = output.lines().map(report).collect();
    assert_eq!(reports[0].1, [(1, 16), (2, 16)]);
    assert_eq!(reports[1].1, [(1, 752), (2, 240)]);
    assert_eq!(reports[2].1, [(1, 120), (2, 256), (4, 120)]);
    assert_eq!(reports[3], reports[2]);
    assert_eq!(reports[4], reports[0]);
    for (order, expectation) in [(1, &reports[1]), (2, &reports[2])] {
        let mut source = setup().to_owned();
        for i in 0..64 {
            for j in i + 1..64 {
                writeln!(
                    source,
                    "show rvKSExteriorDivided[Mixed,{order}](rvKSExteriorBasis[P,{}]);",
                    (1u64 << i) | (1u64 << j)
                )
                .unwrap();
            }
        }
        let actual = stdout(interpret(&source));
        assert_eq!(actual.lines().count(), 2016);
        for prime in [2, 3, 5] {
            let rank = expectation
                .1
                .iter()
                .filter(|p| p.0 % prime as u128 != 0)
                .map(|p| p.1)
                .sum::<u128>();
            assert_eq!(
                modular_rank(&actual, prime) as u128,
                rank,
                "order {order}, prime {prime}"
            );
        }
    }
}

fn blade(active: &[usize], sources: u64) -> u64 {
    active
        .iter()
        .enumerate()
        .map(|(i, &bit)| 1u64 << (bit + usize::from(sources & (1 << i) == 0)))
        .sum()
}

#[test]
fn explicit_pbw_boolean_block_has_integral_torsion_and_divided_powers() {
    // These are three actual adjacent PBW pairs of the integral spin action.
    let active = [0, 6, 10];
    let mut source = setup().to_owned();
    for index in active {
        writeln!(source, "assert rvKSExteriorDivided[A,1](rvKSExteriorBasis[P,{}])==rvKSExteriorScale[3](rvKSExteriorBasis[P,{}]);", 1u64 << index, 1u64 << (index + 1)).unwrap();
    }
    for column in [3, 5, 6] {
        writeln!(
            source,
            "show rvKSExteriorDivided[A,1](rvKSExteriorBasis[P,{}]);",
            blade(&active, column)
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let columns: Vec<_> = output.lines().map(sparse_column).collect();
    let targets = [1, 2, 4].map(|mask| blade(&active, mask));
    let matrix: Vec<Vec<_>> = columns
        .iter()
        .map(|column| {
            assert_eq!(column.len(), 2);
            assert!(column.values().all(|&c| c == 3));
            targets
                .iter()
                .map(|target| column.get(target).copied().unwrap_or(0) / 3)
                .collect()
        })
        .collect();
    assert_eq!(matrix, [vec![1, 1, 0], vec![1, 0, 1], vec![0, 1, 1]]);
    // A unit entry and minor, and determinant ±2, give primitive factors 1,1,2.
    let [a, b, c] = &matrix[..] else {
        unreachable!()
    };
    assert_eq!(a[0], 1);
    assert_eq!(a[0] * b[1] - a[1] * b[0], -1);
    assert_eq!(
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]),
        -2
    );
    let actual = stdout(interpret(&format!("{}rvKSExteriorSmith[A,3,1];", setup())));
    assert_eq!(report(actual.trim()).1, [(1, 18_096), (2, 560)]);

    let all_active = [0, 6, 10, 12, 18, 20, 24, 30, 34, 36, 40, 46, 48, 54, 58, 60];
    let sources = blade(&all_active, (1 << 16) - 1);
    let paired = sources | (sources << 1);
    let inert: u64 = (0..64)
        .map(|i| 1u64 << i)
        .filter(|bit| bit & paired == 0)
        .take(16)
        .sum();
    let source = format!(
        r"{}
        let x=rvKSExteriorBasis[P,{}];
        let top=rvKSExteriorBasis[P,18446744073709551615];
        assert rvKSExteriorDivided[A,0](x)==x;
        assert rvKSExteriorDivided[A,3](x)==rvKSExteriorScale[27](rvKSExteriorBasis[P,{}]);
        assert rvKSExteriorPower[A,3](x)==rvKSExteriorScale[6](rvKSExteriorDivided[A,3](x));
        assert rvKSExteriorPower[A,2](x)==rvKSExteriorScale[2](rvKSExteriorDivided[A,2](x));
        assert rvKSExteriorDivided[A,4](x)==zero[RVKSExterior(P,3)];
        assert rvKSExteriorDivided[A,1](top)==zero[RVKSExterior(P,64)];
        let dual=rvKSExteriorBasis[P,{}];
        assert rvKSExteriorPower[A,3](dual)==rvKSExteriorScale[6](rvKSExteriorDivided[A,3](dual));
        assert rvKSExteriorDivided[A,16](rvKSExteriorBasis[P,{}])==rvKSExteriorScale[{}](rvKSExteriorBasis[P,{}]);
        assert rvKSExteriorDivided[A,17](rvKSExteriorBasis[P,{}])==zero[RVKSExterior(P,32)];
        ",
        setup(),
        blade(&active, 7),
        blade(&active, 0),
        !blade(&active, 0),
        sources,
        3u64.pow(16),
        blade(&all_active, 0),
        sources | inert,
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn boundary_powers_and_nonuniform_scope_have_specific_diagnostics() {
    let out = stdout(interpret(&format!("{}rvKSExteriorSmith[Mixed,32,0];rvKSExteriorSmith[Mixed,32,33];rvKSExteriorSmith[Mixed,64,1];", setup())));
    let mut reports = out.lines().map(report);
    let (middle, factors) = reports.next().unwrap();
    assert_eq!(middle[0], middle[1]);
    assert_eq!(factors, [(1, middle[0])]);
    assert!(reports.next().unwrap().1.is_empty());
    assert_eq!(reports.next().unwrap().0, [1, 0, 1, 1, 1, 0]);
    for (source, expected) in [
        (
            "rvKSExteriorSmith[Mixed,32,1];",
            "nonuniform Kuga–Satake spin Smith factors",
        ),
        (
            "rvKSExteriorSmith[Mixed,65,0];",
            "cohomological degree must be between 0 and 64",
        ),
        (
            "rvog6cusp Q(h,e,f);show rvKSExteriorDivided[A,2](rvKSExteriorBasis[Q,3]);",
            "expected type RVKSExterior(P,2)",
        ),
    ] {
        assert!(diagnostic(interpret(&format!("{}{source}", setup()))).contains(expected));
    }
}
