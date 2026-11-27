mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-ks-mapping-torus.asl")
        .split_once("// Integral mapping")
        .unwrap()
        .0
}

#[derive(Debug, Eq, PartialEq)]
struct Report<'a> {
    fields: BTreeMap<&'a str, u128>,
    factors: Vec<(u128, u128)>,
}

fn report(line: &str) -> Report<'_> {
    let (_, body) = line.split_once(" = ").unwrap();
    let (heading, factors) = body
        .split_once(", factors=")
        .or_else(|| body.split_once(", torsion="))
        .unwrap();
    let fields = heading
        .split(", ")
        .map(|part| {
            let (name, value) = part.split_once('=').unwrap();
            (name, value.parse().unwrap())
        })
        .collect();
    let factors = factors
        .strip_prefix('[')
        .unwrap()
        .strip_suffix(']')
        .unwrap();
    let factors = if factors.is_empty() {
        Vec::new()
    } else {
        factors
            .split("), (")
            .map(|pair| {
                let (value, count) = pair.trim_matches(['(', ')']).split_once(", ").unwrap();
                (value.parse().unwrap(), count.parse().unwrap())
            })
            .collect()
    };
    Report { fields, factors }
}

#[test]
fn example_shows_different_actions_and_exact_large_cover_torsion() {
    let output = stdout(example("hyperelliptic-ks-mapping-torus"));
    let mut lines = output.lines();
    assert_eq!(lines.next().unwrap(), "rvKSExteriorCoeff[130] = 9");
    assert_eq!(lines.next().unwrap(), "rvKSExteriorCoeff[130] = 0");
    let group = report(lines.next().unwrap());
    assert_eq!(group.fields["rank"], 120);
    assert_eq!(group.fields["cokernel_free"], 1896);
    assert_eq!(group.fields["common_base"], 3);
    assert_eq!(group.fields["cover_base"], 2);
    assert_eq!(group.fields["common_power"], 2);
    assert_eq!(group.fields["factorial"], 2);
    assert_eq!(group.factors, [(1, 120)]);
    let torus = report(lines.next().unwrap());
    assert_eq!(torus.fields["cover_base"], 1u128 << 127);
    assert_eq!(torus.fields["common_base"], 3);
    assert_eq!(torus.fields["free"], 1312);
    assert_eq!(torus.factors, [(1, 16)]);
    let product = report(lines.next().unwrap());
    assert_eq!(product.fields["free"], choose(65, 32));
    assert!(product.factors.is_empty());
    assert!(lines.next().is_none());
}

type Column = BTreeMap<u64, i128>;

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

fn active_blade(active: &[usize], f_sources: usize) -> u64 {
    active
        .iter()
        .enumerate()
        .map(|(i, &bit)| 1u64 << (bit + usize::from(f_sources & (1 << i) == 0)))
        .sum()
}

fn exact_pbw_matrix(name: &str, m: i128, order: usize, active: &[usize]) -> Vec<Vec<i128>> {
    let mut source = setup().to_owned();
    let masks: Vec<_> = (0..1 << active.len())
        .map(|i| active_blade(active, i))
        .collect();
    for &mask in &masks {
        writeln!(
            source,
            "show rvKSExteriorDifference[{name},{m},{order}](rvKSExteriorBasis[P,{mask}]);"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    output
        .lines()
        .map(|line| {
            let column = sparse_column(line);
            assert!(column.keys().all(|mask| masks.contains(mask)));
            masks
                .iter()
                .map(|mask| column.get(mask).copied().unwrap_or(0))
                .collect()
        })
        .collect()
}

fn determinant(columns: &[Vec<i128>], col_mask: usize, row_mask: usize) -> i128 {
    if col_mask == 0 {
        return 1;
    }
    let col = col_mask.trailing_zeros() as usize;
    (0..columns.len())
        .filter(|&row| row_mask & (1 << row) != 0)
        .enumerate()
        .map(|(i, row)| {
            let sign = if i.is_multiple_of(2) { 1 } else { -1 };
            sign * columns[col][row]
                * determinant(columns, col_mask ^ (1 << col), row_mask ^ (1 << row))
        })
        .sum()
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn smith_via_minors(matrix: &[Vec<i128>]) -> Vec<u128> {
    let mut previous = 1;
    let mut factors = Vec::new();
    for size in 1..=matrix.len() {
        let subsets: Vec<_> = (1usize..1 << matrix.len())
            .filter(|s| s.count_ones() as usize == size)
            .collect();
        let mut divisor = 0;
        for &rows in &subsets {
            for &cols in &subsets {
                divisor = gcd(divisor, determinant(matrix, cols, rows).unsigned_abs());
            }
        }
        if divisor == 0 {
            break;
        }
        assert!(divisor.is_multiple_of(previous));
        factors.push(divisor / previous);
        previous = divisor;
    }
    factors
}

#[test]
fn independent_integer_minors_on_actual_pbw_group_matrices() {
    // Type II maps three adjacent PBW pairs with common elementary factor 3.
    for (order, expected) in [
        (1, vec![6, 6, 6, 6, 12]),
        (2, vec![72, 72]),
        (3, vec![1296]),
    ] {
        assert_eq!(
            smith_via_minors(&exact_pbw_matrix("A", -2, order, &[0, 6, 10])),
            expected
        );
    }
    // Here an actual Type III marking has one factor 2 and one factor 1.
    assert_eq!(
        smith_via_minors(&exact_pbw_matrix("Mixed", 3, 1, &[4, 6])),
        [3, 3]
    );
    assert_eq!(
        smith_via_minors(&exact_pbw_matrix("Mixed", 3, 2, &[4, 6])),
        [36]
    );
    let output = stdout(interpret(&format!(
        r"{}
        rvKSExteriorGroupSmith[A,3,-2,1];
        rvKSExteriorGroupSmith[Mixed,2,3,1];
        rvKSExteriorGroupSmith[Mixed,2,3,2];
        rvKSExteriorGroupSmith[A,3,1,1];
        rvKSExteriorSmith[A,3,1];
        rvKSExteriorGroupSmith[A,32,-170141183460469231731687303715884105728,16];
        ",
        setup()
    )));
    let groups: Vec<_> = output.lines().map(report).collect();
    assert_eq!(groups[0].factors, [(1, 18096), (2, 560)]);
    assert_eq!(groups[1].factors, [(1, 752), (2, 240)]);
    assert_eq!(groups[2].factors, [(1, 120), (2, 256), (4, 120)]);
    assert_eq!(groups[2].fields["factorial"], 2);
    assert_eq!(groups[3].factors, groups[4].factors);
    for field in ["rank", "cokernel_free", "common_base", "common_power"] {
        assert_eq!(groups[3].fields[field], groups[4].fields[field]);
    }
    assert_eq!(groups[5].fields["common_base"], 3);
    assert_eq!(groups[5].fields["cover_base"], 1u128 << 127);
    assert_eq!(groups[5].fields["common_power"], 16);
    assert_eq!(groups[5].fields["factorial"], (1..=16).product::<u128>());
    assert_eq!(groups[5].factors, [(1, 601080390)]);
}

#[test]
fn finite_differences_agree_with_repeated_actions_and_cover_signs() {
    let mut source = setup().to_owned();
    let active = [0, 6, 10];
    for (i, mask) in [active_blade(&active, 7), !active_blade(&active, 0)]
        .into_iter()
        .enumerate()
    {
        writeln!(
            source,
            r"
            let x{i}=rvKSExteriorBasis[P,{mask}];
            let a{i}=rvKSExteriorDifference[A,-2,1](x{i});
            assert a{i}==add(rvKSExteriorExp[A,-2](x{i}),neg(x{i}));
            let b{i}=rvKSExteriorDifference[A,-2,2](x{i});
            assert b{i}==rvKSExteriorDifference[A,-2,1](a{i});
            assert rvKSExteriorDifference[A,-2,3](x{i})==rvKSExteriorDifference[A,-2,1](b{i});
            assert rvKSExteriorDifference[A,-2,0](x{i})==x{i};
            assert rvKSExteriorDifference[A,0,0](x{i})==x{i};
            "
        )
        .unwrap();
    }
    let all_active = [0, 6, 10, 12, 18, 20, 24, 30, 34, 36, 40, 46, 48, 54, 58, 60];
    let from = active_blade(&all_active, (1 << 16) - 1);
    let to = active_blade(&all_active, 0);
    let coefficient = 6i128.pow(16) * (1..=16i128).product::<i128>();
    writeln!(source, r"
        let large=rvKSExteriorBasis[P,{from}];
        assert rvKSExteriorDifference[A,-2,16](large)==rvKSExteriorScale[{coefficient}](rvKSExteriorBasis[P,{to}]);
        assert rvKSExteriorDifference[A,0,3](large)==zero[RVKSExterior(P,16)];
        assert rvKSExteriorDifference[A,5,17](large)==zero[RVKSExterior(P,16)];
        ").unwrap();
    assert_eq!(stdout(interpret(&source)), "");
}

fn choose(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    (0..k).fold(1, |r, i| r * (n - i) as u128 / (i + 1) as u128)
}

fn mod2_half_rank(k: usize, rank: usize) -> u128 {
    if k > 64 {
        return 0;
    }
    let mut count = 0;
    for singles in 1..=rank.min(k) {
        for doubles in 0..=(rank - singles).min((k - singles) / 2) {
            let inert = k - singles - 2 * doubles;
            let blocks = choose(rank, singles)
                * choose(rank - singles, doubles)
                * choose(64 - 2 * rank, inert);
            // Over F2 a linear triangular change replaces sum(X_i) by X_1.
            count += blocks * (1u128 << (singles - 1));
        }
    }
    count
}

fn p_torsion(group: Option<&Report<'_>>, p: u128) -> u128 {
    let Some(group) = group else {
        return 0;
    };
    group
        .factors
        .iter()
        .filter(|&(factor, _)| {
            (group.fields["common_base"] * group.fields["cover_base"] * factor).is_multiple_of(p)
        })
        .map(|p| p.1)
        .sum()
}

#[test]
fn poincare_duality_and_independent_universal_coefficients_in_every_degree() {
    let mut source = setup().to_owned();
    for (name, m) in [("A", -5), ("B3", 1), ("Mixed", 0)] {
        for degree in 0..=65 {
            writeln!(source, "rvKSMappingTorus[{name},{m},{degree}];").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let reports: Vec<_> = output.lines().map(report).collect();
    assert_eq!(reports.len(), 198);
    for (groups, rank) in reports.chunks_exact(66).zip([16, 32, 0]) {
        for (degree, group) in groups.iter().enumerate() {
            let free = group.fields["free"];
            assert_eq!(free, groups[65 - degree].fields["free"]);
            assert_eq!(
                free,
                group.fields["fiber_invariants"] + group.fields["fiber_coinvariants_free"]
            );
            if let Some(dual) = groups.get(66 - degree) {
                assert_eq!(group.factors, dual.factors);
                assert_eq!(group.fields["common_base"], dual.fields["common_base"]);
                assert_eq!(group.fields["cover_base"], dual.fields["cover_base"]);
            }
            let prior = degree.checked_sub(1).map_or(0, |k| mod2_half_rank(k, rank));
            let mod2_betti = choose(65, degree) - mod2_half_rank(degree, rank) - prior;
            let next = groups.get(degree + 1);
            assert_eq!(
                free + p_torsion(Some(group), 2) + p_torsion(next, 2),
                mod2_betti
            );
            if rank != 32 {
                // For m=-5, or identity monodromy, the group difference vanishes mod 5.
                assert_eq!(
                    free + p_torsion(Some(group), 5) + p_torsion(next, 5),
                    choose(65, degree)
                );
            }
            if rank == 0 {
                assert_eq!(free, choose(65, degree));
                assert!(group.factors.is_empty());
            }
        }
    }
    assert_eq!(reports[0].fields["free"], 1);
    assert_eq!(reports[65].fields["free"], 1);
    assert!(reports[0].factors.is_empty() && reports[65].factors.is_empty());
}

#[test]
fn zero_identity_and_nonuniform_torus_limits() {
    let output = stdout(interpret(&format!(
        r"{}
        rvKSExteriorGroupSmith[Mixed,32,0,1];
        rvKSExteriorGroupSmith[Mixed,32,0,0];
        rvKSExteriorGroupSmith[Mixed,32,7,33];
        rvKSMappingTorus[Mixed,3,3];
        rvKSMappingTorus[Mixed,-3,63];
        ",
        setup()
    )));
    let result: Vec<_> = output.lines().map(report).collect();
    assert!(result[0].factors.is_empty());
    assert_eq!(result[0].fields["rank"], 0);
    assert_eq!(result[1].fields["rank"], choose(64, 32));
    assert_eq!(result[1].fields["factorial"], 1);
    assert_eq!(result[1].fields["cover_base"], 1);
    assert_eq!(result[0], result[2]);
    assert_eq!(result[3].factors, [(1, 752), (2, 240)]);
    assert_eq!(result[3].factors, result[4].factors);
    for (source, expected) in [
        (
            "rvKSMappingTorus[Mixed,3,32];",
            "nonuniform Kuga–Satake spin Smith factors",
        ),
        (
            "rvKSMappingTorus[A,1,66];",
            "degree must be between 0 and 65",
        ),
        (
            "rvKSExteriorGroupSmith[A,65,1,1];",
            "degree must be between 0 and 64",
        ),
        (
            "rvog6cusp Q(h,e,f); show rvKSExteriorDifference[A,1,1](rvKSExteriorBasis[Q,1]);",
            "expected type RVKSExterior(P,1)",
        ),
    ] {
        let actual = diagnostic(interpret(&format!("{}{source}", setup())));
        assert!(actual.contains(expected), "{actual}");
    }
}
