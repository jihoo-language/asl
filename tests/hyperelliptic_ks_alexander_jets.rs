mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-ks-alexander-jets.asl")
        .split_once("// Completed Alexander")
        .unwrap()
        .0
}

fn pairs(text: &str) -> Vec<(usize, u128)> {
    let entries = text.strip_prefix('[').unwrap().strip_suffix(']').unwrap();
    entries
        .split("), (")
        .filter(|entry| !entry.is_empty())
        .map(|entry| {
            let (a, b) = entry.trim_matches(['(', ')']).split_once(", ").unwrap();
            (a.parse().unwrap(), b.parse().unwrap())
        })
        .collect()
}

fn profile(line: &str) -> (u128, Vec<(usize, u128)>) {
    let (_, value) = line.split_once(" = dimension=").unwrap();
    let (dimension, tail) = value.split_once(", exponents=").unwrap();
    let (exponents, differential) = tail.split_once(", differential_ranks=").unwrap();
    assert_eq!(exponents, differential);
    (dimension.parse().unwrap(), pairs(exponents))
}

fn jets(line: &str) -> [u128; 3] {
    let (_, body) = line.split_once(" = ").unwrap();
    let values: Vec<_> = body
        .split(", ")
        .map(|field| field.split_once('=').unwrap().1.parse().unwrap())
        .collect();
    values.try_into().unwrap()
}

#[test]
fn example_distinguishes_modular_blocks_and_certifies_d2() {
    let output = stdout(example("hyperelliptic-ks-alexander-jets"));
    let mut lines = output.lines();
    assert_eq!(
        profile(lines.next().unwrap()),
        (2016, vec![(1, 528), (3, 496)])
    );
    assert_eq!(
        profile(lines.next().unwrap()),
        (2016, vec![(1, 32), (2, 992)])
    );
    assert_eq!(
        profile(lines.next().unwrap()),
        (
            1832624140942590534,
            vec![(1, 601080390), (2, 916312070170755072)]
        )
    );
    assert_eq!(jets(lines.next().unwrap()), [1584, 1520, 64]);
    assert_eq!(jets(lines.next().unwrap()), [2080, 2016, 64]);
    assert_eq!(
        lines.next().unwrap(),
        "rvKSAlexanderWitness[A,1,Q] = order=2, obstruction_pairing=-3, certified_nonzero=true"
    );
    assert_eq!(
        lines.next().unwrap(),
        "rvKSAlexanderWitness[A,1,2] = order=2, obstruction_pairing=1, certified_nonzero=true"
    );
    assert!(lines.next().is_none());
}

fn choose(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    (0..k).fold(1, |acc, i| acc * (n - i) as u128 / (i + 1) as u128)
}

fn sparse_column(line: &str) -> BTreeMap<u64, i128> {
    let (_, terms) = line.split_once("{[").unwrap();
    let terms = terms.strip_suffix("] / 1}").unwrap();
    terms
        .split("), (")
        .filter(|term| !term.is_empty())
        .map(|term| {
            let (i, c) = term.trim_matches(['(', ')']).split_once(", ").unwrap();
            (i.parse().unwrap(), c.parse().unwrap())
        })
        .collect()
}

fn actual_blocks(name: &str, m: i128, slots: [usize; 3]) -> Vec<Vec<Vec<i128>>> {
    let mut source = setup().to_owned();
    let mut masks = Vec::new();
    for size in 0..=slots.len() {
        let group: Vec<_> = (0..1 << size)
            .map(|pattern| {
                slots[..size]
                    .iter()
                    .enumerate()
                    .map(|(i, bit)| 1u64 << (bit + usize::from(pattern & (1 << i) != 0)))
                    .sum()
            })
            .collect();
        for &mask in &group {
            writeln!(
                source,
                "show rvKSExteriorExp[{name},{m}](rvKSExteriorBasis[P,{mask}]);"
            )
            .unwrap();
        }
        masks.push(group);
    }
    let output = stdout(interpret(&source));
    let mut columns = output.lines().map(sparse_column);
    let blocks = masks
        .iter()
        .map(|group| {
            group
                .iter()
                .map(|_| {
                    let column = columns.next().unwrap();
                    assert!(column.keys().all(|key| group.contains(key)));
                    group
                        .iter()
                        .map(|mask| column.get(mask).copied().unwrap_or(0))
                        .collect()
                })
                .collect()
        })
        .collect();
    assert!(columns.next().is_none());
    blocks
}

fn integer_gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

fn matrix_rank(mut matrix: Vec<Vec<i128>>, prime: Option<i128>) -> usize {
    let mut rank = 0;
    for col in 0..matrix.len() {
        let Some(row) = (rank..matrix.len()).find(|&i| matrix[i][col] != 0) else {
            continue;
        };
        matrix.swap(rank, row);
        let pivot = matrix[rank].clone();
        for row in &mut matrix[rank + 1..] {
            let factor = row[col];
            if factor == 0 {
                continue;
            }
            for (entry, &v) in row[col..].iter_mut().zip(&pivot[col..]) {
                *entry = pivot[col] * *entry - factor * v;
                if let Some(p) = prime {
                    *entry = entry.rem_euclid(p);
                }
            }
            if prime.is_none() {
                let content = row[col..].iter().copied().fold(0, integer_gcd);
                if content > 1 {
                    for value in &mut row[col..] {
                        *value /= content;
                    }
                }
            }
        }
        rank += 1;
    }
    rank
}

fn truncated_nullity(matrix: &[Vec<i128>], length: usize, prime: Option<i128>) -> u128 {
    let n = matrix.len();
    let mut truncated = vec![vec![0; n * length]; n * length];
    for power in 0..length {
        for (i, column) in matrix.iter().enumerate() {
            for (j, &entry) in column.iter().enumerate() {
                let mut constant = entry - i128::from(i == j);
                let mut linear = entry;
                if let Some(p) = prime {
                    constant = constant.rem_euclid(p);
                    linear = linear.rem_euclid(p);
                }
                truncated[power * n + j][power * n + i] = constant;
                if power + 1 < length {
                    truncated[(power + 1) * n + j][power * n + i] = linear;
                }
            }
        }
    }
    (n * length - matrix_rank(truncated, prime)) as u128
}

fn actual_fiber_nullity(
    blocks: &[Vec<Vec<i128>>],
    r: usize,
    length: usize,
    p: Option<i128>,
    rank: usize,
) -> u128 {
    (0..=r)
        .map(|s| {
            let copies: u128 = (0..=(r - s) / 2)
                .map(|b| {
                    choose(rank, s) * choose(rank - s, b) * choose(64 - 2 * rank, r - s - 2 * b)
                })
                .sum();
            copies * truncated_nullity(&blocks[s], length, p)
        })
        .sum()
}

#[test]
fn jets_match_direct_truncated_polynomial_matrices_from_actual_pbw_actions() {
    for (name, m, rank, slots, fields) in [
        (
            "A",
            -1,
            16,
            [0, 6, 10],
            &[("Q", None), ("2", Some(2)), ("3", Some(3))][..],
        ),
        ("B3", 1, 32, [6, 12, 20], &[("3", Some(3))][..]),
    ] {
        let blocks = actual_blocks(name, m, slots);
        let mut source = setup().to_owned();
        for (field, _) in fields {
            for r in 0..=3 {
                for length in 1..=5 {
                    writeln!(
                        source,
                        "rvKSAlexanderJets[{name},{m},{r},{length},{field}];"
                    )
                    .unwrap();
                }
            }
        }
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        for &(_, prime) in fields {
            for degree in 0usize..=3 {
                for length in 1..=5 {
                    let kernel = actual_fiber_nullity(&blocks, degree, length, prime, rank);
                    let cokernel = degree
                        .checked_sub(1)
                        .map_or(0, |k| actual_fiber_nullity(&blocks, k, length, prime, rank));
                    assert_eq!(
                        jets(lines.next().unwrap()),
                        [kernel + cokernel, kernel, cokernel]
                    );
                }
            }
        }
        assert!(lines.next().is_none());
    }
}

fn f2_rank(d: usize, pairs: usize) -> u128 {
    (1..=d.min(pairs))
        .map(|s| {
            let others = (0..=(pairs - s).min((d - s) / 2))
                .map(|b| choose(pairs - s, b) * choose(64 - 2 * pairs, d - s - 2 * b))
                .sum::<u128>();
            choose(pairs, s) * others * (1u128 << (s - 1))
        })
        .sum()
}

#[test]
fn all_degrees_obey_duality_and_the_modular_plateau() {
    let modes = [
        ("B3", 1, "Q", 33),
        ("B3", -1, "2", 2),
        ("B3", 1, "3", 3),
        ("Mixed", 1, "2", 2),
        ("B3", 3, "3", 1),
        ("A", 0, "Q", 1),
    ];
    let mut source = setup().to_owned();
    for (name, m, p, plateau) in modes {
        for degree in 0..=64 {
            writeln!(source, "rvKSAlexander[{name},{m},{degree},{p}];").unwrap();
        }
        for degree in 0..=65 {
            writeln!(source, "rvKSAlexanderJets[{name},{m},{degree},1,{p}]; rvKSAlexanderJets[{name},{m},{degree},{plateau},{p}]; rvKSAlexanderJets[{name},{m},{degree},340282366920938463463374607431768211455,{p}];").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    for (name, _, p, plateau) in modes {
        let profiles: Vec<_> = (0..=64).map(|_| profile(lines.next().unwrap())).collect();
        for (degree, (dim, exponents)) in profiles.iter().enumerate() {
            assert_eq!(*dim, choose(64, degree));
            assert_eq!(
                exponents.iter().map(|&(s, b)| s as u128 * b).sum::<u128>(),
                *dim
            );
            assert_eq!(profiles[64 - degree].1, *exponents);
            assert!(exponents.iter().all(|&(size, _)| size <= plateau));
            if p == "2" {
                let rank = f2_rank(degree, if name == "Mixed" { 16 } else { 32 });
                let expected: Vec<_> = [(1, *dim - 2 * rank), (2, rank)]
                    .into_iter()
                    .filter(|entry| entry.1 != 0)
                    .collect();
                assert_eq!(*exponents, expected);
            }
            if plateau == 1 {
                assert_eq!(*exponents, [(1, *dim)]);
            }
        }
        let mut first_jet = Vec::new();
        for degree in 0..=65 {
            let first = jets(lines.next().unwrap());
            let final_jet = jets(lines.next().unwrap());
            let very_long = jets(lines.next().unwrap());
            assert_eq!(final_jet, very_long);
            assert_eq!(final_jet[0], choose(65, degree));
            assert_eq!(first[0], first[1] + first[2]);
            if matches!(degree, 0 | 65) {
                assert_eq!(first[0], 1);
            }
            first_jet.push(first[0]);
        }
        for degree in 0..=65 {
            assert_eq!(first_jet[degree], first_jet[65 - degree]);
        }
    }
    assert!(lines.next().is_none());
}

#[test]
fn lift_certificate_discards_earlier_indeterminacy_and_detects_d3() {
    let input = format!(
        r"{}
        let ev=rvKSExteriorBasis[P,2];
        let topPre=rvKSExteriorBasis[P,18446744073709551614];
        rvKSAlexanderWitness[A,1,Q](ev)(topPre);
        let primal=rvKSExteriorBasis[P,65];
        let dual=rvKSExteriorBasis[P,18446744073709551485];
        let v0=rvKSExteriorDifference[A,1,2](primal);
        let v1=neg(rvKSExteriorExp[A,1](rvKSExteriorDifference[A,1,1](primal)));
        let v2=rvKSExteriorExp[A,2](primal);
        let w0=rvKSExteriorDifference[A,-1,2](dual);
        let w1=neg(rvKSExteriorExp[A,-1](rvKSExteriorDifference[A,-1,1](dual)));
        let w2=rvKSExteriorExp[A,-2](dual);
        rvKSAlexanderWitness[A,1,Q](v0,v1,v2)(w0,w1,w2);
        rvKSAlexanderWitness[A,1,2](v0,v1,v2)(w0,w1,w2);
        let fv=rvKSExteriorBasis[P,1];
        rvKSAlexanderWitness[A,-170141183460469231731687303715884105728,2](fv)(topPre);
        let dualFrom=rvKSExteriorBasis[P,18446744073709551613];
        rvKSAlexanderWitness[A,-170141183460469231731687303715884105728,5](ev,neg(fv))(topPre,dualFrom);
        ",
        setup()
    );
    let output = stdout(interpret(&input));
    assert_eq!(output, concat!(
        "rvKSAlexanderWitness[A,1,Q] = order=1, obstruction_pairing=0, certified_nonzero=false\n",
        "rvKSAlexanderWitness[A,1,Q] = order=3, obstruction_pairing=-18, certified_nonzero=true\n",
        "rvKSAlexanderWitness[A,1,2] = order=3, obstruction_pairing=0, certified_nonzero=false\n",
        "rvKSAlexanderWitness[A,-170141183460469231731687303715884105728,2] = order=1, obstruction_pairing=1, certified_nonzero=true\n",
        "rvKSAlexanderWitness[A,-170141183460469231731687303715884105728,5] = order=2, obstruction_pairing=4, certified_nonzero=true\n",
    ));
}

#[test]
fn invalid_fields_degrees_and_real_lift_equations_are_diagnosed() {
    let invalid = [
        ("rvKSAlexander[A,1,2,4];", "characteristic must be prime"),
        ("rvKSAlexander[A,1,65,Q];", "fiber degree must be between 0 and 64"),
        ("rvKSAlexanderJets[A,1,66,1,Q];", "jet degree must be between 0 and 65"),
        ("rvKSAlexanderJets[A,1,2,0,Q];", "jet length must be positive"),
        ("rvKSAlexanderWitness[A,1,Q]()();", "same length between 1 and 33"),
        ("rvKSAlexanderWitness[A,1,Q](rvKSExteriorBasis[P,1])(rvKSExteriorBasis[P,18446744073709551614]);", "forward Alexander lift: coefficient 0"),
        ("rvKSAlexanderWitness[A,1,Q](rvKSExteriorBasis[P,2])(rvKSExteriorBasis[P,18446744073709551613]);", "adjoint Alexander lift: coefficient 0"),
        ("rvKSAlexanderWitness[A,1,Q](rvKSExteriorBasis[P,2])(rvKSExteriorBasis[P,2]);", "expected type RVKSExterior(P,63)"),
    ];
    for (command, expected) in invalid {
        let actual = diagnostic(interpret(&format!("{}{command}", setup())));
        assert!(actual.contains(expected), "{actual}");
    }
}
