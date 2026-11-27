mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-ks-exterior-monodromy.asl")
        .split_once("// Exterior monodromy")
        .unwrap()
        .0
}

#[test]
fn runnable_example_uses_both_spin_monodromies() {
    let report = stdout(example("hyperelliptic-ks-exterior-monodromy"));
    let lines: Vec<_> = report.lines().collect();
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[0], "rvKSExteriorData[A,1] = dimension=64, nilpotency_index=2, graded=[(0, 16), (1, 32), (2, 16)], jordan=[(1, 32), (2, 16)]");
    assert_eq!(lines[1], "rvKSExteriorData[B3,2] = dimension=2016, nilpotency_index=3, graded=[(0, 496), (2, 1024), (4, 496)], jordan=[(1, 528), (3, 496)]");
    assert!(lines[2].starts_with(
        "rvKSExteriorData[A,32] = dimension=1832624140942590534, nilpotency_index=17"
    ));
}

fn choose(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    (0..k).fold(1, |value, i| value * (n - i) as u128 / (i + 1) as u128)
}

#[test]
fn all_degrees_match_an_independent_binomial_character_and_jordan_dimension() {
    let mut source = setup().to_owned();
    let cases = [("A", 16), ("B3", 32)];
    for (name, _) in cases {
        for degree in 0..=64 {
            writeln!(source, "rvKSExteriorData[{name},{degree}];").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut actual = output.lines();
    for (name, rank) in cases {
        for degree in 0usize..=64 {
            let mut weight = [0; 65];
            for positive in 0..=rank {
                for negative in 0..=rank.min(degree.saturating_sub(positive)) {
                    if let Some(trivial) = degree.checked_sub(positive + negative) {
                        weight[32 + positive - negative] += choose(rank, positive)
                            * choose(rank, negative)
                            * choose(64 - 2 * rank, trivial);
                    }
                }
            }
            let graded: Vec<_> = weight
                .iter()
                .enumerate()
                .filter_map(|(j, &count)| {
                    (count != 0).then_some((degree as isize + j as isize - 32, count))
                })
                .collect();
            let jordan: Vec<_> = (0..=rank)
                .filter_map(|j| {
                    let count = weight[32 + j] - weight.get(34 + j).unwrap_or(&0);
                    (count != 0).then_some((j + 1, count))
                })
                .collect();
            let dimension = choose(64, degree);
            let index = 1 + rank.min(degree).min(64 - degree);
            assert_eq!(graded.iter().map(|(_, n)| n).sum::<u128>(), dimension);
            assert_eq!(
                jordan
                    .iter()
                    .map(|&(size, n)| size as u128 * n)
                    .sum::<u128>(),
                dimension
            );
            assert_eq!(jordan.last().unwrap().0, index);
            assert_eq!(actual.next().unwrap(), format!("rvKSExteriorData[{name},{degree}] = dimension={dimension}, nilpotency_index={index}, graded={graded:?}, jordan={jordan:?}"));
        }
    }
    assert!(actual.next().is_none());
}

#[test]
fn grade_one_is_the_exact_spin_representation_under_both_maps() {
    let mut source = setup().to_owned();
    for (index, mask) in (0u32..128)
        .filter(|mask| mask.count_ones().is_multiple_of(2))
        .enumerate()
    {
        writeln!(
            source,
            "let s{index}=rvOG6CliffordBasis[P,{mask}]; let x{index}=rvKSExteriorBasis[P,{}];",
            1u64 << index
        )
        .unwrap();
        writeln!(source, "assert rvKSExteriorEmbed(s{index})==x{index}; assert rvKSExteriorToSpin(x{index})==s{index};").unwrap();
        for name in ["A", "B3"] {
            writeln!(source, "assert rvKSExteriorToSpin(rvKSExteriorPower[{name},1](x{index}))==rvOG6CliffordMul(rvOG6SpinLog[{name}],s{index});").unwrap();
            writeln!(source, "assert rvKSExteriorToSpin(rvKSExteriorExp[{name},-5](x{index}))==rvOG6CliffordMul(rvOG6SpinExp[{name},-5],s{index});").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), "");
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

fn insert_mod_101(pivots: &mut BTreeMap<u64, Column>, mut column: Column) -> bool {
    column.retain(|_, value| {
        *value = value.rem_euclid(101);
        *value != 0
    });
    while let Some((&position, &coefficient)) = column.first_key_value() {
        let Some(prior) = pivots.get(&position) else {
            pivots.insert(position, column);
            return true;
        };
        let pivot = prior[&position];
        for entry in column.values_mut() {
            *entry *= pivot;
        }
        for (&row, &entry) in prior {
            *column.entry(row).or_default() -= coefficient * entry;
        }
        column.retain(|_, value| {
            *value = value.rem_euclid(101);
            *value != 0
        });
    }
    false
}

#[test]
fn full_degree_two_actions_have_the_reported_ranks_for_both_types() {
    for (name, rank, middle) in [("A", 16, 32), ("B3", 32, 0)] {
        for (order, expected) in [
            (1, rank * (rank - 1) + rank * middle),
            (2, rank * (rank - 1) / 2),
        ] {
            let mut source = setup().to_owned();
            for i in 0..64 {
                for j in i + 1..64 {
                    let mask = (1u64 << i) | (1u64 << j);
                    writeln!(
                        source,
                        "show rvKSExteriorPower[{name},{order}](rvKSExteriorBasis[P,{mask}]);"
                    )
                    .unwrap();
                }
            }
            let output = stdout(interpret(&source));
            assert_eq!(output.lines().count(), 2016);
            let mut pivots = BTreeMap::new();
            for line in output.lines() {
                insert_mod_101(&mut pivots, sparse_column(line));
            }
            assert_eq!(pivots.len(), expected, "{name} order {order}");
        }
    }
}

#[test]
fn higher_wedge_actions_match_leibniz_and_the_exact_exponential() {
    let mut source = setup().to_owned();
    for index in 0..64 {
        writeln!(
            source,
            "show rvKSExteriorPower[B3,1](rvKSExteriorBasis[P,{}]);",
            1u64 << index
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let mut pivots = BTreeMap::new();
    let indices: Vec<_> = output
        .lines()
        .enumerate()
        .filter_map(|(index, line)| {
            insert_mod_101(&mut pivots, sparse_column(line)).then_some(index)
        })
        .take(3)
        .collect();
    assert_eq!(indices.len(), 3);
    let source = format!(
        r"{}
        let x=rvKSExteriorBasis[P,{}];
        let y=rvKSExteriorBasis[P,{}];
        let third=rvKSExteriorBasis[P,{}];
        let dx=rvKSExteriorPower[B3,1](x);
        let dy=rvKSExteriorPower[B3,1](y);
        let dz=rvKSExteriorPower[B3,1](third);
        let xyz=rvKSExteriorCup(rvKSExteriorCup(x,y),third);
        let one=rvKSExteriorPower[B3,1](xyz);
        let two=rvKSExteriorPower[B3,2](xyz);
        let three=rvKSExteriorPower[B3,3](xyz);
        let assigned=rvKSExteriorCup(rvKSExteriorCup(dx,dy),dz);
        let direct=add(add(rvKSExteriorCup(rvKSExteriorCup(dx,y),third),rvKSExteriorCup(rvKSExteriorCup(x,dy),third)),rvKSExteriorCup(rvKSExteriorCup(x,y),dz));
        assert one==direct;
        assert three==rvKSExteriorScale[6](assigned);
        assert rvKSExteriorPower[B3,1](three)==zero[RVKSExterior(P,3)];
        let Taylor=add(add(xyz,rvKSExteriorScale[6](one)),add(rvKSExteriorScale[18](two),rvKSExteriorScale[36](three)));
        let exp=rvKSExteriorExp[B3,6](xyz);
        assert exp==Taylor;
        assert rvKSExteriorExp[B3,-6](exp)==xyz;
        assert rvKSExteriorExp[B3,3](exp)==rvKSExteriorExp[B3,9](xyz);
        show three;
        ",
        setup(),
        1u64 << indices[0],
        1u64 << indices[1],
        1u64 << indices[2]
    );
    assert!(!sparse_column(stdout(interpret(&source)).trim()).is_empty());
}

#[test]
fn high_degree_exponential_matches_the_dual_pairing_and_direct_derivation() {
    let mut source = format!(
        r"{}
        let omega=rvKSExteriorBasis[P,18446744073709551614];
        let omega2=rvKSExteriorBasis[P,9223372036854775807];
        let input=add(omega,rvKSExteriorScale[-2](omega2));
        let d=rvKSExteriorPower[B3,1](input);
        let high=rvKSExteriorExp[B3,-170141183460469231731687303715884105728](input);
        assert high==add(input,rvKSExteriorScale[-170141183460469231731687303715884105728](d));
        let moderate=rvKSExteriorExp[B3,7](input);
        assert rvKSExteriorExp[B3,-7](moderate)==input;
        ",
        setup()
    );
    for index in 0..64 {
        writeln!(
            source,
            "let x{index}=rvKSExteriorBasis[P,{}]; assert rvKSExteriorCup(moderate,rvKSExteriorExp[B3,7](x{index}))==rvKSExteriorCup(input,x{index});",
            1u64 << index,
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn full_width_masks_formal_zero_and_nominal_types() {
    let source = format!(
        r"{}
        let unit=rvKSExteriorBasis[P,0];
        let x=rvKSExteriorBasis[P,1];
        let high=rvKSExteriorBasis[P,9223372036854775808];
        let xz=rvKSExteriorBasis[P,9223372036854775809];
        let top=rvKSExteriorBasis[P,18446744073709551615];
        assert rvKSExteriorCup(x,high)==xz;
        assert rvKSExteriorCup(high,x)==neg(xz);
        assert rvKSExteriorCup(high,high)==zero[RVKSExterior(P,2)];
        assert rvKSExteriorCup(unit,high)==high;
        assert rvKSExteriorCup(top,high)==zero[RVKSExterior(P,65)];
        assert rvKSExteriorCup(high,zero[RVKSExterior(P,-1)])==zero[RVKSExterior(P,0)];
        assert rvKSExteriorPower[A,2](unit)==zero[RVKSExterior(P,0)];
        assert rvKSExteriorPower[B3,1](top)==zero[RVKSExterior(P,64)];
        assert rvKSExteriorPower[B3,0](top)==top;
        assert rvKSExteriorExp[B3,-170141183460469231731687303715884105728](top)==top;
        rvKSExteriorCoeff[18446744073709551615](top);
        rvKSExteriorCoeff[9223372036854775809](rvKSExteriorCup(high,x));
        ",
        setup()
    );
    assert_eq!(stdout(interpret(&source)), "rvKSExteriorCoeff[18446744073709551615] = 1\nrvKSExteriorCoeff[9223372036854775809] = -1\n");

    for (source, expected) in [
        ("rvKSExteriorData[A,65];", "between 0 and 64"),
        ("show rvKSExteriorBasis[P,18446744073709551616];", "unsigned 64-bit integer"),
        ("show rvKSExteriorEmbed(rvOG6CliffordBasis[P,1]);", "expected type RVOG6Clifford(P,0)"),
        ("show rvKSExteriorToSpin(rvKSExteriorBasis[P,3]);", "expected type RVKSExterior(P,1)"),
        ("rvog6cusp Q(h,e,f); show rvKSExteriorPower[A,1](rvKSExteriorBasis[Q,0]);", "expected type RVKSExterior(P,0)"),
        ("rvog6cusp Q(h,e,f); show rvKSExteriorCup(rvKSExteriorBasis[P,1],rvKSExteriorBasis[Q,2]);", "expected type RVKSExterior(P,1)"),
        ("show rvKSExteriorCup(zero[RVKSExterior(P,170141183460469231731687303715884105727)],rvKSExteriorBasis[P,1]);", "exterior degree exceeds the signed i128 range"),
    ] {
        let actual = diagnostic(interpret(&format!("{}{source}", setup())));
        assert!(actual.contains(expected), "{actual}");
    }
}
