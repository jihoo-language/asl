mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn example_reports_both_actual_massey_witnesses_and_the_period() {
    let output = stdout(example("deligne-yoneda-massey"));
    assert!(output.starts_with(concat!(
        "rvPartitionResolution[2,0,2,1] = syzygies=[(2, 1), (0, 2), (0, 1), (2, 2)], repeatFrom=0, period=4\n",
        "rvPartitionGlobalDimension[2,0] = infinity\n",
        "rvPartitionGlobalDimension[2,2] = 2\n",
    )), "{output}");
    assert!(output.contains("rvPartitionMassey[2,0,0,1,0] = left=[([(2, 2)], 1)], right=[([(0, 2)], 1)], cycle=[([(2, 1), (0, 2)], 1), ([(2, 2), (2, 1)], 1)], class=RVPartitionExt(2,0,2,0,2){-1}, indeterminacy=0\n"), "{output}");
    assert!(
        output.contains("class=RVPartitionExt(2,0,0,2,2){-1}, indeterminacy=0\n"),
        "{output}"
    );
    assert!(
        output.ends_with("rvPartitionBarExtDim[2,0,2,0,2] = 1\n"),
        "{output}"
    );
}

fn cyclic_head(start: usize, degree: usize) -> usize {
    if degree.div_ceil(2).is_multiple_of(2) {
        start
    } else {
        2 - start
    }
}

#[test]
fn ordinary_and_bar_ext_agree_in_every_low_degree_small_fiber() {
    let mut source = String::new();
    let mut expected = Vec::new();
    for (n, vertices) in [(0, 3), (1, 4), (2, 4), (3, 4)] {
        for i in 0..vertices {
            let finite = match (n, i) {
                (1, 3) => vec![3, 2, 3],
                (1, 2) => vec![2, 3],
                (2, 2) => vec![2, 1, 2],
                (2, 1) => vec![1, 2],
                _ => vec![i],
            };
            for j in 0..vertices {
                for degree in 0..=6 {
                    let head = if n == 0 && i != 1 {
                        Some(cyclic_head(i, degree))
                    } else {
                        finite.get(degree).copied()
                    };
                    writeln!(source, "rvPartitionBarExtDim[2,{n},{i},{j},{degree}];").unwrap();
                    expected.push(usize::from(head == Some(j)));
                }
            }
        }
    }
    let output = stdout(interpret(&source));
    assert_eq!(
        output
            .lines()
            .map(|line| line.rsplit_once(" = ").unwrap().1.parse::<usize>().unwrap())
            .collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn yoneda_products_match_the_periodic_cyclic_pattern_and_rational_scalars() {
    let mut source = String::new();
    for i in [0, 2] {
        for n in 0..=5 {
            let j = cyclic_head(i, n);
            for m in 0..=5 {
                let k = cyclic_head(j, m);
                let grade = format!("2,0,{i},{k},{}", n + m);
                let rhs = if n.is_multiple_of(2) || m.is_multiple_of(2) {
                    format!("rvPartitionExtClass[{grade}]")
                } else {
                    format!("zero[RVPartitionExt({grade})]")
                };
                writeln!(source,"assert rvPartitionYoneda(rvPartitionExtClass[2,0,{j},{k},{m}],rvPartitionExtClass[2,0,{i},{j},{n}]) == {rhs};").unwrap();
            }
        }
    }
    source.push_str(concat!(
        "let a=rvPartitionExtClass[2,0,2,0,1];\n",
        "let b=rvPartitionExtClass[2,0,0,0,300000000000000000000000000000000000000];\n",
        "let c=rvPartitionYoneda(rvPartitionExtScale[-3,2](b),rvPartitionExtScale[2,5](a));\n",
        "assert c == rvPartitionExtScale[-3,5](rvPartitionExtClass[2,0,2,0,300000000000000000000000000000000000001]);\n",
        "assert add(c,neg(c)) == zero[RVPartitionExt(2,0,2,0,300000000000000000000000000000000000001)];\n",
        "rvPartitionExtCoefficient(c);\n",
    ));
    assert_eq!(
        stdout(interpret(&source)),
        "rvPartitionExtCoefficient = -3/5\n"
    );
}

#[test]
fn actual_diagram_comparison_lifts_commute_with_differentials() {
    let mut source = String::from("let a=rvPartitionArrow[2,0,0];\nlet b=rvPartitionArrow[2,0,1];\nlet zero0=zero[RVPartitionFiber(2,0)];\n");
    source.push_str("assert rvPartitionDifferential[2,0,2,1,1] == a;\nassert rvPartitionDifferential[2,0,2,1,2] == rvPartitionCompose(a,b);\n");
    for k in 1..=8 {
        writeln!(source, "let d{k}=rvPartitionDifferential[2,0,2,1,{}];\nlet e{k}=rvPartitionDifferential[2,0,0,1,{k}];",k+1).unwrap();
        writeln!(source, "let f{k}=rvPartitionComparisonLift[2,0,2,0,1,{k}];").unwrap();
        writeln!(source, "assert rvPartitionCompose(d{k},rvPartitionComparisonLift[2,0,2,0,1,{}]) == rvPartitionCompose(f{k},e{k});",k-1).unwrap();
        if k > 1 {
            writeln!(
                source,
                "assert rvPartitionCompose(d{k},d{}) == zero0;",
                k - 1
            )
            .unwrap();
        }
    }
    source.push_str("assert rvPartitionComparisonLift[2,0,2,0,1,300000000000000000000000000000000000001] == f1;\n");
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn finite_linked_components_and_dual_numbers_have_actual_products() {
    let mut source = String::new();
    for (n, i, j) in [(1, 3, 2), (2, 2, 1)] {
        writeln!(source,"let a{n}=rvPartitionExtClass[2,{n},{i},{j},1];\nlet b{n}=rvPartitionExtClass[2,{n},{j},{i},1];").unwrap();
        writeln!(
            source,
            "assert rvPartitionYoneda(a{n},b{n}) == zero[RVPartitionExt(2,{n},{j},{j},2)];"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvPartitionYoneda(b{n},a{n}) == rvPartitionExtClass[2,{n},{i},{i},2];"
        )
        .unwrap();
        writeln!(
            source,
            "rvPartitionResolution[2,{n},{i},1];\nrvPartitionGlobalDimension[2,{n}];"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvPartitionDifferential[2,{n},{i},1,3] == zero[RVPartitionFiber(2,{n})];"
        )
        .unwrap();
    }
    source.push_str(concat!(
        "let u=rvPartitionExtClass[1,0,0,0,1];\n",
        "let w=rvPartitionExtClass[1,0,0,0,12345678901234567890];\n",
        "assert rvPartitionYoneda(w,u) == rvPartitionExtClass[1,0,0,0,12345678901234567891];\n",
        "rvPartitionResolution[1,0,0,1];\n",
        "rvPartitionGlobalDimension[1,0];\n",
        "rvPartitionGlobalDimension[2,3];\n",
    ));
    assert_eq!(stdout(interpret(&source)), concat!(
        "rvPartitionResolution[2,1,3,1] = syzygies=[(3, 1), (2, 1), (3, 2)], projectiveDimension=2\n",
        "rvPartitionGlobalDimension[2,1] = 2\n",
        "rvPartitionResolution[2,2,2,1] = syzygies=[(2, 1), (1, 1), (2, 2)], projectiveDimension=2\n",
        "rvPartitionGlobalDimension[2,2] = 2\n",
        "rvPartitionResolution[1,0,0,1] = syzygies=[(0, 1)], repeatFrom=0, period=1\n",
        "rvPartitionGlobalDimension[1,0] = infinity\n",
        "rvPartitionGlobalDimension[2,3] = 0\n",
    ));
}

#[test]
fn malformed_types_periods_and_massey_inputs_have_diagnostics() {
    for (source, expected) in [
        ("rvPartitionResolution[3,0,0,1];","limited to rank 2"),
        ("rvPartitionResolution[2,0,0,4];","module length must be in"),
        ("rvPartitionExtDim[2,0,9,0,2];","simple vertex is out of range"),
        ("rvPartitionBarExtDim[2,0,0,2,7];","limited to 6"),
        ("rvPartitionMassey[2,0,0,0,0];","must be composable"),
        ("rvPartitionMassey[2,0,1,0,8];","arrow index is out of range"),
        ("rvPartitionMassey[2,2,0,1,0];","both adjacent Yoneda products to vanish"),
        ("show rvPartitionExtClass[2,0,2,2,2];","Ext group has dimension zero"),
        ("show rvPartitionDifferential[2,0,2,1,0];","index must be positive"),
        ("show rvPartitionComparisonLift[2,0,2,2,2,0];","requires a nonzero Ext group"),
        ("show rvPartitionYoneda(rvPartitionExtClass[2,0,2,0,1],rvPartitionExtClass[2,0,2,0,1]);","matching intermediate simple"),
        ("show rvPartitionYoneda(rvPartitionExtClass[1,0,0,0,1],rvPartitionExtClass[1,0,0,0,340282366920938463463374607431768211455]);","Ext degree overflow"),
    ] {
        let message = diagnostic(interpret(source));
        assert!(message.contains(expected),"{source}: {message}");
    }
}
