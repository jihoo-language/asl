mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const HUGE: u128 = u128::MAX;

#[test]
fn example_exhibits_loops_orbits_noncommutativity_and_the_finite_centralizer() {
    assert_eq!(
        stdout(example("deligne-partition-algebra")),
        concat!(
            "rvPartitionTrace = [(2,1)]\n",
            "rvPartitionEntry[3] = 3\n",
            "rvPartitionEntry[3] = 1\n",
            "rvPartitionEntry[3] = 0\n",
            "RVPartitionAlgebra(2){([0, 1, 1, 2],0,1)}\n",
            "RVPartitionAlgebra(2){([0, 1, 2, 0],0,1)}\n",
            "rvPartitionDimensions[2, 2] = generic=15, finite=8\n",
            "rvPartitionCentralizerCheck[2, 2] = generic=15, finite=8, multiplicity squares=8\n",
        )
    );
}

#[test]
fn concatenation_juxtaposition_and_reflection_obey_the_algebra_laws() {
    let output = stdout(interpret(
        r"
        let unit = rvPartitionIdentity[0];
        let e = rvPartitionIdentity[1];
        let j = rvPartitionDiagram[1]{888,999};
        let p = rvPartitionDiagram[2]{0,1,1,0};
        let a = rvPartitionTensor(e,j);
        let b = rvPartitionTensor(j,e);
        let c = rvPartitionDiagram[2]{0,0,1,2};
        assert a == rvPartitionDiagram[2]{0,1,0,2};
        assert rvPartitionReflect(c) == rvPartitionDiagram[2]{0,1,2,2};
        assert rvPartitionTensor(e,e) == rvPartitionIdentity[2];
        assert rvPartitionTensor(unit,c) == c;
        assert rvPartitionTensor(c,unit) == c;
        assert rvPartitionCompose(p,p) == rvPartitionIdentity[2];
        assert rvPartitionCompose(a,rvPartitionIdentity[2]) == a;
        assert rvPartitionCompose(rvPartitionIdentity[2],a) == a;
        assert rvPartitionCompose(p,a) == rvPartitionDiagram[2]{0,1,1,2};
        assert rvPartitionCompose(a,p) == rvPartitionDiagram[2]{0,1,2,0};
        let ap = add(rvPartitionScale[2,1](a),neg(p));
        assert rvPartitionCompose(rvPartitionCompose(ap,b),c) == rvPartitionCompose(ap,rvPartitionCompose(b,c));
        assert rvPartitionCompose(add(a,b),c) == add(rvPartitionCompose(a,c),rvPartitionCompose(b,c));
        assert rvPartitionReflect(rvPartitionCompose(ap,c)) == rvPartitionCompose(rvPartitionReflect(c),rvPartitionReflect(ap));
        assert rvPartitionCompose(rvPartitionTensor(e,j),rvPartitionTensor(j,e)) == rvPartitionTensor(rvPartitionCompose(e,j),rvPartitionCompose(j,e));
        assert rvPartitionTensor(rvPartitionTensor(e,j),e) == rvPartitionTensor(e,rvPartitionTensor(j,e));
        assert rvPartitionCompose(zero[RVPartitionAlgebra(2)],ap) == zero[RVPartitionAlgebra(2)];
        assert add(ap,neg(ap)) == zero[RVPartitionAlgebra(2)];
        rvPartitionTrace(rvPartitionCompose(ap,c));
        rvPartitionTrace(rvPartitionCompose(c,ap));
        rvPartitionTrace(rvPartitionIdentity[3]);
        rvPartitionTrace(rvPartitionIdentity[0]);
    ",
    ));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines[0], lines[1]);
    assert_eq!(
        &lines[2..],
        ["rvPartitionTrace = [(3,1)]", "rvPartitionTrace = [(0,1)]"]
    );
}

fn set_partitions(n: usize) -> Vec<Vec<usize>> {
    fn visit(n: usize, blocks: usize, labels: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if labels.len() == n {
            out.push(labels.clone());
            return;
        }
        for b in 0..=blocks {
            labels.push(b);
            visit(n, blocks.max(b + 1), labels, out);
            labels.pop();
        }
    }
    let mut out = Vec::new();
    visit(n, 0, &mut Vec::new(), &mut out);
    out
}

fn literal(labels: &[usize], operation: &str) -> String {
    let rank = labels.len() / 2;
    let labels = labels
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",");
    format!("{operation}[{rank}]{{{labels}}}")
}

fn matches(labels: &[usize], upper: &[usize], lower: &[usize], exact: bool) -> bool {
    let tuple: Vec<_> = upper.iter().chain(lower).copied().collect();
    (0..labels.len()).all(|i| {
        (0..i).all(|j| {
            if exact {
                (labels[i] == labels[j]) == (tuple[i] == tuple[j])
            } else {
                labels[i] != labels[j] || tuple[i] == tuple[j]
            }
        })
    })
}

fn tuples(n: usize) -> Vec<[usize; 2]> {
    (0..n).flat_map(|i| (0..n).map(move |j| [i, j])).collect()
}

fn value(line: &str) -> i128 {
    line.rsplit_once(" = ").unwrap().1.parse().unwrap()
}

fn trace_value(line: &str, n: i128) -> i128 {
    let terms = line
        .strip_prefix("rvPartitionTrace = [")
        .unwrap()
        .strip_suffix(']')
        .unwrap();
    if terms.is_empty() {
        return 0;
    }
    terms
        .trim_start_matches('(')
        .trim_end_matches(')')
        .split("), (")
        .map(|term| {
            let (power, coefficient) = term.split_once(',').unwrap();
            coefficient.parse::<i128>().unwrap() * n.pow(power.parse().unwrap())
        })
        .sum()
}

#[test]
fn every_rank_two_product_and_trace_agrees_with_explicit_finite_matrices() {
    let diagrams = set_partitions(4);
    for n in [2, 3] {
        let coordinates = tuples(n);
        let mut script = String::new();
        for (i, d) in diagrams.iter().enumerate() {
            writeln!(script, "let d{i} = {};", literal(d, "rvPartitionDiagram")).unwrap();
        }
        let pairs: Vec<_> = if n == 2 {
            (0..diagrams.len())
                .flat_map(|i| (0..diagrams.len()).map(move |j| (i, j)))
                .collect()
        } else {
            vec![(0, 3), (3, 0), (3, 7), (5, 11), (11, 5), (13, 14), (14, 14)]
        };
        let mut expected = Vec::new();
        for (k, (i, j)) in pairs.into_iter().enumerate() {
            writeln!(script, "let p{k} = rvPartitionCompose(d{i},d{j});").unwrap();
            let mut trace = 0;
            for upper in &coordinates {
                for lower in &coordinates {
                    let entry = coordinates
                        .iter()
                        .filter(|middle| {
                            matches(&diagrams[i], upper, *middle, false)
                                && matches(&diagrams[j], *middle, lower, false)
                        })
                        .count() as i128;
                    writeln!(script, "rvPartitionEntry[{n},{upper:?},{lower:?}](p{k});").unwrap();
                    expected.push((false, entry));
                    if upper == lower {
                        trace += entry;
                    }
                }
            }
            writeln!(script, "rvPartitionTrace(p{k});").unwrap();
            expected.push((true, trace));
        }
        let output = stdout(interpret(&script));
        let lines: Vec<_> = output.lines().collect();
        assert_eq!(lines.len(), expected.len());
        for (line, (is_trace, count)) in lines.into_iter().zip(expected) {
            assert_eq!(
                if is_trace {
                    trace_value(line, n as i128)
                } else {
                    value(line)
                },
                count,
                "{line}"
            );
        }
    }
}

#[test]
fn mobius_orbit_elements_are_indicators_and_the_kernel_is_counted_by_blocks() {
    let diagrams = set_partitions(4);
    let mut script = String::new();
    let mut expected = Vec::new();
    for (i, d) in diagrams.iter().enumerate() {
        writeln!(script, "let a{i} = {};", literal(d, "rvPartitionOrbit")).unwrap();
        let blocks = d.iter().max().unwrap() + 1;
        for n in [1, 2, 3, 4] {
            let coordinates = tuples(n);
            let mut nonzero = false;
            for upper in &coordinates {
                for lower in &coordinates {
                    writeln!(script, "rvPartitionEntry[{n},{upper:?},{lower:?}](a{i});").unwrap();
                    let entry = matches(d, upper, lower, true);
                    expected.push(i128::from(entry));
                    nonzero |= entry;
                }
            }
            assert_eq!(nonzero, blocks <= n);
        }
    }
    script.push_str("rvPartitionCoefficient[[0,0,0,0],0](rvPartitionOrbit[2]{0,1,2,3});");
    expected.push(-6);
    assert_eq!(
        stdout(interpret(&script))
            .lines()
            .map(value)
            .collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn centralizer_dimensions_match_independent_orbit_counts_and_stable_multiplicities() {
    let mut script = String::new();
    let mut expected = String::new();
    for r in 0..=3 {
        let diagrams = set_partitions(2 * r);
        for n in 0..=2 * r + 1 {
            let degree = diagrams
                .iter()
                .filter(|labels| labels.iter().max().is_none_or(|&label| label < n))
                .count();
            writeln!(
                script,
                "rvPartitionDimensions[{r},{n}]; rvPartitionCentralizerCheck[{r},{n}];"
            )
            .unwrap();
            writeln!(
                expected,
                "rvPartitionDimensions[{r}, {n}] = generic={}, finite={degree}",
                diagrams.len()
            )
            .unwrap();
            writeln!(expected,"rvPartitionCentralizerCheck[{r}, {n}] = generic={}, finite={degree}, multiplicity squares={degree}",diagrams.len()).unwrap();
        }
    }
    script.push_str("rvPartitionDimensions[32,64];");
    let bell = "172134143357358850934369963665272571125557575184049758045339873395";
    writeln!(
        expected,
        "rvPartitionDimensions[32, 64] = generic={bell}, finite={bell}"
    )
    .unwrap();
    assert_eq!(stdout(interpret(&script)), expected);
}

#[test]
fn sparse_polynomial_entries_handle_zero_large_axes_and_exact_cancellation() {
    let almost = HUGE - 1;
    let script = format!(
        r"
        let u = rvPartitionIdentity[0];
        let c = add(rvPartitionScale[1,{HUGE}](u),rvPartitionScale[-2,{almost}](u));
        rvPartitionEntry[2,[],[]](c);
        rvPartitionEntry[1,[],[]](c);
        rvPartitionEntry[0,[],[]](u);
        rvPartitionEntry[0,[],[]](c);
        rvPartitionEntry[{HUGE},[],[]](rvPartitionScale[1,2](u));
        rvPartitionCoefficient[[9,8,9,8],0](rvPartitionIdentity[2]);
    "
    );
    let output = stdout(interpret(&script));
    let results = output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1)
        .collect::<Vec<_>>();
    assert_eq!(
        results,
        [
            "0",
            "-1",
            "1",
            "0",
            "115792089237316195423570985008687907852589419931798687112530834793049593217025",
            "1"
        ]
    );
}

#[test]
fn wrong_rank_overflow_and_orbit_budgets_have_direct_diagnostics() {
    for (source,expected) in [
        ("show rvPartitionDiagram[1]{0};","exactly 2r"),
        ("show rvPartitionIdentity[33];","rank exceeds"),
        ("show rvPartitionCompose(rvPartitionIdentity[1],rvPartitionIdentity[2]);","equal ranks"),
        ("show rvPartitionTensor(rvPartitionIdentity[17],rvPartitionIdentity[16]);","rank exceeds"),
        ("show rvPartitionOrbit[5]{0,1,2,3,4,5,6,7,8,9};","coarsenings exceed"),
        ("rvPartitionCentralizerCheck[9,2];","exceeds rank 8"),
        ("rvPartitionTrace(S.hi);","RVPartitionAlgebra"),
        ("rvPartitionEntry[3,[],[0]](rvPartitionIdentity[1]);","wrong tensor rank"),
        ("rvPartitionEntry[0,[0],[0]](rvPartitionIdentity[1]);","outside 0..N"),
        ("rvPartitionCoefficient[[0],0](rvPartitionIdentity[1]);","exactly 2r"),
        ("rvPartitionTrace(rvPartitionScale[1,340282366920938463463374607431768211455](rvPartitionIdentity[1]));","t degree exceeds"),
        ("show rvPartitionScale[2,0](rvPartitionScale[170141183460469231731687303715884105727,0](rvPartitionIdentity[1]));","signed 128-bit range"),
    ] {
        let got=diagnostic(interpret(source));
        assert!(got.contains(expected),"{got}");
    }
}
