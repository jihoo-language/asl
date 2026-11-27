mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

type Tail = Vec<usize>;
type Path = Vec<Tail>;

#[test]
fn runnable_example_distinguishes_mixed_paths_with_the_same_shape() {
    let output = stdout(example("deligne-vacillating-paths"));
    let fixed = concat!(
        "rvDelignePathEigenvalues[3, 6] = [15, 10, 9, 5, 5, 2, 0]\n",
        "rvDelignePathEigenvalues[3, 6] = [15, 10, 9, 5, 3, 0, 0]\n",
        "rvPartitionTraceAt[6] = 16\n",
        "rvPartitionTraceAt[6] = 16\n",
        "rvKaroubiDimension[Mixed] = numerator=[(1,8), (2,-6), (3,1)], denominator=3\n",
        "rvPartitionEntry[6] = 1\n",
        "rvPartitionEntry[6] = 0\n",
    );
    assert!(output.starts_with(fixed), "{output}");
    assert!(output.contains("rvDelignePathCounts[3] = paths=31, dimension=203"));
    assert!(output.contains("rvDelignePathCounts[4] = paths=164, dimension=4140"));
}

fn partitions(weight: usize) -> Vec<Tail> {
    if weight == 0 {
        return vec![Vec::new()];
    }
    (1..=weight)
        .flat_map(|first| {
            partitions(weight - first)
                .into_iter()
                .filter(move |tail| tail.first().is_none_or(|&x| x <= first))
                .map(move |tail| std::iter::once(first).chain(tail).collect())
        })
        .collect()
}

fn contains(larger: &[usize], smaller: &[usize]) -> bool {
    smaller
        .iter()
        .enumerate()
        .all(|(i, &row)| larger.get(i).is_some_and(|&b| b >= row))
}

fn paths(rank: usize) -> Vec<Path> {
    let mut paths = vec![vec![Vec::new()]];
    for level in 1..=2 * rank {
        let candidates: Vec<_> = (0..=level / 2).flat_map(partitions).collect();
        let mut next = Vec::new();
        for path in paths {
            let parent = path.last().unwrap();
            let parent_weight: usize = parent.iter().sum();
            for child in &candidates {
                let size: usize = child.iter().sum();
                let neighbor = if level.is_multiple_of(2) {
                    size == parent_weight + 1 && contains(child, parent)
                } else {
                    size + 1 == parent_weight && contains(parent, child)
                };
                if parent == child || neighbor {
                    let mut successor = path.clone();
                    successor.push(child.clone());
                    next.push(successor);
                }
            }
        }
        paths = next;
    }
    paths
}

fn literal(path: &Path, rank: usize) -> String {
    let tails = path
        .iter()
        .map(|tail| format!("{tail:?}"))
        .collect::<Vec<_>>()
        .join(",");
    format!("rvDelignePathProjector[{rank}]{{{tails}}}")
}

fn hook_dimension(tail: &[usize], n: usize) -> usize {
    let partition: Vec<_> = std::iter::once(n - tail.iter().sum::<usize>())
        .chain(tail.iter().copied())
        .collect();
    let product: usize = partition
        .iter()
        .enumerate()
        .flat_map(|(i, &row)| {
            let partition = &partition;
            (0..row).map(move |column| {
                row - column
                    + partition[i + 1..]
                        .iter()
                        .filter(|&&lower| lower > column)
                        .count()
            })
        })
        .product();
    (1..=n).product::<usize>() / product
}

#[test]
fn the_complete_rank_three_paths_are_symbolically_orthogonal_and_have_the_expected_ranks() {
    let paths = paths(3);
    assert_eq!(paths.len(), 31);
    let mut input = String::from("let zeroLocal = zero[RVPartitionLocal(3)];\n");
    for (i, path) in paths.iter().enumerate() {
        writeln!(input, "let p{i} = {};", literal(path, 3)).unwrap();
        writeln!(input, "rvPartitionTraceAt[6](p{i});").unwrap();
        writeln!(input, "assert rvPartitionCompose(p{i},p{i}) == p{i};").unwrap();
        writeln!(input, "assert rvPartitionReflect(p{i}) == p{i};").unwrap();
        for j in 0..i {
            writeln!(input, "assert rvPartitionCompose(p{i},p{j}) == zeroLocal;").unwrap();
        }
    }
    let sum = (0..paths.len())
        .map(|i| format!("p{i}"))
        .reduce(|a, b| format!("add({a},{b})"))
        .unwrap();
    writeln!(
        input,
        "assert {sum} == rvPartitionLocal(rvPartitionIdentity[3]);"
    )
    .unwrap();
    let output = stdout(interpret(&input));
    let ranks: Vec<usize> = output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect();
    let expected: Vec<_> = paths
        .iter()
        .map(|path| hook_dimension(path.last().unwrap(), 6))
        .collect();
    assert_eq!(ranks, expected);
    assert_eq!(ranks.iter().sum::<usize>(), 216);
}

#[test]
fn compressed_counts_match_the_independent_vacillating_paths() {
    let input = "rvDelignePathCounts[0]; rvDelignePathCounts[1]; rvDelignePathCounts[2]; rvDelignePathCounts[3]; rvDelignePathCounts[4];";
    let mut expected = String::new();
    for r in 0..=4 {
        let paths = paths(r);
        let mut counts = BTreeMap::<Tail, usize>::new();
        for path in &paths {
            *counts.entry(path.last().unwrap().clone()).or_default() += 1;
        }
        let dim: usize = counts.values().map(|n| n * n).sum();
        let terms = counts
            .into_iter()
            .map(|(tail, count)| format!("({tail:?},{count})"))
            .collect::<Vec<_>>();
        writeln!(
            expected,
            "rvDelignePathCounts[{r}] = paths={}, dimension={dim}, multiplicities=[{}]",
            paths.len(),
            terms.join(", ")
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(input)), expected);
}

#[test]
fn finite_tower_operators_match_direct_transposition_matrices() {
    let n = 6;
    let mut input = String::new();
    let mut expected = Vec::new();
    let lowers = [[0, 1, 2], [0, 0, 1], [2, 4, 2], [3, 4, 5]];
    for level in 0..=6usize {
        writeln!(input, "let op{level} = rvDeligneTowerOperator[3,{level}];").unwrap();
        let active = level / 2;
        for lower in lowers {
            let mut counts = BTreeMap::new();
            for a in 0..n {
                for b in a + 1..n {
                    if !level.is_multiple_of(2) && [a, b].contains(&lower[active]) {
                        continue;
                    }
                    let mut upper = lower;
                    for letter in &mut upper[..active] {
                        *letter = if *letter == a {
                            b
                        } else if *letter == b {
                            a
                        } else {
                            *letter
                        };
                    }
                    *counts.entry(upper).or_insert(0) += 1;
                }
            }
            for a in 0..n {
                for b in 0..n {
                    for c in 0..n {
                        let upper = [a, b, c];
                        writeln!(
                            input,
                            "rvPartitionEntry[{n},{upper:?},{lower:?}](op{level});"
                        )
                        .unwrap();
                        expected.push(counts.get(&upper).copied().unwrap_or(0));
                    }
                }
            }
        }
    }
    let output = stdout(interpret(&input));
    assert_eq!(
        output
            .lines()
            .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
            .collect::<Vec<usize>>(),
        expected
    );
}

#[test]
fn every_doubled_level_has_the_claimed_symbolic_eigenvalue_on_mixed_paths() {
    let first = vec![
        vec![],
        vec![],
        vec![1],
        vec![1],
        vec![2],
        vec![2],
        vec![2, 1],
    ];
    let second = vec![
        vec![],
        vec![],
        vec![1],
        vec![1],
        vec![1, 1],
        vec![1, 1],
        vec![2, 1],
    ];
    let mut source = String::new();
    for j in 0..=6 {
        writeln!(source, "let op{j} = rvDeligneTowerOperator[3,{j}];").unwrap();
    }
    for (i, path) in [first, second].iter().enumerate() {
        writeln!(source, "let p{i} = {};", literal(path, 3)).unwrap();
        for (j, tail) in path.iter().enumerate() {
            let ordinary_content = |t: i128| {
                let m = tail.iter().sum::<usize>() as i128;
                let first = t - m - i128::from(!j.is_multiple_of(2));
                std::iter::once(first)
                    .chain(tail.iter().map(|&v| v as i128))
                    .enumerate()
                    .map(|(row, width)| {
                        (0..width)
                            .map(|column| 2 * (column - row as i128))
                            .sum::<i128>()
                    })
                    .sum::<i128>()
            };
            let y = ordinary_content(10);
            let w = (ordinary_content(12) - 2 * ordinary_content(11) + y) / 2;
            let v = ordinary_content(11) - y - 21 * w;
            let u = y - 10 * v - 100 * w;
            writeln!(source,"assert rvPartitionCompose(rvPartitionScale[2,0](op{j}),p{i}) == add(rvPartitionScale[{u},0](p{i}),add(rvPartitionScale[{v},1](p{i}),rvPartitionScale[{w},2](p{i}))); ").unwrap();
        }
    }
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn rank_two_unique_paths_recover_the_previous_primitive_projectors() {
    let output = stdout(interpret(
        r"
        let h = rvDelignePathProjector[2]{[],[],[1],[1],[2]};
        let a = rvDelignePathProjector[2]{[],[],[1],[1],[1,1]};
        assert h == rvDeligneSquareProjector[5];
        assert a == rvDeligneSquareProjector[6];
        let root = rvDelignePathProjector[0]{[]};
        assert root == rvPartitionLocal(rvPartitionIdentity[0]);
        rvPartitionTraceAt[6](rvDeligneTowerOperator[0,0]);
        let one = rvPartitionLocal(rvPartitionIdentity[0]);
        let tMinus4 = add(rvPartitionScale[1,1](one),rvPartitionScale[-4,0](one));
        let rational = rvPartitionLocalDivide[0,0,0,0,1,1,1](tMinus4);
        rvPartitionEntry[4,[],[]](rational);
        rvPartitionEntry[5,[],[]](rational);
    ",
    ));
    assert_eq!(
        output,
        "rvPartitionTraceAt[6] = 15\nrvPartitionEntry[4] = 1/3\nrvPartitionEntry[5] = 1/6\n"
    );
}

#[test]
fn path_shape_branches_rank_and_extended_localization_are_checked() {
    for (query,expected) in [
        ("show rvDeligneTowerOperator[3,7];","level exceeds"),
        ("show rvDelignePathProjector[4]{[]};","limited to rank 3"),
        ("rvDelignePathCounts[17];","limited to rank 16"),
        ("show rvDelignePathProjector[2]{[],[],[1]};","2r+1 vertices"),
        ("show rvDelignePathProjector[1]{[1],[],[1]};","start at the empty tail"),
        ("show rvDelignePathProjector[2]{[],[],[1],[2],[2]};","invalid Deligne path branch"),
        ("show rvDelignePathProjector[2]{[],[],[1],[1],[1,2]};","invalid Deligne path tail"),
        ("show rvPartitionLocalDivide[0,0,0,0,1](rvPartitionLocal(rvPartitionIdentity[0]));","four or seven exponents"),
        ("rvPartitionEntry[4,[],[]](rvPartitionLocalDivide[0,0,0,0,0,0,1](rvPartitionLocal(rvPartitionIdentity[0])));","pole at t=4"),
    ] {
        let error=diagnostic(interpret(query));
        assert!(error.contains(expected),"{error}");
    }
}
