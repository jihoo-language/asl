mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const V: &str = r"
    let u = rvStableSchurCharacter[[]];
    let w = rvStableSchurCharacter[[1]];
    let v = add(u,w);
";
const HUGE: u128 = u128::MAX;

#[test]
fn example_recovers_symmetric_exterior_and_inner_schur_decompositions() {
    assert_eq!(
        stdout(example("stable-inner-schur")),
        concat!(
            "RVStableSchurCharacter{([3],1), ([2],2), ([1, 1],1), ([1],4), ([],3)}\n",
            "RVStableSchurCharacter{([2, 1],1), ([2],2), ([1, 1],2), ([1],3), ([],1)}\n",
            "RVStableSchurCharacter{([1, 1, 1],1), ([1, 1],1)}\n",
            "rvSchurCharacterValue[[3, 2, 1]] = 3\n",
            "rvSchurCharacterValue[[3, 2, 1]] = -1\n",
        )
    );
}

fn values(output: &str) -> Vec<i128> {
    output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

fn partitions(n: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![Vec::new()];
    }
    (1..=n)
        .flat_map(|first| {
            partitions(n - first)
                .into_iter()
                .filter(move |tail| tail.first().is_none_or(|&size| size <= first))
                .map(move |tail| std::iter::once(first).chain(tail).collect())
        })
        .collect()
}

fn permutation(cycles: &[usize]) -> Vec<usize> {
    let mut out = Vec::new();
    for &cycle in cycles {
        let start = out.len();
        out.extend(start + 1..start + cycle);
        out.push(start);
    }
    out
}

fn power(p: &[usize], d: usize) -> Vec<usize> {
    (0..p.len())
        .map(|i| (0..d).fold(i, |current, _| p[current]))
        .collect()
}

fn cycles(p: &[usize]) -> Vec<usize> {
    let mut seen = vec![false; p.len()];
    let mut out = Vec::new();
    for i in 0..p.len() {
        if seen[i] {
            continue;
        }
        let mut current = i;
        let mut length = 0;
        while !seen[current] {
            seen[current] = true;
            length += 1;
            current = p[current];
        }
        out.push(length);
    }
    out.sort_unstable_by(|a, b| b.cmp(a));
    out
}

#[test]
fn finite_permutations_independently_verify_adams_values_and_tensor_eigen_characters() {
    let mut source = format!(
        r"{V}
        let y = add(rvStableSchurCharacter[[2]],add(rvStableSchurCharacterScale[-2](rvStableSchurCharacter[[1,1]]),w));
        let ay = rvStableSchurAdams[3](y);
        let av = rvStableSchurAdams[3](v);
        let exterior = rvStableSchurExterior[3](v);
        let symmetric = rvStableSchurSymmetric[3](v);
        let schur = rvStableInnerSchur[[2,1]](v);
        let negative = rvStableSchurExterior[3](neg(v));
    "
    );
    let mut expected = Vec::new();
    for n in [4, 6, 8] {
        for label in [
            "y",
            "ay",
            "av",
            "exterior",
            "symmetric",
            "schur",
            "negative",
        ] {
            writeln!(
                source,
                "let {label}{n} = rvStableSchurSpecialize[{n}]({label});"
            )
            .unwrap();
        }
        for cycle in partitions(n) {
            let p = permutation(&cycle);
            let p3 = power(&p, 3);
            let p3_type = cycles(&p3);
            writeln!(source,"rvSchurCharacterValue[{cycle:?}](ay{n}); rvSchurCharacterValue[{p3_type:?}](y{n});").unwrap();
            for label in ["av", "exterior", "symmetric", "schur", "negative"] {
                writeln!(source, "rvSchurCharacterValue[{cycle:?}]({label}{n});").unwrap();
            }
            // Each l-cycle contributes det(1+zP)=1-(-z)^l. It contributes
            // the inverse of (1-z^l) to the symmetric-power series.
            let (mut exterior, mut symmetric) = ([1, 0, 0, 0], [1, 0, 0, 0]);
            for &length in &cycle {
                for i in (length..4).rev() {
                    exterior[i] += if length.is_multiple_of(2) {
                        -exterior[i - length]
                    } else {
                        exterior[i - length]
                    };
                }
                for i in length..4 {
                    symmetric[i] += symmetric[i - length];
                }
            }
            let fixed = |permutation: &[usize]| {
                permutation
                    .iter()
                    .enumerate()
                    .filter(|(i, &j)| *i == j)
                    .count() as i128
            };
            let f1 = fixed(&p);
            let f3 = fixed(&p3);
            expected.push([
                f3,
                exterior[3],
                symmetric[3],
                (f1.pow(3) - f3) / 3,
                -symmetric[3],
            ]);
        }
    }
    let output = values(&stdout(interpret(&source)));
    assert_eq!(output.len(), expected.len() * 7);
    for (values, expected) in output.chunks_exact(7).zip(expected) {
        assert_eq!(values[0], values[1]);
        assert_eq!(values[2..], expected);
    }
}

#[test]
fn adams_are_internal_ring_maps_and_virtual_lambda_series_are_inverse() {
    let source = format!(
        r"{V}
        let a = rvStableSchurCharacter[[2]];
        let b = rvStableSchurCharacter[[1,1]];
        let y = add(a,neg(b));
        assert rvStableSchurAdams[2](v) == add(add(rvStableSchurCharacterScale[2](u),w),y);
        assert rvStableSchurAdams[2](rvStableSchurKronecker(y,v)) == rvStableSchurKronecker(rvStableSchurAdams[2](y),rvStableSchurAdams[2](v));
        assert rvStableSchurAdams[2](add(v,neg(y))) == add(rvStableSchurAdams[2](v),neg(rvStableSchurAdams[2](y)));
        assert rvStableSchurAdams[2](rvStableSchurAdams[3](y)) == rvStableSchurAdams[6](y);
        assert rvStableInnerSchur[[3]](v) == rvStableSchurSymmetric[3](v);
        assert rvStableInnerSchur[[1,1,1]](v) == rvStableSchurExterior[3](v);
        assert rvStableInnerSchur[[2,1,1]](neg(y)) == rvStableInnerSchur[[3,1]](y);
        assert rvStableInnerSchur[[2,1]](neg(y)) == neg(rvStableInnerSchur[[2,1]](y));
        assert rvStableSchurExterior[3](add(w,y)) == add(rvStableSchurExterior[3](w),add(rvStableSchurExterior[3](y),add(rvStableSchurKronecker(rvStableSchurExterior[2](w),y),rvStableSchurKronecker(w,rvStableSchurExterior[2](y)))));
        let h2 = rvStableSchurSymmetric[2](v);
        let h3 = rvStableSchurSymmetric[3](v);
        let h4 = rvStableSchurSymmetric[4](v);
        let h5 = rvStableSchurSymmetric[5](v);
        let term1 = rvStableSchurKronecker(h3,add(rvStableSchurKronecker(h2,v),neg(h3)));
        let term2 = rvStableSchurKronecker(h4,rvStableSchurKronecker(v,v));
        let term3 = rvStableSchurKronecker(h5,v);
        assert rvStableInnerSchur[[3,2,1]](v) == add(add(term1,neg(term2)),term3);
        let mixed = rvStableInnerSchur[[2,1]](a);
        rvStableSchurCharacterCoeff[[2,1]](mixed);
        rvStableSchurCharacterCoeff[[3,1]](mixed);
    "
    );
    let mut source = source;
    for m in 1..=4 {
        let inverse=(0..=m).map(|i|format!("rvStableSchurKronecker(rvStableSchurExterior[{i}](y),rvStableSchurExterior[{}](neg(y)))",m-i))
            .reduce(|a,b|format!("add({a},{b})")).unwrap();
        writeln!(source, "assert {inverse} == zero[RVStableSchurCharacter];").unwrap();
    }
    assert_eq!(values(&stdout(interpret(&source))), [9, 7]);
}

fn path_children(tail: &[usize], remove: bool) -> Vec<Vec<usize>> {
    let mut children = vec![tail.to_vec()];
    for (index, &row) in tail.iter().enumerate() {
        if remove && tail.get(index + 1).is_none_or(|&width| row > width) {
            let mut child = tail.to_vec();
            child[index] -= 1;
            if row == 1 {
                child.pop();
            }
            children.push(child);
        }
        if !remove && (index == 0 || tail[index - 1] > row) {
            let mut child = tail.to_vec();
            child[index] += 1;
            children.push(child);
        }
    }
    if !remove {
        let mut child = tail.to_vec();
        child.push(1);
        children.push(child);
    }
    children
}

fn hook_dimension(tail: &[usize], n: usize) -> i128 {
    let rows: Vec<_> = std::iter::once(n - tail.iter().sum::<usize>())
        .chain(tail.iter().copied())
        .collect();
    let hooks: usize = rows
        .iter()
        .enumerate()
        .flat_map(|(i, &width)| {
            let rows = &rows;
            (0..width)
                .map(move |j| width - j + rows[i + 1..].iter().filter(|&&row| row > j).count())
        })
        .product();
    ((1..=n).product::<usize>() / hooks) as i128
}

#[test]
fn deligne_place_permutations_recover_every_rank_three_inner_schur_coefficient() {
    let mut paths = vec![vec![Vec::new()]];
    for level in 1..=6usize {
        paths = paths
            .into_iter()
            .flat_map(|path| {
                path_children(path.last().unwrap(), !level.is_multiple_of(2))
                    .into_iter()
                    .map(move |child| {
                        let mut next = path.clone();
                        next.push(child);
                        next
                    })
            })
            .collect();
    }
    let mut groups = BTreeMap::<Vec<usize>, Vec<String>>::new();
    let mut source = format!(
        r"{V}
        let swap=rvPartitionLocal(rvPartitionDiagram[3]{{0,1,2,1,0,2}});
        let rotate=rvPartitionLocal(rvPartitionDiagram[3]{{0,1,2,1,2,0}});
        let sym=rvStableSchurSymmetric[3](v);
        let sch=rvStableInnerSchur[[2,1]](v);
        let alt=rvStableSchurExterior[3](v);
    "
    );
    for (i, path) in paths.iter().enumerate() {
        let literal = path
            .iter()
            .map(|tail| format!("{tail:?}"))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(source, "let p{i} = rvDelignePathProjector[3]{{{literal}}};").unwrap();
        groups
            .entry(path.last().unwrap().clone())
            .or_default()
            .push(format!("p{i}"));
    }
    for (i, (tail, paths)) in groups.iter().enumerate() {
        let isotype = paths
            .iter()
            .cloned()
            .reduce(|a, b| format!("add({a},{b})"))
            .unwrap();
        writeln!(source,"let t{i} = {isotype}; rvPartitionTraceAt[6](t{i}); rvPartitionTraceAt[6](rvPartitionCompose(t{i},swap)); rvPartitionTraceAt[6](rvPartitionCompose(t{i},rotate));").unwrap();
        for functor in ["sym", "sch", "alt"] {
            writeln!(source, "rvStableSchurCharacterCoeff[{tail:?}]({functor});").unwrap();
        }
    }
    let output = values(&stdout(interpret(&source)));
    assert_eq!(output.len(), groups.len() * 6);
    for ((shape, _), vals) in groups.iter().zip(output.chunks_exact(6)) {
        let dim = hook_dimension(shape, 6);
        assert!(vals[..3].iter().all(|value| value % dim == 0));
        let (a, b, c) = (vals[0] / dim, vals[1] / dim, vals[2] / dim);
        assert_eq!(
            [a + 3 * b + 2 * c, 2 * a - 2 * c, a - 3 * b + 2 * c],
            [6 * vals[3], 6 * vals[4], 6 * vals[5]]
        );
    }
}

#[test]
fn scalar_and_zero_cases_avoid_irrelevant_weight_work() {
    let source = format!(
        r"{V}
        let z=zero[RVStableSchurCharacter];
        let c=rvStableSchurCharacterScale[3](u);
        assert rvStableSchurAdams[{HUGE}](c) == c;
        assert rvStableSchurExterior[{HUGE}](u) == z;
        assert rvStableSchurSymmetric[{HUGE}](u) == u;
        assert rvStableSchurExterior[{HUGE}](neg(u)) == neg(u);
        assert rvStableSchurSymmetric[{HUGE}](neg(u)) == z;
        assert rvStableSchurExterior[0](z) == u;
        assert rvStableSchurExterior[{HUGE}](z) == z;
        assert rvStableInnerSchur[[]](z) == u;
        assert rvStableSchurExterior[{HUGE}](c) == z;
        assert rvStableSchurSymmetric[{HUGE}](neg(c)) == z;
        assert rvStableSchurExterior[2](c) == c;
        assert rvStableSchurSymmetric[3](c) == rvStableSchurCharacterScale[10](u);
        assert rvStableInnerSchur[[2,1]](c) == rvStableSchurCharacterScale[8](u);
        assert rvStableSchurSymmetric[1](v) == v;
        assert rvStableSchurAdams[1](rvStableSchurCharacter[[{HUGE}]]) == rvStableSchurCharacter[[{HUGE}]];
        rvSchurCharacterValue[[]](rvStableSchurSpecialize[0](u));
    "
    );
    assert_eq!(values(&stdout(interpret(&source))), [1]);
}

#[test]
fn bad_adams_indices_shapes_and_weight_bounds_have_direct_diagnostics() {
    for (query, expected) in [
        ("show rvStableSchurAdams[0](v);", "must be positive"),
        ("show rvStableInnerSchur[[1,2]](v);", "weakly decreasing"),
        (
            "show rvStableSchurAdams[3](S.hi);",
            "RVStableSchurCharacter",
        ),
        ("show rvStableSchurAdams[17](v);", "weighted cycle limit"),
        (
            "show rvStableInnerSchur[[9]](rvStableSchurCharacter[[2]]);",
            "weighted cycle limit",
        ),
        (
            "rvSchurCharacterValue[[3]](rvStableSchurSpecialize[2](v));",
            "RVSchurCharacter(3)",
        ),
        (
            "rvSchurCharacterValue[[17]](rvStableSchurSpecialize[17](v));",
            "weight limit of 16",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{V} {query}")));
        assert!(error.contains(expected), "{error}");
    }
}
