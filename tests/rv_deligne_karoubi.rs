mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const ANCHORS: [usize; 7] = [0, 1, 1, 0, 1, 5, 6];

#[test]
fn runnable_example_splits_the_square_with_typed_corner_isomorphisms() {
    assert_eq!(
        stdout(example("deligne-karoubi-square")),
        concat!(
            "rvKaroubiDimension[A] = numerator=[(0,1)], denominator=1\n",
            "rvKaroubiDimension[H] = numerator=[(0,-1), (1,1)], denominator=1\n",
            "rvKaroubiDimension[K] = numerator=[(1,-3), (2,1)], denominator=2\n",
            "rvKaroubiDimension[L] = numerator=[(0,2), (1,-3), (2,1)], denominator=2\n",
            "rvKaroubiTrace = numerator=[(0,-1), (1,1)], denominator=1\n",
            "rvPartitionZeroAt[3] = true\n",
            "rvPartitionTraceAt[4] = 2\n",
            "rvPartitionTraceAt[5] = 5\n",
        )
    );
}

#[test]
fn normalization_cancels_exactly_and_specialization_keeps_rational_signs() {
    let script = r"
        let u = rvPartitionLocal(rvPartitionIdentity[0]);
        let t = rvPartitionScale[1,1](u);
        let tm1 = add(t,neg(u));
        let tm2 = add(t,rvPartitionScale[-2,0](u));
        let d = rvPartitionCompose(rvPartitionScale[2,1](u),rvPartitionCompose(tm1,tm2));
        let q = rvPartitionLocalDivide[1,1,1,1](d);
        assert q == u;
        let remainder = rvPartitionLocalDivide[0,0,0,1](add(rvPartitionScale[1,2](u),rvPartitionScale[-4,0](u)));
        assert remainder == add(t,rvPartitionScale[2,0](u));
        let a = rvPartitionLocalDivide[0,0,1,0](u);
        let b = rvPartitionLocalDivide[0,0,0,1](u);
        let difference = add(a,neg(b));
        assert difference == rvPartitionLocalDivide[0,0,1,1](neg(u));
        assert rvPartitionCompose(a,tm1) == u;
        assert rvPartitionCompose(difference,rvPartitionCompose(tm1,tm2)) == neg(u);
        assert add(a,neg(a)) == zero[RVPartitionLocal(0)];
        assert rvPartitionTensor(q,q) == u;
        rvPartitionEntry[1,[],[]](q);
        rvPartitionEntry[2,[],[]](remainder);
        rvPartitionEntry[0,[],[]](a);
        rvPartitionEntry[1,[],[]](rvPartitionLocalDivide[1,0,0,0](b));
        rvPartitionEntry[1,[],[]](add(a,neg(a)));
        let j = rvPartitionLocal(rvPartitionDiagram[1]{0,1});
        let jOverT = rvPartitionLocalDivide[0,1,0,0](j);
        rvPartitionTraceAt[0](jOverT);
    ";
    assert_eq!(
        stdout(interpret(script)),
        concat!(
            "rvPartitionEntry[1] = 1\n",
            "rvPartitionEntry[2] = 4\n",
            "rvPartitionEntry[0] = -1\n",
            "rvPartitionEntry[1] = -1/2\n",
            "rvPartitionEntry[1] = 0\n",
            "rvPartitionTraceAt[0] = 1\n",
        )
    );
}

#[test]
fn square_projectors_preserve_specializations_where_denominator_factors_cancel() {
    let output = stdout(interpret(
        r"
        rvPartitionZeroAt[0](rvDeligneSquareProjector[5]);
        rvPartitionZeroAt[1](rvDeligneSquareProjector[4]);
        rvPartitionZeroAt[1](rvDeligneSquareProjector[6]);
        rvPartitionZeroAt[2](rvDeligneSquareProjector[6]);
        rvPartitionZeroAt[1](rvDeligneSquareProjector[0]);
    ",
    ));
    assert_eq!(
        output,
        concat!(
            "rvPartitionZeroAt[0] = true\n",
            "rvPartitionZeroAt[1] = true\n",
            "rvPartitionZeroAt[1] = true\n",
            "rvPartitionZeroAt[2] = true\n",
            "rvPartitionZeroAt[1] = false\n",
        )
    );
}

#[test]
fn all_fifteen_matrix_units_obey_the_symbolic_product_table() {
    let mut source = String::from("let z = zero[RVPartitionLocal(2)];\n");
    let mut units = Vec::new();
    for (i, &a) in ANCHORS.iter().enumerate() {
        for (j, &b) in ANCHORS.iter().enumerate() {
            if a == b {
                units.push((i, j));
                writeln!(source, "let e{i}x{j} = rvDeligneSquareUnit[{i},{j}];").unwrap();
            } else {
                writeln!(source, "assert rvDeligneSquareUnit[{i},{j}] == z;").unwrap();
            }
        }
    }
    assert_eq!(units.len(), 15);
    for &(i, j) in &units {
        for &(k, l) in &units {
            let right = if j == k {
                format!("e{i}x{l}")
            } else {
                "z".into()
            };
            writeln!(
                source,
                "assert rvPartitionCompose(e{i}x{j},e{k}x{l}) == {right};"
            )
            .unwrap();
        }
    }
    let total = (0..7)
        .map(|i| format!("e{i}x{i}"))
        .reduce(|a, b| format!("add({a},{b})"))
        .unwrap();
    writeln!(
        source,
        "assert {total} == rvPartitionLocal(rvPartitionIdentity[2]);"
    )
    .unwrap();
    assert!(stdout(interpret(&source)).is_empty());
}

fn rational_mod(line: &str, p: i128) -> i128 {
    let value = line.rsplit_once(" = ").unwrap().1;
    let (a, b) = value.split_once('/').unwrap_or((value, "1"));
    let a = a.parse::<i128>().unwrap().rem_euclid(p);
    let b = b.parse::<i128>().unwrap().rem_euclid(p);
    assert_ne!(b, 0);
    a * power_mod(b, p - 2, p) % p
}

fn power_mod(mut value: i128, mut exponent: i128, p: i128) -> i128 {
    let mut out = 1;
    while exponent != 0 {
        if exponent % 2 != 0 {
            out = out * value % p;
        }
        value = value * value % p;
        exponent /= 2;
    }
    out
}

fn rank(mut matrix: Vec<Vec<i128>>, p: i128) -> usize {
    let mut rank = 0;
    for j in 0..matrix.len() {
        let Some(row) = (rank..matrix.len()).find(|&i| matrix[i][j] != 0) else {
            continue;
        };
        matrix.swap(rank, row);
        let inverse = power_mod(matrix[rank][j], p - 2, p);
        let pivot = matrix[rank].clone();
        for row in matrix.iter_mut().skip(rank + 1) {
            let multiplier = row[j] * inverse % p;
            for (old, &top) in row.iter_mut().zip(&pivot).skip(j) {
                *old = (*old - multiplier * top).rem_euclid(p);
            }
        }
        rank += 1;
    }
    rank
}

#[test]
fn actual_finite_matrix_ranks_and_traces_match_the_seven_summands() {
    for n in [3, 4, 5] {
        let tuples: Vec<_> = (0..n).flat_map(|a| (0..n).map(move |b| [a, b])).collect();
        let degree = tuples.len();
        let mut source = String::new();
        for i in 0..7 {
            writeln!(source, "let p{i} = rvDeligneSquareProjector[{i}];").unwrap();
            writeln!(source, "rvPartitionTraceAt[{n}](p{i});").unwrap();
            for upper in &tuples {
                for lower in &tuples {
                    writeln!(source, "rvPartitionEntry[{n},{upper:?},{lower:?}](p{i});").unwrap();
                }
            }
        }
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        let w = n - 1;
        for expected in [1, w, w, 1, w, n * (n - 3) / 2, (n - 1) * (n - 2) / 2] {
            let trace: usize = lines
                .next()
                .unwrap()
                .rsplit_once(" = ")
                .unwrap()
                .1
                .parse()
                .unwrap();
            assert_eq!(trace, expected);
            let p = 1009;
            let matrix: Vec<Vec<_>> = (0..degree)
                .map(|_| {
                    (0..degree)
                        .map(|_| {
                            let line = lines.next().unwrap();
                            if expected == 0 {
                                assert_eq!(line.rsplit_once(" = ").unwrap().1, "0");
                            }
                            rational_mod(line, p)
                        })
                        .collect()
                })
                .collect();
            assert_eq!(rank(matrix, p), expected);
        }
        assert!(lines.next().is_none());
    }
}

#[test]
fn nominal_corner_arithmetic_and_underlying_maps_are_exact() {
    let output = stdout(interpret(
        r"
        let a = rvDeligneSquareProjector[0];
        let b = rvDeligneSquareProjector[3];
        rvkaroubi A(a);
        rvkaroubi B(b);
        rvkaroubi Z(zero[RVPartitionLocal(2)]);
        let ab = rvKaroubiMap[A,B](rvDeligneSquareUnit[3,0]);
        let ba = rvKaroubiMap[B,A](rvDeligneSquareUnit[0,3]);
        assert rvKaroubiCompose(rvKaroubiIdentity[B],ab) == ab;
        assert rvKaroubiCompose(ab,rvKaroubiIdentity[A]) == ab;
        assert rvKaroubiCompose(ba,ab) == rvKaroubiIdentity[A];
        assert rvKaroubiCompose(ab,ba) == rvKaroubiIdentity[B];
        assert add(ab,neg(ab)) == zero[RVKaroubiMap(A,B)];
        assert rvKaroubiCompose(ab,zero[RVKaroubiMap(Z,A)]) == zero[RVKaroubiMap(Z,B)];
        assert rvKaroubiIdentity[Z] == zero[RVKaroubiMap(Z,Z)];
        assert rvKaroubiUnderlying(rvKaroubiIdentity[A]) == a;
        rvKaroubiTrace(rvKaroubiIdentity[Z]);
    ",
    ));
    assert_eq!(output, "rvKaroubiTrace = numerator=[], denominator=1\n");
}

#[test]
fn nonidempotents_wrong_corners_nominality_and_true_poles_are_checked() {
    let setup = r"
        let unit = rvPartitionLocal(rvPartitionIdentity[0]);
        let a = rvDeligneSquareProjector[0];
        let g = rvDeligneSquareProjector[3];
        rvkaroubi A(a);
        rvkaroubi Alias(a);
        rvkaroubi G(g);
        rvkaroubi Root(unit);
    ";
    for (query, expected) in [
        (
            "rvkaroubi X(rvPartitionLocal(rvPartitionDiagram[1]{0,1}));",
            "requires an idempotent",
        ),
        ("rvkaroubi A(a);", "already defined"),
        ("show rvKaroubiIdentity[Missing];", "unknown Karoubi object"),
        ("show zero[RVKaroubiMap(A,Root)];", "same partition rank"),
        ("show rvKaroubiMap[A,G](a);", "f = e_target f e_source"),
        (
            "show rvKaroubiCompose(rvKaroubiIdentity[A],rvKaroubiIdentity[Alias]);",
            "same intermediate object",
        ),
        (
            "rvKaroubiTrace(rvKaroubiMap[A,G](rvDeligneSquareUnit[3,0]));",
            "endomorphism",
        ),
        (
            "show rvPartitionCompose(a,rvPartitionIdentity[2]);",
            "explicit rvPartitionLocal",
        ),
        ("show rvDeligneSquareProjector[7];", "index must be in 0..7"),
        (
            "show rvPartitionLocalDivide[65,0,0,0](unit);",
            "exceeds 64 factors",
        ),
        ("show rvPartitionScale[1,129](unit);", "degree 128"),
        (
            "rvPartitionEntry[1,[],[]](rvPartitionLocalDivide[0,0,1,0](unit));",
            "pole at t=1",
        ),
        ("rvPartitionZeroAt[0](a);", "pole at t=0"),
        (
            "rvPartitionZeroAt[2](rvDeligneSquareProjector[4]);",
            "pole at t=2",
        ),
        (
            "rvPartitionZeroAt[1](rvDeligneSquareProjector[5]);",
            "pole at t=1",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {query}")));
        assert!(error.contains(expected), "{error}");
    }
}
