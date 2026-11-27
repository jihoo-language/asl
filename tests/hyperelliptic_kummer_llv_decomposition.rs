mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

type Components = BTreeMap<Vec<usize>, (usize, u128)>;

fn decomposition(line: &str) -> (Components, u128, Vec<u128>) {
    let (_, rest) = line.split_once("components=[").unwrap();
    let (components, rest) = rest.split_once("] total=").unwrap();
    let (total, betti) = rest.split_once(" betti=[").unwrap();
    let entries = components
        .split("), (")
        .map(|entry| {
            let entry = entry.trim_start_matches('(').trim_end_matches(')');
            let (shape, sizes) = entry.split_once("], ").unwrap();
            let shape = shape.strip_prefix('[').unwrap();
            let shape = if shape.is_empty() {
                Vec::new()
            } else {
                shape
                    .split(", ")
                    .map(|part| part.parse().unwrap())
                    .collect()
            };
            let (multiplicity, dimension) = sizes.split_once(", ").unwrap();
            (
                shape,
                (multiplicity.parse().unwrap(), dimension.parse().unwrap()),
            )
        })
        .collect();
    let betti = betti
        .strip_suffix(']')
        .unwrap()
        .split(", ")
        .map(|x| x.parse().unwrap())
        .collect();
    (entries, total.parse().unwrap(), betti)
}

#[test]
fn example_identifies_the_first_harmonic_and_exterior_representations() {
    let output = stdout(example("hyperelliptic-kummer-llv-decomposition"));
    let mut lines = output.lines();
    assert_eq!(lines.next(), Some("rvHilbLLVMultiplicity[2][2] = 1"));
    assert_eq!(lines.next(), Some("rvHilbLLVDimension[2] = 324"));
    let (cube, cube_total, cube_grades) = decomposition(lines.next().unwrap());
    assert_eq!(
        cube,
        BTreeMap::from([(vec![3], (1, 2900)), (vec![1, 1], (1, 300))])
    );
    assert_eq!(cube_total, 3200);
    assert_eq!(cube_grades, [1, 23, 299, 2554, 299, 23, 1]);
    let (four, four_total, _) = decomposition(lines.next().unwrap());
    assert_eq!(
        four,
        BTreeMap::from([
            (vec![4], (1, 20150)),
            (vec![2, 1], (1, 5175)),
            (vec![2], (1, 324)),
            (vec![], (1, 1))
        ])
    );
    assert_eq!(four_total, 25650);
    assert_eq!(lines.next(), None);
}

#[test]
fn high_decompositions_reproduce_every_independent_goettsche_betti_number() {
    let mut source = String::new();
    for n in 4..=9 {
        writeln!(source, "rvHilbLLVDecomposition[{n}];").unwrap();
        for i in 0..=2 * n {
            writeln!(source, "rvKummerHilbBetti[{n},{}];", 2 * i).unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    for (n, (count, expected_total)) in [
        (4, 25650),
        (6, 176256),
        (11, 1073720),
        (13, 5930496),
        (20, 30178575),
        (25, 143184000),
    ]
    .into_iter()
    .enumerate()
    {
        let n = n + 4;
        let (entries, total, betti) = decomposition(lines.next().unwrap());
        assert_eq!(entries.len(), count, "n={n}");
        assert_eq!(total, expected_total, "n={n}");
        assert_eq!(betti.len(), 2 * n + 1);
        assert_eq!(betti.iter().copied().sum::<u128>(), total);
        for (index, actual) in betti.into_iter().enumerate() {
            let expected: u128 = lines
                .next()
                .unwrap()
                .split_once(" = ")
                .unwrap()
                .1
                .parse()
                .unwrap();
            assert_eq!(actual, expected, "n={n}, i={}", 2 * index);
        }
    }
    assert_eq!(lines.next(), None);
}

#[test]
fn complete_published_tables_through_seven_have_the_exact_multiplicities() {
    type Table = &'static [(&'static [usize], usize)];
    let tables: &[(usize, Table)] = &[
        (2, &[(&[2], 1)]),
        (
            5,
            &[
                (&[5], 1),
                (&[3, 1], 1),
                (&[3], 1),
                (&[2, 1], 1),
                (&[1, 1], 1),
                (&[1], 1),
            ],
        ),
        (
            6,
            &[
                (&[6], 1),
                (&[4, 1], 1),
                (&[4], 1),
                (&[3, 1], 1),
                (&[3], 1),
                (&[2, 2], 1),
                (&[2, 1], 1),
                (&[2], 2),
                (&[1, 1, 1], 1),
                (&[1], 1),
                (&[], 1),
            ],
        ),
        (
            7,
            &[
                (&[7], 1),
                (&[5, 1], 1),
                (&[5], 1),
                (&[4, 1], 1),
                (&[4], 1),
                (&[3, 2], 1),
                (&[3, 1], 2),
                (&[3], 2),
                (&[2, 1, 1], 1),
                (&[2, 1], 2),
                (&[2], 1),
                (&[1, 1], 2),
                (&[1], 2),
            ],
        ),
    ];
    let mut source = String::new();
    for (n, _) in tables {
        writeln!(source, "rvHilbLLVDecomposition[{n}];").unwrap();
    }
    let output = stdout(interpret(&source));
    for ((n, table), line) in tables.iter().zip(output.lines()) {
        let (entries, total, betti) = decomposition(line);
        if *n == 2 {
            assert_eq!(total, 324);
            assert_eq!(betti, [1, 23, 276, 23, 1]);
        }
        let actual = entries
            .into_iter()
            .map(|(lambda, (multiplicity, _))| (lambda, multiplicity))
            .collect::<BTreeMap<_, _>>();
        let expected = table
            .iter()
            .map(|&(lambda, multiplicity)| (lambda.to_vec(), multiplicity))
            .collect();
        assert_eq!(actual, expected, "n={n}");
    }
}

#[test]
fn individual_coefficients_and_nominal_types_are_exact() {
    let out = stdout(interpret(
        r"
        rvHilbLLVMultiplicity[0]{};
        rvHilbLLVMultiplicity[0]{1};
        rvHilbLLVMultiplicity[3]{2};
        rvHilbLLVMultiplicity[6]{2,0,0};
        rvHilbLLVMultiplicity[7]{3,1};
        rvHilbLLVMultiplicity[3]{184467440737095516160};
        rvHilbLLVDimension{};
        rvHilbLLVDimension{1};
        rvHilbLLVDimension{1,1};
        rvHilbLLVDimension{1,1,1,1,1,1,1,1,1,1,1,1};
        rvHilbLLVDecomposition[0];
    ",
    ));
    let lines = out.lines().collect::<Vec<_>>();
    assert_eq!(
        &lines[..10],
        &[
            "rvHilbLLVMultiplicity[0][] = 1",
            "rvHilbLLVMultiplicity[0][1] = 0",
            "rvHilbLLVMultiplicity[3][2] = 0",
            "rvHilbLLVMultiplicity[6][2] = 2",
            "rvHilbLLVMultiplicity[7][3, 1] = 2",
            "rvHilbLLVMultiplicity[3][184467440737095516160] = 0",
            "rvHilbLLVDimension[] = 1",
            "rvHilbLLVDimension[1] = 25",
            "rvHilbLLVDimension[1, 1] = 300",
            "rvHilbLLVDimension[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] = 5200300",
        ]
    );
    assert_eq!(
        decomposition(lines[10]),
        (BTreeMap::from([(vec![], (1, 1))]), 1, vec![1])
    );
}

#[test]
fn weyl_dimension_uses_exact_wide_cancellation_and_decimal_output() {
    assert_eq!(
        stdout(interpret("rvHilbLLVDimension{500};")),
        "rvHilbLLVDimension[500] = 15532816129550473797602715976873765373376\n"
    );
}

#[test]
fn rejects_nondominant_or_excessive_queries() {
    for (source, expected) in [
        ("rvHilbLLVMultiplicity[1]{1};", "whose LLV type is D12"),
        ("rvHilbLLVDecomposition[1];", "whose LLV type is D12"),
        (
            "rvHilbLLVDimension{1,2};",
            "highest weight must be weakly decreasing",
        ),
        (
            "rvHilbLLVMultiplicity[3]{1,1,1,1,1,1,1,1,1,1,1,1,1};",
            "at most twelve nonzero coordinates",
        ),
        (
            "rvHilbLLVMultiplicity[1000000000000000000]{};",
            "exceeds the polynomial computation limit",
        ),
        (
            "rvHilbLLVDimension{340282366920938463463374607431768211455};",
            "exceeds the exact integer range",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{error}");
    }
}
