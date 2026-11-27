mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn low_rank_layers_distinguish_regular_radical_and_negligible_directions() {
    let output = stdout(example("deligne-singular-loewy"));
    assert!(output.starts_with(concat!(
        "rvPartitionLoewy[1,0] = dimension=2, negligible=2, radical=1, semisimple=1, tensorImage=0, powers=[2, 1, 0]\n",
        "rvPartitionLoewy[2,0] = dimension=15, negligible=15, radical=9, semisimple=6, tensorImage=0, powers=[15, 9, 5, 0]\n",
        "rvPartitionLoewy[2,1] = dimension=15, negligible=14, radical=3, semisimple=12, tensorImage=1, powers=[15, 3, 1, 0]\n",
        "rvPartitionLoewy[2,2] = dimension=15, negligible=7, radical=5, semisimple=10, tensorImage=8, powers=[15, 5, 1, 0]\n",
        "rvPartitionLoewy[2,3] = dimension=15, negligible=1, radical=0, semisimple=15, tensorImage=14, powers=[15, 0]\n",
        "rvPartitionLoewy[3,4] = dimension=203, negligible=16, radical=11, semisimple=192, tensorImage=187, powers=[203, 11] (through J^1)\n",
        "rvPartitionNegligibleContains = true\n",
        "rvPartitionRadicalContains[1] = false\n",
    )), "{output}");
    assert!(
        output.ends_with("rvPartitionRegularPair = 1\nrvPartitionCategoricalPair = 0\n"),
        "{output}"
    );
}

#[test]
fn rank_three_rational_nullspaces_and_both_quotient_dimensions_are_certified() {
    let expected = [
        (203, 111),
        (202, 44),
        (171, 93),
        (81, 39),
        (16, 11),
        (1, 0),
        (0, 0),
    ];
    let mut source = String::new();
    for (n, &(_, radical)) in expected.iter().enumerate() {
        writeln!(source, "rvPartitionLoewy[3,{n}];").unwrap();
        if radical != 0 {
            writeln!(
                source,
                "let a{n}=rvPartitionRadicalBasis[3,{n},1,{}];",
                radical - 1
            )
            .unwrap();
            writeln!(
                source,
                "assert rvPartitionSemisimple(a{n}) == zero[RVPartitionSemisimple(3,{n})];"
            )
            .unwrap();
            writeln!(
                source,
                "assert rvPartitionTensorImage(a{n}) == zero[RVPartitionTensorImage(3,{n})];"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&source));
    for (n, ((negligible, radical), line)) in expected.into_iter().zip(output.lines()).enumerate() {
        let prefix = format!("rvPartitionLoewy[3,{n}] = dimension=203, negligible={negligible}, radical={radical}, semisimple={}, tensorImage={}, powers=[203, {radical}]", 203-radical, 203-negligible);
        assert!(line.starts_with(&prefix), "{line}");
    }
    assert_eq!(output.lines().count(), 7);
}

#[test]
fn a_bad_modular_prime_cannot_supply_a_spurious_radical() {
    let output = stdout(interpret(
        "rvPartitionLoewy[0,0];\n\
         rvPartitionLoewy[1,1000000007];\n",
    ));
    assert_eq!(output, concat!(
        "rvPartitionLoewy[0,0] = dimension=1, negligible=0, radical=0, semisimple=1, tensorImage=1, powers=[1, 0]\n",
        "rvPartitionLoewy[1,1000000007] = dimension=2, negligible=0, radical=0, semisimple=2, tensorImage=2, powers=[2, 0]\n",
    ));
}

#[test]
fn every_nonzero_square_layer_has_actual_product_witnesses_and_cubes_vanish() {
    let mut oracle = String::new();
    let counts = [5, 1, 1];
    for (n, size) in counts.into_iter().enumerate() {
        for j in 0..size {
            writeln!(oracle, "rvPartitionLoewyWitness[2,{n},2,{j}];").unwrap();
        }
    }
    let witnesses = stdout(interpret(&oracle));
    let mut lines = witnesses.lines();
    let mut source = String::new();
    for (n, (size, radical)) in counts.into_iter().zip([9, 3, 5]).enumerate() {
        for j in 0..size {
            let report = lines.next().unwrap().split_once(" = [(1,").unwrap().1;
            let (left, right) = report.trim_end_matches(")]").split_once("), (1,").unwrap();
            writeln!(source, "let p{n}{j}=rvPartitionCompose(rvPartitionRadicalBasis[2,{n},1,{left}],rvPartitionRadicalBasis[2,{n},1,{right}]);").unwrap();
            writeln!(
                source,
                "assert p{n}{j} == rvPartitionRadicalBasis[2,{n},2,{j}];"
            )
            .unwrap();
            writeln!(
                source,
                "rvPartitionRadicalContains[2](p{n}{j});\nrvPartitionRadicalContains[3](p{n}{j});"
            )
            .unwrap();
            for k in 0..radical {
                writeln!(source, "assert rvPartitionCompose(p{n}{j},rvPartitionRadicalBasis[2,{n},1,{k}]) == zero[RVPartitionFiber(2,{n})];").unwrap();
            }
        }
    }
    assert!(lines.next().is_none());
    let output = stdout(interpret(&source));
    assert_eq!(
        output,
        "rvPartitionRadicalContains[2] = true\nrvPartitionRadicalContains[3] = false\n".repeat(7)
    );
}

fn labels_equivalent(diagram: &[usize; 4], tuple: &[usize; 4]) -> bool {
    (0..4).all(|i| (0..i).all(|j| diagram[i] != diagram[j] || tuple[i] == tuple[j]))
}

#[test]
fn quotient_products_match_explicit_tensor_matrices_and_rational_arithmetic() {
    let mut input = String::new();
    let cases = [
        ([0, 0, 1, 2], [0, 1, 1, 0]),
        ([0, 1, 2, 2], [0, 0, 0, 0]),
        ([0, 1, 0, 2], [0, 0, 1, 2]),
    ];
    let orbits = [
        [0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, 0, 1, 1],
        [0, 1, 0, 0],
        [0, 1, 0, 1],
        [0, 1, 1, 0],
        [0, 1, 1, 1],
    ];
    let mut expected = Vec::new();
    for (i, (a, b)) in cases.into_iter().enumerate() {
        writeln!(
            input,
            "let a{i} = rvPartitionSpecialize[2](rvPartitionDiagram[2]{{{},{},{},{}}});",
            a[0], a[1], a[2], a[3]
        )
        .unwrap();
        writeln!(
            input,
            "let b{i} = rvPartitionSpecialize[2](rvPartitionDiagram[2]{{{},{},{},{}}});",
            b[0], b[1], b[2], b[3]
        )
        .unwrap();
        writeln!(input, "let ab{i}=rvPartitionCompose(a{i},b{i});").unwrap();
        writeln!(input, "let im{i}=rvPartitionCompose(rvPartitionTensorImage(a{i}),rvPartitionTensorImage(b{i}));").unwrap();
        writeln!(input, "assert im{i} == rvPartitionTensorImage(ab{i});").unwrap();
        writeln!(input, "assert rvPartitionCompose(rvPartitionSemisimple(a{i}),rvPartitionSemisimple(b{i})) == rvPartitionSemisimple(ab{i});").unwrap();
        for [a0, a1, b0, b1] in orbits {
            writeln!(
                input,
                "rvPartitionFiberCoefficient[[{a0},{a1},{b0},{b1}]](im{i});"
            )
            .unwrap();
            expected.push(
                (0..4)
                    .filter(|m| {
                        let (m0, m1) = (m / 2, m % 2);
                        labels_equivalent(&a, &[a0, a1, m0, m1])
                            && labels_equivalent(&b, &[m0, m1, b0, b1])
                    })
                    .count(),
            );
        }
    }
    input.push_str("let half = rvPartitionFiberScale[1,2](im0);\nassert add(half,half) == im0;\nassert add(half,neg(half)) == zero[RVPartitionTensorImage(2,2)];\n");
    let output = stdout(interpret(&input));
    assert_eq!(
        output
            .lines()
            .map(|line| line.split_once(" = ").unwrap().1.parse::<usize>().unwrap())
            .collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn categorical_trace_form_counts_exact_orbits_by_falling_factorials() {
    let output = stdout(interpret(
        "let x=rvPartitionSpecialize[3](rvPartitionOrbit[2]{0,1,0,1});\n\
         let y=rvPartitionSpecialize[3](rvPartitionOrbit[2]{0,1,2,0});\n\
         let z=rvPartitionNegligibleBasis[2,3,0];\n\
         rvPartitionCategoricalPair(x,x);\n\
         rvPartitionCategoricalPair(y,y);\n\
         rvPartitionCategoricalPair(x,y);\n\
         rvPartitionCategoricalPair(z,z);\n\
         rvPartitionNegligibleContains(z);\n\
         rvPartitionRadicalContains[1](z);\n\
         rvPartitionLoewy[2,4];\n",
    ));
    assert_eq!(output, concat!(
        "rvPartitionCategoricalPair = 6\n",
        "rvPartitionCategoricalPair = 6\n",
        "rvPartitionCategoricalPair = 0\n",
        "rvPartitionCategoricalPair = 0\n",
        "rvPartitionNegligibleContains = true\n",
        "rvPartitionRadicalContains[1] = false\n",
        "rvPartitionLoewy[2,4] = dimension=15, negligible=0, radical=0, semisimple=15, tensorImage=15, powers=[15, 0]\n",
    ));
}

#[test]
fn poles_actual_finite_types_and_out_of_range_queries_are_checked() {
    assert_eq!(
        stdout(interpret(
            "let d=rvPartitionDiagram[1]{0,1};\n\
         let rem=rvPartitionLocalDivide[0,1,0,0](rvPartitionLocal(rvPartitionScale[1,1](d)));\n\
         assert rvPartitionSpecialize[0](rem) == rvPartitionSpecialize[0](d);\n\
         let half=rvPartitionLocalDivide[1,0,0,0](rvPartitionLocal(d));\n\
         rvPartitionFiberCoefficient[[0,1]](rvPartitionSpecialize[0](half));\n\
         let du=rvPartitionScale[2,18446744073709551615](d);\n\
         rvPartitionFiberCoefficient[[0,1]](rvPartitionSpecialize[1](du));\n",
        )),
        "rvPartitionFiberCoefficient[[0, 1]] = 1/2\nrvPartitionFiberCoefficient[[0, 1]] = 2\n"
    );
    for (source, expected) in [
        ("rvPartitionLoewy[4,3];", "limited to rank 3"),
        ("let u=rvPartitionSpecialize[0](rvDeligneSquareProjector[0]);", "pole at t=0"),
        ("let u=rvPartitionFiberScale[1,0](zero[RVPartitionFiber(2,2)]);", "denominator must be positive"),
        ("let x=rvPartitionCompose(zero[RVPartitionFiber(1,0)],zero[RVPartitionFiber(1,1)]);", "same rank, parameter, and quotient"),
        ("let x=rvPartitionCompose(zero[RVPartitionFiber(1,0)],zero[RVPartitionSemisimple(1,0)]);", "same rank, parameter, and quotient"),
        ("let x=rvPartitionTensorImage(zero[RVPartitionSemisimple(1,0)]);", "requires an RVPartitionFiber"),
        ("let x=rvPartitionRadicalBasis[2,3,1,0];", "radical basis index is out of range"),
        ("let x=rvPartitionNegligibleBasis[2,4,0];", "negligible basis index is out of range"),
        ("rvPartitionLoewyWitness[3,2,2,0];", "higher Loewy powers"),
    ] {
        assert!(diagnostic(interpret(source)).contains(expected), "{source}");
    }
}
