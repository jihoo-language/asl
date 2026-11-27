mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn runnable_example_shows_nonzero_gauge_classes_and_baer_rank_addition() {
    let output = stdout(example("deligne-extension-groupoid"));
    assert!(output.starts_with("rvPartitionExtensionGroupoid[E] = q=2, arrows=[(0, 0, 4)], vertexSplittings=[4], Z1=4, B1=0, Ext1=4, Hom=4, classes=16, stabilizer=16, mass=1\n"),"{output}");
    assert!(output.contains("total=16, checkedIsotypes=3"));
    assert!(output.contains("([(0, 1, 4)], classes=1, hall=35, aut=20160, mass=1/16)"));
    assert!(output.contains("([(0, 1, 2), (0, 2, 1)], classes=9, hall=3"));
    assert!(output.contains("([(0, 2, 2)], classes=6, hall=1"));
    assert!(
        output.ends_with(concat!(
            "rvPartitionExtensionLift[E] = [1, 0, 0, 1]\n",
            "rvPartitionExtensionCoboundary[G] = cocycle=[1, 0], class=[]\n",
            "rvPartitionHallHom[q=2] = dimension=2, cardinality=4\n",
            "rvPartitionHallAut[q=2] = endDimension=2, radicalDimension=1, cardinality=2\n",
        )),
        "{output}"
    );
}

#[test]
fn all_matrix_ranks_match_independent_hall_flags_over_four_fields() {
    let mut source = String::from("rvunram F[2,2]{1,1,1};\n");
    let samples = [
        ("2", 2, [1, 9, 6], [35, 3, 1]),
        ("3", 3, [1, 32, 48], [130, 4, 1]),
        ("5", 5, [1, 144, 480], [806, 6, 1]),
        ("F", 4, [1, 75, 180], [357, 5, 1]),
    ];
    for (i, (field, _, _, _)) in samples.into_iter().enumerate() {
        writeln!(
            source,
            "let m{i}=rvPartitionHallModule[1,0,{field}]{{(0,1,2)}};"
        )
        .unwrap();
        writeln!(
            source,
            "rvpartitionextensions E{i}(m{i},m{i});\nrvPartitionExtensionDistribution[E{i}];"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    for (line, (_, q, ext, hall)) in output.lines().zip(samples) {
        assert!(
            line.contains(&format!("q={q}, total={}", ext.iter().sum::<usize>())),
            "{line}"
        );
        for ((shape, ext), hall) in ["[(0, 1, 4)]", "[(0, 1, 2), (0, 2, 1)]", "[(0, 2, 2)]"]
            .into_iter()
            .zip(ext)
            .zip(hall)
        {
            assert!(
                line.contains(&format!("({shape}, classes={ext}, hall={hall},")),
                "{line}"
            );
        }
    }
    assert_eq!(output.lines().count(), 4);
}

#[test]
fn nonzero_changes_of_splitting_preserve_the_canonical_nontrivial_class() {
    let output = stdout(interpret(concat!(
        "let M=rvPartitionHallModule[1,0,3]{(0,1,1),(0,2,1)};\n",
        "let N=rvPartitionHallModule[1,0,3]{(0,1,1)};\n",
        "rvpartitionextensions G(M,N);\n",
        "rvPartitionExtensionGroupoid[G];\n",
        "let x=rvPartitionExtensionCocycle[G]{1,0,0};\n",
        "let y=rvPartitionExtensionCocycle[G]{1,2,0};\n",
        "assert x==y;\n",
        "assert x==rvPartitionExtension[G]{1};\n",
        "assert rvPartitionExtensionMiddle(x)==rvPartitionHallModule[1,0,3]{(0,2,2)};\n",
        "rvPartitionExtensionLift(y);\n",
        "rvPartitionExtensionCoboundary[G]{0,0,1};\n",
        "rvPartitionExtensionDistribution[G];\n",
    )));
    assert!(
        output.contains("Z1=2, B1=1, Ext1=1, Hom=2, classes=3, stabilizer=9, mass=1/3\n"),
        "{output}"
    );
    assert!(output.contains("rvPartitionExtensionLift[G] = [1, 0, 0]\n"));
    assert!(output.contains("rvPartitionExtensionCoboundary[G] = cocycle=[0, 2, 0], class=[0]\n"));
    assert!(output.contains("classes=2, hall=4"), "{output}");
}

#[test]
fn nonhereditary_hom_minus_ext_is_not_fixed_by_dimension_vectors() {
    let output = stdout(interpret(concat!(
        "let M=rvPartitionHallModule[1,0,3]{(0,1,2)};\n",
        "let J=rvPartitionHallModule[1,0,3]{(0,2,1)};\n",
        "rvpartitionextensions E(M,M);\n",
        "rvpartitionextensions F(J,M);\n",
        "rvPartitionExtensionGroupoid[E];\n",
        "rvPartitionExtensionGroupoid[F];\n",
        "rvPartitionHallAut(M);\n",
        "rvPartitionHallAut(J);\n",
        "rvPartitionExtensionDistribution[F];\n",
    )));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].contains("Ext1=4, Hom=4, classes=81"));
    assert!(lines[1].contains("Ext1=0, Hom=2, classes=1"));
    assert_eq!(
        lines[2],
        "rvPartitionHallAut[q=3] = endDimension=4, radicalDimension=0, cardinality=48"
    );
    assert_eq!(
        lines[3],
        "rvPartitionHallAut[q=3] = endDimension=2, radicalDimension=1, cardinality=6"
    );
    assert!(lines[4].contains("q=3, total=1"));
}

#[test]
fn linked_block_orientation_and_named_field_addition_are_genuine() {
    let output=stdout(interpret(concat!(
        "rvunram F[2,1]{1,1,1};\n",
        "let S=rvPartitionHallModule[1,0,F]{(0,1,1)};\n",
        "rvpartitionextensions E(S,S);\n",
        "assert add(rvPartitionExtension[E]{2},rvPartitionExtension[E]{3})==rvPartitionExtension[E]{1};\n",
        "card \\x : RVPartitionExtension(E) => S.hi;\n",
        "let M=rvPartitionHallModule[2,0,5]{(0,2,1)};\n",
        "let N=rvPartitionHallModule[2,0,5]{(2,1,1)};\n",
        "rvpartitionextensions G(M,N);\n",
        "rvpartitionextensions H(N,M);\n",
        "rvPartitionExtensionDistribution[G];\n",
        "rvPartitionExtensionDistribution[H];\n",
    )));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines[0], "card = 4 / 4");
    assert!(
        lines[1].contains("total=1") && lines[1].contains("classes=1, hall=5"),
        "{output}"
    );
    assert!(
        lines[2].contains("total=5") && lines[2].contains("classes=4, hall=1"),
        "{output}"
    );
}

#[test]
fn invalid_cocycles_and_nominal_space_mismatches_are_rejected() {
    let fixture="let S=rvPartitionHallModule[1,0,3]{(0,1,1)}; let J=rvPartitionHallModule[1,0,3]{(0,2,1)}; rvpartitionextensions E(S,S); rvpartitionextensions G(J,S);";
    for (suffix,expected) in [
        ("show rvPartitionExtensionCocycle[G]{0,1};","is not a cocycle"),
        ("show rvPartitionExtension[E]{3};","outside the finite field codes"),
        ("show rvPartitionExtension[E]{};","requires 1 coordinates"),
        ("rvPartitionExtensionCoboundary[G]{1};","requires 2 coordinates"),
        ("rvpartitionextensions E(J,S);","already defined"),
        ("rvPartitionExtensionDistribution[Unknown];","unknown partition extension space"),
        ("rvpartitionextensions H(S,S); show add(rvPartitionExtension[E]{1},rvPartitionExtension[H]{1});","RVPartitionExtension"),
        ("rvpartitionextensions H(S,rvPartitionHallModule[1,0,5]{(0,1,1)});","same rank, parameter, and coefficient field"),
        ("let B=rvPartitionHallModule[1,0,3]{(0,1,4)}; rvpartitionextensions H(B,B);","limited to total length 6"),
        ("let B=rvPartitionHallModule[1,0,17]{(0,1,2)}; rvpartitionextensions H(B,B); rvPartitionExtensionDistribution[H];","exceeds 50,000"),
    ] {
        let message=diagnostic(interpret(&format!("{fixture} {suffix}")));
        assert!(message.contains(expected),"{suffix}: {message}");
    }
}
