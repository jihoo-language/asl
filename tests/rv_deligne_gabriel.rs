mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn example_has_nonzero_roundtrip_and_required_cubic_relations() {
    let output = stdout(example("deligne-gabriel-quiver"));
    let mut lines = output.lines();
    let data = lines.next().unwrap();
    assert!(data.contains("simples=[1, 1, 2], primitives=[[0], [1], [2, 3]], projectives=[4, 1, 5], basicDimension=7"), "{data}");
    assert!(data.contains("ext1=[[0, 0, 1], [0, 0, 0], [1, 0, 0]], ext2=[[0, 0, 1], [0, 0, 0], [1, 0, 0]], arrows=[(2, 0), (0, 2)]"), "{data}");
    assert_eq!(lines.next().unwrap(), "rvPartitionGabrielRelations[2,0] = quadraticPaths=[[0, 1], [1, 0]], products=[[(0,1)], [(0,-1)]], quadraticRelations=[], zeroCubicPaths=[[0, 1, 0], [1, 0, 1]], minimalCubicRelations=[[0, 1, 0], [1, 0, 1]]");
    assert_eq!(
        lines.next().unwrap(),
        "rvPartitionRadicalContains[2] = true"
    );
    assert_eq!(lines.next().unwrap(), "rvPartitionGabrielRelations[2,2] = quadraticPaths=[[0, 1], [1, 0]], products=[[], [(0,1)]], quadraticRelations=[[(0,1)]], zeroCubicPaths=[[0, 1, 0], [1, 0, 1]], minimalCubicRelations=[]");
    assert!(lines.next().is_none());
}

#[test]
fn all_primitive_idempotents_split_the_full_specialized_identity() {
    let mut source = String::new();
    for (n, count) in [4, 6, 6, 7, 7].into_iter().enumerate() {
        writeln!(source, "let zero{n}=zero[RVPartitionFiber(2,{n})];").unwrap();
        for i in 0..count {
            writeln!(source, "let e{n}{i}=rvPartitionPrimitive[2,{n},{i}];").unwrap();
            writeln!(
                source,
                "assert rvPartitionCompose(e{n}{i},e{n}{i}) == e{n}{i};"
            )
            .unwrap();
            for j in 0..i {
                writeln!(source, "assert rvPartitionCompose(e{n}{i},e{n}{j}) == zero{n};\nassert rvPartitionCompose(e{n}{j},e{n}{i}) == zero{n};").unwrap();
            }
        }
        let complete = (0..count)
            .map(|i| format!("e{n}{i}"))
            .reduce(|a, b| format!("add({a},{b})"))
            .unwrap();
        writeln!(
            source,
            "assert {complete} == rvPartitionSpecialize[{n}](rvPartitionIdentity[2]);"
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn arrows_are_actual_corner_classes_and_products_recover_the_square_basis() {
    let mut source = String::new();
    for (n, (i, j)) in [(2, 0), (3, 2), (2, 1)].into_iter().enumerate() {
        writeln!(source, "let e{n}=rvPartitionBasicIdempotent[2,{n},{i}];\nlet f{n}=rvPartitionBasicIdempotent[2,{n},{j}];").unwrap();
        writeln!(
            source,
            "let a{n}=rvPartitionArrow[2,{n},0];\nlet b{n}=rvPartitionArrow[2,{n},1];"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvPartitionCompose(f{n},rvPartitionCompose(a{n},e{n})) == a{n};"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvPartitionCompose(e{n},rvPartitionCompose(b{n},f{n})) == b{n};"
        )
        .unwrap();
        writeln!(source, "let aba{n}=rvPartitionCompose(a{n},b{n});\nlet bab{n}=rvPartitionCompose(b{n},a{n});\nlet zero{n}=zero[RVPartitionFiber(2,{n})];").unwrap();
        writeln!(source, "assert rvPartitionCompose(a{n},bab{n}) == zero{n};\nassert rvPartitionCompose(b{n},aba{n}) == zero{n};").unwrap();
        if n == 0 {
            writeln!(
                source,
                "assert bab{n} == rvPartitionBasicSquare[2,{n},{i},{i},0];"
            )
            .unwrap();
            writeln!(
                source,
                "assert aba{n} == neg(rvPartitionBasicSquare[2,{n},{j},{j},0]);"
            )
            .unwrap();
        } else {
            writeln!(source, "assert bab{n} == zero{n};\nassert aba{n} == rvPartitionBasicSquare[2,{n},{j},{j},0];").unwrap();
        }
        writeln!(
            source,
            "rvPartitionRadicalContains[1](a{n});\nrvPartitionRadicalContains[2](a{n});"
        )
        .unwrap();
    }
    assert_eq!(
        stdout(interpret(&source)),
        "rvPartitionRadicalContains[1] = true\nrvPartitionRadicalContains[2] = false\n".repeat(3)
    );
}

#[test]
fn cartan_and_projective_layers_distinguish_all_three_singular_rank_two_cases() {
    let output = stdout(interpret("rvPartitionGabriel[2,1]; rvPartitionGabriel[2,2]; rvPartitionGabriel[2,3]; rvPartitionGabriel[2,257];"));
    let lines: Vec<_> = output.lines().collect();
    for (line, (simple, projective, basic)) in lines.iter().zip([
        ("[1, 3, 1, 1]", "[1, 3, 3, 2]", 7),
        ("[1, 1, 2, 2]", "[1, 4, 3, 2]", 7),
        ("[1, 1, 3, 2]", "[1, 1, 3, 2]", 4),
        ("[1, 1, 3, 2]", "[1, 1, 3, 2]", 4),
    ]) {
        assert!(line.contains(&format!("simples={simple},")), "{line}");
        assert!(
            line.contains(&format!("projectives={projective}, basicDimension={basic}")),
            "{line}"
        );
    }
    assert_eq!(lines.len(), 4);
    assert!(lines[0].contains("cartan=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 1], [0, 0, 1, 1]]"));
    assert!(lines[1].contains("cartan=[[1, 0, 0, 0], [0, 2, 1, 0], [0, 1, 1, 0], [0, 0, 0, 1]]"));
    assert!(lines[1].contains("ext2=[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]]"));
    assert!(lines[2].ends_with("arrows=[]"));
}

#[test]
fn the_dual_number_fiber_has_one_loop_and_one_minimal_quadratic_relation() {
    let output = stdout(interpret("rvPartitionGabriel[0,0]; rvPartitionGabriel[1,0]; rvPartitionGabrielRelations[1,0]; rvPartitionGabrielRelations[1,1];"));
    assert_eq!(output, concat!(
        "rvPartitionGabriel[0,0] = simples=[1], primitives=[[0]], projectives=[1], basicDimension=1, cartan=[[1]], layers=[[[1, 0, 0]]], ext1=[[0]], ext2=[[0]], arrows=[]\n",
        "rvPartitionGabriel[1,0] = simples=[1], primitives=[[0]], projectives=[2], basicDimension=2, cartan=[[2]], layers=[[[1, 1, 0]]], ext1=[[1]], ext2=[[1]], arrows=[(0, 0)]\n",
        "rvPartitionGabrielRelations[1,0] = quadraticPaths=[[0, 0]], products=[[]], quadraticRelations=[[(0,1)]], zeroCubicPaths=[[0, 0, 0]], minimalCubicRelations=[]\n",
        "rvPartitionGabrielRelations[1,1] = quadraticPaths=[], products=[], quadraticRelations=[], zeroCubicPaths=[], minimalCubicRelations=[]\n",
    ));
}

#[test]
fn geometry_and_type_boundaries_produce_direct_diagnostics() {
    for (source, expected) in [
        ("rvPartitionGabriel[3,0];", "limited to rank 2"),
        (
            "show rvPartitionPrimitive[2,0,4];",
            "generator index is out of range",
        ),
        (
            "show rvPartitionArrow[2,3,0];",
            "generator index is out of range",
        ),
        (
            "show rvPartitionBasicSquare[2,2,2,2,0];",
            "generator index is out of range",
        ),
        (
            "show rvPartitionCompose(rvPartitionArrow[2,1,0],rvPartitionArrow[2,2,1]);",
            "same rank, parameter, and quotient",
        ),
    ] {
        let output = diagnostic(interpret(source));
        assert!(output.contains(expected), "{source}: {output}");
    }
}
