mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn two_middle_summands_give_actual_nonzero_cancelling_composites() {
    let output = stdout(example("deligne-auslander-reiten"));
    let mut lines = output.lines();
    assert_eq!(lines.next().unwrap(), "rvPartitionIndecomposables[2,0] = modules=[(0, 1, 1, false), (0, 2, 3, false), (0, 3, 4, true), (1, 1, 1, true), (2, 1, 2, false), (2, 2, 3, false), (2, 3, 5, true)]");
    assert!(lines
        .next()
        .unwrap()
        .contains("powers=[29, 22, 14, 6, 2, 0]"));
    let mesh = lines.next().unwrap();
    assert!(
        mesh.starts_with("rvPartitionAR[2,0,2,2] = translate=[0, 2], middle=[[2, 3], [0, 1]],"),
        "{mesh}"
    );
    assert!(mesh.contains("incomingFactors="));
    assert!(mesh.contains("outgoingFactors="));
    assert_eq!(lines.next().unwrap(), "rvPartitionHomApply[0] = [(1,1)]");
    assert_eq!(lines.next().unwrap(), "rvPartitionHomApply[0] = [(1,-1)]");
    assert_eq!(
        lines.next().unwrap(),
        "rvPartitionAuslanderRadicalContains[2] = true"
    );
    assert!(lines.next().is_none());
}

#[test]
fn all_auslander_powers_and_mesh_quotient_layers_match_small_fixtures() {
    let output = stdout(interpret(concat!(
        "rvPartitionAuslander[2,0];rvPartitionAuslander[2,1];rvPartitionAuslander[2,2];",
        "rvPartitionAuslander[2,3];rvPartitionAuslander[1,0];rvPartitionAuslander[0,0];",
    )));
    let expected = [
        (7, 29, "[29, 22, 14, 6, 2, 0]", 8, 4, "[7, 8, 8, 4, 2, 0]"),
        (7, 20, "[20, 13, 7, 3, 1, 0]", 6, 3, "[7, 6, 4, 2, 1, 0]"),
        (7, 20, "[20, 13, 7, 3, 1, 0]", 6, 3, "[7, 6, 4, 2, 1, 0]"),
        (4, 4, "[4, 0]", 0, 0, "[4, 0]"),
        (2, 5, "[5, 3, 1, 0]", 2, 1, "[2, 2, 1, 0]"),
        (1, 1, "[1, 0]", 0, 0, "[1, 0]"),
    ];
    assert_eq!(output.lines().count(), expected.len());
    for (line, (modules, dimension, powers, arrows, meshes, quotient)) in
        output.lines().zip(expected)
    {
        assert!(
            line.contains(&format!(
                "indecomposables={modules}, dimension={dimension}, powers={powers}"
            )),
            "{line}"
        );
        assert!(
            line.ends_with(&format!("meshes={meshes}, pathQuotient={quotient}")),
            "{line}"
        );
        let arrow_data = line
            .split_once("arrows=")
            .unwrap()
            .1
            .split_once(", meshes=")
            .unwrap()
            .0;
        assert_eq!(arrow_data.matches('(').count(), arrows, "{line}");
    }
}

#[test]
fn the_two_middle_sequence_is_nonsplit_and_actual_paths_enforce_truncation() {
    let output = stdout(interpret(concat!(
        "let f1=rvPartitionARMap[2,0,2,2,0,0];\n",
        "let g1=rvPartitionARMap[2,0,2,2,1,0];\n",
        "let p=rvPartitionHomCompose(g1,f1);\n",
        "let sec=rvPartitionHomBasis[2,0,2,2,2,3,2];\n",
        "let ret=rvPartitionHomBasis[2,0,2,3,0,2,1];\n",
        "assert rvPartitionHomCompose(g1,sec)==zero[RVPartitionHom(2,0,2,2,2,2)];\n",
        "assert rvPartitionHomCompose(ret,f1)==zero[RVPartitionHom(2,0,0,2,0,2)];\n",
        "rvPartitionHomDim[2,0,2,2,0,1];\n",
        "rvPartitionHomDim[2,0,0,1,0,2];\n",
        "assert rvPartitionHomPath[1](p)==zero[RVPartitionFiber(2,0)];\n",
        "assert rvPartitionCompose(rvPartitionModulePath[2,0,0,2,1],rvPartitionHomPath[0](p))==rvPartitionModulePath[2,0,2,3,2];\n",
        "rvPartitionAR[2,0,2,3];\n",
    )));
    assert_eq!(
        output,
        concat!(
            "rvPartitionHomDim[2,0,2,2,0,1] = dimension=0, layers=[]\n",
            "rvPartitionHomDim[2,0,0,1,0,2] = dimension=0, layers=[]\n",
            "rvPartitionAR[2,0,2,3] = projective\n",
        )
    );
}

#[test]
fn independently_enumerated_hom_layers_and_cartan_entries_agree() {
    let fixtures = [
        (0, vec![3, 1, 3], vec![2, 1, 0]),
        (1, vec![1, 1, 3, 2], vec![0, 1, 3, 2]),
        (2, vec![1, 3, 2, 1], vec![0, 2, 1, 3]),
    ];
    let mut source = String::new();
    let mut expected = Vec::new();
    for (n, lengths, successors) in fixtures {
        let modules: Vec<_> = lengths
            .iter()
            .enumerate()
            .flat_map(|(i, &l)| (1..=l).map(move |m| (i, m)))
            .collect();
        for &(i, l) in &modules {
            for &(j, m) in &modules {
                writeln!(source, "rvPartitionHomDim[2,{n},{i},{l},{j},{m}];").unwrap();
                let layers: Vec<_> = (0..m)
                    .filter(|&s| {
                        let head = (0..s).fold(j, |v, _| successors[v]);
                        head == i && s + l >= m
                    })
                    .collect();
                expected.push(format!(
                    "rvPartitionHomDim[2,{n},{i},{l},{j},{m}] = dimension={}, layers={layers:?}",
                    layers.len()
                ));
            }
        }
    }
    assert_eq!(
        stdout(interpret(&source)).lines().collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn hom_arithmetic_and_actual_radical_basis_positions_are_reusable() {
    let output = stdout(interpret(concat!(
        "let e=rvPartitionHomBasis[2,0,0,3,0,3,0];\n",
        "let nil=rvPartitionHomBasis[2,0,0,3,0,3,2];\n",
        "let f=add(e,rvPartitionHomScale[1,2](nil));\n",
        "assert rvPartitionHomCompose(f,f)==add(e,nil);\n",
        "assert add(f,neg(f))==zero[RVPartitionHom(2,0,0,3,0,3)];\n",
        "rvPartitionAuslanderRadicalContains[1](e);\n",
        "rvPartitionAuslanderRadicalContains[1](nil);\n",
        "let a=rvPartitionAuslanderArrow[1,0,0];\n",
        "rvPartitionAuslanderRadicalContains[1](a);\n",
        "rvPartitionAuslanderRadicalContains[2](a);\n",
        "let d=rvPartitionAuslanderRadicalBasis[1,0,2,0];\n",
        "rvPartitionAuslanderRadicalContains[2](d);\n",
        "rvPartitionAuslanderRadicalContains[3](d);\n",
    )));
    assert_eq!(
        output,
        concat!(
            "rvPartitionAuslanderRadicalContains[1] = false\n",
            "rvPartitionAuslanderRadicalContains[1] = true\n",
            "rvPartitionAuslanderRadicalContains[1] = true\n",
            "rvPartitionAuslanderRadicalContains[2] = false\n",
            "rvPartitionAuslanderRadicalContains[2] = true\n",
            "rvPartitionAuslanderRadicalContains[3] = false\n",
        )
    );
}

#[test]
fn bad_module_lengths_grades_and_mesh_indices_are_checked() {
    for (source,expected) in [
        ("rvPartitionAuslander[3,0];","limited to rank 2"),
        ("rvPartitionHomDim[2,0,0,4,0,1];","module top or length is out of range"),
        ("show rvPartitionHomBasis[2,0,0,1,0,3,0];","not in this Hom basis"),
        ("show rvPartitionHomCompose(rvPartitionAuslanderArrow[2,1,0],rvPartitionAuslanderArrow[2,0,0]);","matching intermediate module"),
        ("show rvPartitionARMap[2,0,2,3,0,0];","projective module has no almost split sequence"),
        ("show rvPartitionARMap[2,0,2,1,0,1];","middle index is out of range"),
        ("show rvPartitionARMap[2,0,2,1,3,0];","side must be 0 or 1"),
        ("rvPartitionHomApply[1](rvPartitionHomBasis[2,0,0,1,0,1,0]);","module basis layer is out of range"),
        ("show rvPartitionAuslanderRadicalBasis[1,0,9,0];","radical basis index is out of range"),
    ] {
        let message=diagnostic(interpret(source));
        assert!(message.contains(expected),"{source}: {message}");
    }
}
