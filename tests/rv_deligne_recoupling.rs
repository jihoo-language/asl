mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn runnable_example_reports_block_and_normalized_squared_overlaps() {
    let output = stdout(example("deligne-racah-blocks"));
    assert_eq!(output, concat!(
        "rvDeligneBlockPaths[3, [2, 1]] = [[[], [], [1], [1], [1, 1], [1, 1], [2, 1]], [[], [], [1], [1], [2], [2], [2, 1]]]\n",
        "rvDeligneBlock[3, [2, 1]] = [[1/2, 3/4], [1, -1/2]]\n",
        "rvDeligneBlockAt[3, [2, 1], 6] = [[0, 1], [0, 0]]\n",
        "rvDeligneRacah[3, [2, 1]] = [[1/4, 3/4], [3/4, 1/4]]\n",
        "rvPartitionTrace = numerator=[(0,1)], denominator=(t-1)^2\n",
        "rvPartitionTraceAt[6] = 1/25\n",
    ));
}

#[test]
fn exact_block_entries_multiply_and_reconstruct_a_noncentral_corner() {
    let mut input = String::from(
        "let a = rvPartitionLocal(rvPartitionDiagram[3]{0,1,2,1,0,2});\n\
         let b = rvPartitionLocal(rvPartitionDiagram[3]{0,1,2,0,2,1});\n\
         let ab = rvPartitionCompose(a,b);\n\
         let one = rvPartitionLocal(rvPartitionIdentity[0]);\n\
         let zero0 = zero[RVPartitionLocal(0)];\n\
         let zero3 = zero[RVPartitionLocal(3)];\n",
    );
    for i in 0..2 {
        for j in 0..2 {
            writeln!(input, "let u{i}{j} = rvDeligneMatrixUnit[3,[2,1],{i},{j}];").unwrap();
            for value in ["a", "b", "ab"] {
                writeln!(
                    input,
                    "let {value}{i}{j} = rvDeligneBlockEntry[3,[2,1],{i},{j}]({value});"
                )
                .unwrap();
            }
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            writeln!(input, "assert ab{i}{j} == add(rvPartitionTensor(a{i}0,b0{j}),rvPartitionTensor(a{i}1,b1{j}));").unwrap();
        }
        for (j, k) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            for l in 0..2 {
                let result = if j == k {
                    format!("u{i}{l}")
                } else {
                    "zero3".into()
                };
                writeln!(
                    input,
                    "assert rvPartitionCompose(u{i}{j},u{k}{l}) == {result};"
                )
                .unwrap();
            }
        }
    }
    input.push_str(
        "assert rvDeligneBlockEntry[3,[2,1],0,1](u01) == one;\n\
         assert rvDeligneBlockEntry[3,[2,1],1,0](u01) == zero0;\n\
         assert rvPartitionCompose(u01,rvDeligneMatrixUnit[3,[3],0,0]) == zero3;\n\
         let central = add(u00,u11);\n\
         let first = add(rvPartitionTensor(ab00,u00),rvPartitionTensor(ab01,u01));\n\
         let last = add(rvPartitionTensor(ab10,u10),rvPartitionTensor(ab11,u11));\n\
         assert add(first,last) == rvPartitionCompose(ab,central);\n\
         rvDeligneBlockAt[3,[2,1],6](ab);\n",
    );
    assert_eq!(
        stdout(interpret(&input)),
        "rvDeligneBlockAt[3, [2, 1], 6] = [[-1/2, -3/4], [1, -1/2]]\n"
    );
}

#[test]
fn all_rank_two_blocks_reconstruct_an_arbitrary_inhomogeneous_diagram() {
    let mut input = String::from(
        "let top = rvPartitionDiagram[2]{0,0,0,1};\n\
         let bottom = rvPartitionOrbit[2]{0,1,2,2};\n\
         let d = rvPartitionLocal(add(top,rvPartitionScale[-2,3](bottom)));\n",
    );
    let mut terms = Vec::new();
    for (tail, size) in [("[]", 2), ("[1]", 3), ("[2]", 1), ("[1,1]", 1)] {
        for i in 0..size {
            for j in 0..size {
                let term = format!("rvPartitionTensor(rvDeligneBlockEntry[2,{tail},{i},{j}](d),rvDeligneMatrixUnit[2,{tail},{i},{j}])");
                terms.push(term);
            }
        }
    }
    let reconstruction = terms
        .into_iter()
        .reduce(|a, b| format!("add({a},{b})"))
        .unwrap();
    writeln!(input, "assert {reconstruction} == d;").unwrap();
    assert_eq!(stdout(interpret(&input)), "");
}

#[test]
fn every_rank_three_bridge_and_overlap_normalization_is_exercised() {
    let mut input = String::new();
    let sizes = [
        ("[]", 5),
        ("[1]", 10),
        ("[2]", 6),
        ("[1,1]", 6),
        ("[2,1]", 2),
        ("[3]", 1),
        ("[1,1,1]", 1),
    ];
    for (tail, _) in sizes {
        writeln!(
            input,
            "rvDeligneBlockAt[3,{tail},6](rvPartitionIdentity[3]);"
        )
        .unwrap();
        writeln!(input, "rvDeligneRacah[3,{tail}];").unwrap();
    }
    let output = stdout(interpret(&input));
    let mut lines = output.lines().map(|line| line.split_once(" = ").unwrap().1);
    for (_, size) in sizes {
        let rows: Vec<_> = (0..size)
            .map(|i| {
                let entries: Vec<_> = (0..size).map(|j| if i == j { "1" } else { "0" }).collect();
                format!("[{}]", entries.join(", "))
            })
            .collect();
        assert_eq!(lines.next().unwrap(), format!("[{}]", rows.join(", ")));
        assert!(lines.next().is_some());
    }
    assert!(lines.next().is_none());
}

#[test]
fn trivial_pair_is_the_contraction_of_two_standard_factors() {
    let output = stdout(interpret(
        "let w = rvDelignePathProjector[1]{[],[],[1]};\n\
         let g = rvDeligneSquareProjector[3];\n\
         assert rvDeligneMatrixUnit[3,[1],3,3] == rvPartitionTensor(g,w);\n\
         let c = rvDeligneRacahEntry[3,[1],3,3];\n\
         rvPartitionTraceAt[4](c);\n\
         rvPartitionTraceAt[5](c);\n\
         rvPartitionTraceAt[6](c);\n\
         rvDeligneRacahAt[0,[],0];\n",
    ));
    assert_eq!(
        output,
        concat!(
            "rvPartitionTraceAt[4] = 1/9\n",
            "rvPartitionTraceAt[5] = 1/16\n",
            "rvPartitionTraceAt[6] = 1/25\n",
            "rvDeligneRacahAt[0, [], 0] = [[1]]\n",
        )
    );
}

#[test]
fn rank_shape_index_input_and_actual_specialization_are_checked() {
    for (source, expected) in [
        ("rvDeligneRacah[4,[1]];", "limited to rank 3"),
        ("rvDeligneRacah[3,[1,2]];", "invalid Deligne terminal tail"),
        ("rvDeligneRacah[3,[2,2]];", "invalid Deligne terminal tail"),
        (
            "let a = rvDeligneMatrixUnit[3,[2,1],0,2];",
            "path index is out of range",
        ),
        (
            "let a = rvDeligneBlockEntry[3,[2,1],0,0](rvPartitionIdentity[2]);",
            "diagram of the same rank",
        ),
        ("rvDeligneRacahAt[3,[3],5];", "terminal has dimension zero"),
        ("rvDeligneRacahAt[3,[1],2];", "pole at t=2"),
        (
            "rvDeligneBlockAt[2,[1],5](rvPartitionIdentity[3]);",
            "diagram of the same rank",
        ),
        (
            "rvPartitionTraceAt[1](rvDeligneRacahEntry[3,[1],3,3]);",
            "pole at t=1",
        ),
    ] {
        assert!(diagnostic(interpret(source)).contains(expected), "{source}");
    }
}
