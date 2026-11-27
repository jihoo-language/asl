mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn hilbert_square_and_cube_match_classical_branching() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-reduced-llv")),
        concat!(
            "rvHilbReducedLLVDecomposition[2, 2] = components=[([1], 1, 23)] total=23\n",
            "rvHilbReducedLLVDecomposition[2, 4] = components=[([2], 1, 275), ([], 1, 1)] total=276\n",
            "rvHilbPrimitiveLLVDecomposition[2, 2] = components=[([1], 1, 22)] total=22\n",
            "rvHilbPrimitiveLLVDecomposition[2, 4] = components=[([2], 1, 252), ([], 1, 1)] total=253\n",
            "rvHilbReducedLLVDecomposition[3, 4] = components=[([2], 1, 275), ([1], 1, 23), ([], 1, 1)] total=299\n",
            "rvHilbPrimitiveLLVDecomposition[3, 6] = components=[([3], 1, 2002), ([1, 1], 1, 231), ([1], 1, 22)] total=2255\n",
            "rvHilbReducedLLVMultiplicity[3, 4][1] = 1\n",
            "rvHilbPrimitiveLLVMultiplicity[3, 6][1, 1] = 1\n",
        )
    );
    // In degree six on the Hilbert cube, the harmonic cubic of the
    // standard 25 leaves Sym^3(23), and the exterior square leaves
    // exterior^2(23) + 1. Branch to SO(22) and subtract degree four.
    let result = stdout(interpret(
        "rvHilbReducedLLVDecomposition[3,6]; rvHilbPrimitiveLLVDecomposition[3,4];",
    ));
    assert_eq!(
        result,
        concat!(
            "rvHilbReducedLLVDecomposition[3, 6] = components=[([3], 1, 2277), ([1, 1], 1, 253), ([1], 1, 23), ([], 1, 1)] total=2554\n",
            "rvHilbPrimitiveLLVDecomposition[3, 4] = components=[([2], 1, 252), ([1], 1, 22), ([], 2, 1)] total=276\n",
        )
    );
}

fn total(line: &str) -> u128 {
    line.rsplit_once(" total=").unwrap().1.parse().unwrap()
}

fn hodge_through_degree(line: &str, n: usize) -> Vec<u128> {
    let matrix = line
        .split_once("hodge=[[")
        .unwrap()
        .1
        .strip_suffix("]]")
        .unwrap();
    let mut grades = vec![0; 2 * n + 1];
    for (p, row) in matrix.split("], [").enumerate() {
        for (q, value) in row.split(", ").enumerate() {
            if p + q <= 2 * n {
                grades[p + q] += value.parse::<u128>().unwrap();
            }
        }
    }
    grades
}

#[test]
fn every_grade_matches_goettsche_and_the_polarized_hodge_difference() {
    for n in 2..=9 {
        let mut program = format!("rvHilbLLVHodgeTable[{n}];\n");
        for d in 0..=4 * n {
            writeln!(program, "rvKummerHilbBetti[{n},{d}];").unwrap();
            writeln!(program, "rvHilbReducedLLVDecomposition[{n},{d}];").unwrap();
            writeln!(program, "rvHilbPrimitiveLLVDecomposition[{n},{d}];").unwrap();
        }
        let output = stdout(interpret(&program));
        let mut lines = output.lines();
        let hodge = hodge_through_degree(lines.next().unwrap(), n);
        for d in 0..=4 * n {
            let betti = lines.next().unwrap().split_once(" = ").unwrap().1;
            assert_eq!(
                total(lines.next().unwrap()),
                betti.parse::<u128>().unwrap(),
                "n={n}, d={d}"
            );
            let primitive = if d <= 2 * n {
                hodge[d] - d.checked_sub(2).map_or(0, |lower| hodge[lower])
            } else {
                0
            };
            assert_eq!(total(lines.next().unwrap()), primitive, "n={n}, d={d}");
        }
        assert_eq!(lines.next(), None);
    }
}

#[test]
fn coefficient_reports_handle_empty_grades_weights_and_limits() {
    let got = stdout(interpret(
        r"
        rvHilbReducedLLVDecomposition[0,0];
        rvHilbPrimitiveLLVDecomposition[0,0];
        rvHilbReducedLLVDecomposition[3,-2];
        rvHilbReducedLLVDecomposition[3,14];
        rvHilbReducedLLVDecomposition[1000000000000000000,0];
        rvHilbPrimitiveLLVDecomposition[1000000000000000000,9];
        rvHilbPrimitiveLLVDecomposition[3,8];
        rvHilbReducedLLVMultiplicity[3,4]{1,0,0};
        rvHilbReducedLLVMultiplicity[3,4]{3};
        rvHilbPrimitiveLLVMultiplicity[3,4]{};
        rvHilbPrimitiveLLVMultiplicity[3,4]{3};
        rvHilbPrimitiveLLVMultiplicity[3,6]{1,1,1,1,1,1,1,1,1,1,-1};
        ",
    ));
    assert_eq!(
        got,
        concat!(
            "rvHilbReducedLLVDecomposition[0, 0] = components=[([], 1, 1)] total=1\n",
            "rvHilbPrimitiveLLVDecomposition[0, 0] = components=[([], 1, 1)] total=1\n",
            "rvHilbReducedLLVDecomposition[3, -2] = components=[] total=0\n",
            "rvHilbReducedLLVDecomposition[3, 14] = components=[] total=0\n",
            "rvHilbReducedLLVDecomposition[1000000000000000000, 0] = components=[([], 1, 1)] total=1\n",
            "rvHilbPrimitiveLLVDecomposition[1000000000000000000, 9] = components=[] total=0\n",
            "rvHilbPrimitiveLLVDecomposition[3, 8] = components=[] total=0\n",
            "rvHilbReducedLLVMultiplicity[3, 4][1] = 1\n",
            "rvHilbReducedLLVMultiplicity[3, 4][3] = 0\n",
            "rvHilbPrimitiveLLVMultiplicity[3, 4][] = 2\n",
            "rvHilbPrimitiveLLVMultiplicity[3, 4][3] = 0\n",
            "rvHilbPrimitiveLLVMultiplicity[3, 6][1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1] = 0\n",
        )
    );
    for (program, message) in [
        (
            "rvHilbReducedLLVDecomposition[1,2];",
            "whose LLV type is D12",
        ),
        (
            "rvHilbPrimitiveLLVDecomposition[1,4];",
            "whose LLV type is D12",
        ),
        (
            "rvHilbReducedLLVMultiplicity[3,4]{1,2};",
            "type B11 highest weight",
        ),
        (
            "rvHilbPrimitiveLLVMultiplicity[3,4]{1,-1};",
            "type D11 highest weight",
        ),
        (
            "rvHilbPrimitiveLLVMultiplicity[3,4]{1,1,1,1,1,1,1,1,1,1,-2};",
            "type D11 highest weight",
        ),
        (
            "rvHilbReducedLLVDecomposition[1000000000000000000,2];",
            "polynomial computation limit",
        ),
    ] {
        assert!(
            diagnostic(interpret(program)).contains(message),
            "{program}"
        );
    }
}

#[test]
fn both_chiral_constituents_occur_at_their_first_allowed_length() {
    let got = stdout(interpret(
        r"
        rvHilbReducedLLVMultiplicity[66,132]{1,1,1,1,1,1,1,1,1,1,1};
        rvHilbPrimitiveLLVMultiplicity[66,132]{1,1,1,1,1,1,1,1,1,1,1};
        rvHilbPrimitiveLLVMultiplicity[66,132]{1,1,1,1,1,1,1,1,1,1,-1};
        rvHilbPrimitiveLLVMultiplicity[66,130]{1,1,1,1,1,1,1,1,1,1,-1};
        ",
    ));
    assert_eq!(
        got,
        concat!(
            "rvHilbReducedLLVMultiplicity[66, 132][1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] = 1\n",
            "rvHilbPrimitiveLLVMultiplicity[66, 132][1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] = 1\n",
            "rvHilbPrimitiveLLVMultiplicity[66, 132][1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1] = 1\n",
            "rvHilbPrimitiveLLVMultiplicity[66, 130][1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1] = 0\n",
        )
    );
}
