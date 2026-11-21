mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cusp_jets_with_the_same_residue_have_different_jordan_blocks() {
    assert_eq!(
        stdout(example("rank-jets")),
        concat!(
            "rvdim[One, Cusp, [0, 1, 2], 7, [0, 0, 1]] = 4\n",
            "rvjet[One, Cusp, [0, 1, 2], 7, 7, [[0, 1], [], [1]]] = dimension 4, blocks {1: 4}\n",
            "rvjet[One, Cusp, [0, 1, 2], 7, 7, [[], [0, 1], [1]]] = dimension 6, blocks {1: 2, 2: 2}\n",
            "rvjet[One, Cusp, [0, 1, 2], 7, 7, [[0, 0, 0, 1], [0, 0, 1], [1]]] = dimension 18, blocks {2: 2, 7: 2}\n",
            "rvjet[Cusp, Cusp, [0, 1, 2], 7, 4, [[], [0, 1], [1]]] = dimension 20, blocks {1: 12, 2: 4}\n",
            "rvjet[Cusp, Cusp, [0, 1, 2], 7, 4, [[0, 0, 0, 1], [0, 0, 1], [1]]] = dimension 40, blocks {2: 12, 4: 4}\n",
        )
    );
}

#[test]
fn residue_points_and_unit_rescaling_preserve_the_expected_homology() {
    let cusp = include_str!("../examples/rank-jets.asl")
        .split("\n// All three")
        .next()
        .unwrap();
    let program = format!(
        "{cusp}\n{}",
        r"
        rvdim[Cusp, Cusp, [0,1,2], 7, [1,3,1]];
        rvjet[Cusp, Cusp, [0,1,2], 7, 1, [[1],[3],[1,0,0]]];
        rvjet[Cusp, Cusp, [0,1,2], 7, 1, [[2],[3],[1]]];
        rvjet[Cusp, Cusp, [0,1,2], 7, 4, [[],[0,2,1],[2,1]]];
        "
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvdim[Cusp, Cusp, [0, 1, 2], 7, [1, 3, 1]] = 4\n",
            "rvjet[Cusp, Cusp, [0, 1, 2], 7, 1, [[1], [3], [1, 0, 0]]] = dimension 4, blocks {1: 4}\n",
            "rvjet[Cusp, Cusp, [0, 1, 2], 7, 1, [[2], [3], [1]]] = dimension 0, blocks {}\n",
            "rvjet[Cusp, Cusp, [0, 1, 2], 7, 4, [[], [0, 2, 1], [2, 1]]] = dimension 20, blocks {1: 12, 2: 4}\n",
        )
    );
}

#[test]
fn zero_actions_and_shifts_in_large_endpoint_degrees_are_exact() {
    let max = usize::MAX;
    let program = format!(
        "qmodule T[0] {{a:{max}, b:{}}} {{0(a)=>[b]}};\n\
         rvjet[T, One, [0], 283, 6, [[0,255]]];\n\
         rvjet[T, One, [0], 3, 3, [[]]];\n\
         rvjet[One, One, [], 3, {max}, []];",
        max - 1
    );
    assert_eq!(
        stdout(interpret(&program)),
        format!(
            "rvjet[T, One, [0], 283, 6, [[0, 255]]] = dimension 2, blocks {{1: 2}}\n\
             rvjet[T, One, [0], 3, 3, [[]]] = dimension 6, blocks {{3: 2}}\n\
             rvjet[One, One, [], 3, {max}, []] = dimension {max}, blocks {{{max}: 1}}\n"
        )
    );
}

#[test]
fn invalid_coefficients_fields_truncations_and_oversized_matrices_are_rejected() {
    for (program, expected) in [
        (
            "rvjet[One, One, [0], 7, 0, [[1]]];",
            "positive truncation order",
        ),
        (
            "rvjet[One, One, [0], 7, 2, []];",
            "one jet per selected primitive",
        ),
        (
            "rvjet[One, One, [0], 7, 2, [[4]]];",
            "field coefficient 4 is out of range",
        ),
        (
            "rvjet[One, One, [0], 7, 2, [[1,0,1]]];",
            "nonzero coefficient beyond the truncation",
        ),
        ("rvjet[One, One, [0], 9, 2, [[1]]];", "reducible over F2"),
        (
            "qmodule T[0] {a:1,b:0} {0(a)=>[b]}; rvjet[T,One,[0],3,300000,[[1]]];",
            "Margolis jet calculation exceeds the storage limit",
        ),
        (
            "qmodule U[0] {a:0} {}; rvjet[U, One, [1], 7, 2, [[]]];",
            "not declared for this module",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
    let error = diagnostic(interpret(&format!(
        "qmodule T[0] {{a:1,b:0}} {{0(a)=>[b]}}; rvjet[T,One,[0],3,{},[[]]];",
        usize::MAX
    )));
    assert!(
        error.contains("Margolis jet dimension is out of range"),
        "{error}"
    );
}
