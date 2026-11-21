mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn hochster_faces_carry_augmented_differentials_and_detect_the_flag_triple() {
    assert_eq!(
        stdout(example("rank-hochster")),
        concat!(
            "rvhochster[Flag, [], -1] = 1\n",
            "rvhochster[Flag, [1, 2], 0] = 1\n",
            "rvhochster[Flag, [1, 2, 3, 4, 5, 6], 1] = 1\n",
            "rvhochster[Flag, [1, 2, 3, 4, 5, 6], 2] = 0\n",
            "RVMKoszul(Flag, 4, [1, 1, 1, 1, 1, 1]){x1*x5 e2 e3 e4 e6}\n",
            "RVMTor(Flag, 4, [1, 1, 1, 1, 1, 1]){g0}\n"
        )
    );
}

#[test]
fn missing_vertices_void_complex_and_non_squarefree_fibers_are_exact() {
    let source = "
        rvideal Missing[[0,1]] {rvmon[[0,1]]{1,0}};
        rvideal Void[[0,1]] {rvmon[[0,1]]{0,0}};
        rvideal Edge[[0,1]] {rvmon[[0,1]]{1,1}};
        rvhochster[Missing,[0],-1];
        rvhochster[Missing,[0],0];
        rvhochster[Missing,[1],-1];
        rvhochster[Void,[],-1];
        rvhochster[Void,[0,1],-2];
        rvhochster[Edge,[0,1],0];
        rvmbetti[Edge,0,[2,1]];
        rvmbetti[Edge,1,[2,1]];
        rvmbetti[Edge,2,[2,1]];
        let raw=rvhface[Missing,[0],[]];
        assert rvTorClass(raw)==rvmtorgen[Missing,1,[1,0],0];
        ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvhochster[Missing, [0], -1] = 1\n",
            "rvhochster[Missing, [0], 0] = 0\n",
            "rvhochster[Missing, [1], -1] = 0\n",
            "rvhochster[Void, [], -1] = 0\n",
            "rvhochster[Void, [0, 1], -2] = 0\n",
            "rvhochster[Edge, [0, 1], 0] = 1\n",
            "rvmbetti[Edge, 0, [2, 1]] = 0\n",
            "rvmbetti[Edge, 1, [2, 1]] = 0\n",
            "rvmbetti[Edge, 2, [2, 1]] = 0\n"
        )
    );
}

#[test]
fn induced_vertex_sets_faces_and_ideal_hypotheses_are_checked() {
    let setup = "
        rvideal Square[[0,1]] {rvmon[[0,1]]{2,0}};
        rvideal Edge[[0,1]] {rvmon[[0,1]]{1,1}};
        let local=rvTorClass(rvhface[Edge,[0,1],[0]]);
    ";
    for (input, expected) in [
        ("rvhochster[Square,[0],0];", "squarefree"),
        ("show rvhface[Edge,[0,1],[0,1]];", "not a face"),
        ("show rvhface[Edge,[0],[1]];", "contained"),
        ("rvhochster[Edge,[1,0],0];", "strictly increasing"),
        ("rvhochster[Edge,[2],0];", "outside"),
        (
            "show rvTHigher(local,rvMForget(local));",
            "either all total or all multigraded",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{input}")));
        assert!(error.contains(expected), "{input}: {error}");
    }
}

#[test]
fn the_void_complex_has_no_large_wedge_slice_to_enumerate() {
    let width = 40.min(usize::BITS as usize);
    let indices: Vec<_> = (0..width).collect();
    let zero = vec!["0"; width].join(",");
    let multidegree = vec![1; width];
    let homological = width / 2;
    let source = format!(
        "rvideal Void[{indices:?}] {{rvmon[{indices:?}]{{{zero}}}}};
         rvmbetti[Void,{homological},{multidegree:?}];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!("rvmbetti[Void, {homological}, {multidegree:?}] = 0\n")
    );
}
