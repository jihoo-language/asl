mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn large_frame_changes_and_residue_chart_denominators_cancel() {
    assert_eq!(
        stdout(example("rv-projective-gauge")),
        concat!(
            "rvStableComponents = [[([0], infinity), ([1], 0), ([2], 1), ([3], 6)]]\n",
            "rvStableComponents = [[([0], infinity), ([1], 0), ([2], 1), ([3], 6)]]\n",
            "RVResSpace(7, 2){[1, 1]}\n",
            "RVResSpace(7, 2){[1, 2]}\n",
            "RVResSpace(7, 3){[1, 0, 0], [0, 1, 2]}\n",
        )
    );
}

#[test]
fn stable_brackets_can_be_oversized_as_field_elements_at_both_split_embeddings() {
    let program = r"
        let R = rvStableLine[7]{[1,0],[0,1],[10000000000000000000000000,1],[1,10000000000000000000000000]};
        rvStableComponents(R);
        rvqpadic K[7,2,3];
        let z = zero[RVQpQuad(K)];
        let one = rvQpQuad[K,1,0,1];
        let h = rvQpQuad[K,85070591730234615865843651857942052864,0,1];
        let j = rvQpQuad[K,-85070591730234615865843651857942052864,3,1];
        let Q = rvStableLineQuad[K]{[one,z],[z,one],[h,one],[j,one]};
        rvStableSplits(rvStableTree(Q));
        rvqpadic C[7,2,4];
        let cz = zero[RVQpQuad(C)];
        let co = rvQpQuad[C,1,0,1];
        let ch = rvQpQuad[C,85070591730234615865843651857942052864,0,1];
        let cj = rvQpQuad[C,-85070591730234615865843651857942052864,3,1];
        rvStableSplits(rvStableTree(rvStableLineQuad[C]{[co,cz],[cz,co],[ch,co],[cj,co]}));
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvStableComponents = [[([0], infinity), ([1], 0), ([2], 1), ([3], 4)]]\n",
            "rvStableSplits = [([0, 1], [2, 3], 1)]\n",
            "rvStableSplits = []\n",
        )
    );
}

#[test]
fn exact_gauges_preserve_collision_and_unrepresentable_scalar_diagnostics() {
    let collision = diagnostic(interpret(
        "show rvStableLine[7]{[1,0],[0,1],[10000000000000000000000000,10000000000000000000000000],[2,2]};",
    ));
    assert!(
        collision.contains("marks 2 and 3 collide exactly"),
        "{collision}"
    );
    let inverse = diagnostic(interpret(
        "rvqpadic K[7,2,3]; show rvQpQuadInv(rvQpQuadPow[46](rvQpQuad[K,3,-1,1]));",
    ));
    assert!(
        inverse.contains("quadratic exact integer arithmetic is out of range"),
        "{inverse}"
    );
}
