mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn mixed_hopf_matches_independent_fusion_and_both_stabilizations() {
    let output = stdout(example("rv-tl-framed"));
    assert!(output.starts_with("RVCycloQ(28)"), "{output}");
    assert!(
        output.ends_with("rvTLLinking = [[8, 1], [1, 9]]\n"),
        "{output}"
    );
}

#[test]
fn nonuniform_blocks_transport_colors_through_inverse_and_braid_relation() {
    let program = r"
        let unlink = rvTLFramed[4,2,[],[1,2],[0,0]];
        assert rvTLFramed[4,2,[1,-1],[1,2],[0,0]] == unlink;
        assert rvTLFramed[4,2,[-1,1],[1,2],[0,0]] == unlink;
        let a = rvTLFramed[4,3,[1,2,1],[1,2],[1,-2]];
        let b = rvTLFramed[4,3,[2,1,2],[1,2],[1,-2]];
        assert a == b;
        let c = rvTLFramed[4,3,[1,2,1,-2,-1,-2],[1,2,0],[1,-2,4]];
        assert c == rvTLFramed[4,3,[],[1,2,0],[1,-2,4]];
        assert rvTLFramed[4,2,[1,1],[0,2],[9,-1]] == rvTLFramed[4,1,[],[2],[-1]];
        assert rvTLFramed[4,0,[],[],[]] == rvCycloRat[1](rvCyclo[24,0]);
        assert rvTLFramed[0,2,[1,1],[0,0],[-6,12]] == rvCycloRat[1](rvCyclo[8,0]);
        let d2 = rvCycloRat[1](add(rvCyclo[28,4],add(rvCyclo[28,0],rvCyclo[28,-4])));
        let reduced = rvTLColoredJones[5,2,2,[1,1,1]];
        assert rvTLFramed[5,2,[1,1,1],[2],[0]] == rvCycloQMul(d2,reduced);
        assert rvTLFramed[3,1,[],[1],[-170141183460469231731687303715884105728]]
            == rvTLFramed[3,1,[],[1],[-8]];
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn framing_matrices_keep_orientation_component_order_and_diagonal() {
    let program = r"
        rvTLLinking[2,[1,1],[2,-3]];
        rvTLLinking[2,[-1,-1],[-2,3]];
        rvTLLinking[0,[],[]];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvTLLinking = [[2, 1], [1, -3]]\n",
            "rvTLLinking = [[-2, -1], [-1, 3]]\n",
            "rvTLLinking = []\n",
        )
    );
}

#[test]
fn framed_links_check_components_colors_and_elementary_crossing_budget() {
    for (program, fragment) in [
        (
            "show rvTLFramed[3,2,[1],[1,2],[0]];",
            "requires one color per component",
        ),
        (
            "show rvTLFramed[3,2,[],[1,2],[0]];",
            "requires one framing per component",
        ),
        (
            "show rvTLFramed[3,2,[],[1,4],[0,0]];",
            "colors must be at most the level",
        ),
        ("rvTLLinking[2,[2],[0,0]];", "Artin index must satisfy"),
        (
            "show rvTLFramed[50,2,[1,1],[50,50],[0,0]];",
            "limit of 4096 elementary crossings",
        ),
        (
            "rvTLLinking[18446744073709551615,[],[]];",
            "calculation exceeds the storage limit",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(fragment), "{program}\n{error}");
    }
}
