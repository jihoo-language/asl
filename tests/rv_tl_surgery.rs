mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn kirby_moves_and_an_actual_trefoil_surgery() {
    assert_eq!(
        stdout(example("rv-tl-surgery")),
        concat!(
            "RVCycloQ(16){[0, 1, 0, -1] / 1}\n",
            "RVCycloQ(16){[-1] / 1}\n",
            "RVCycloQ(16){[2, 0, -1, 0, 0, 0, 1] / 4}\n",
        )
    );
}

#[test]
fn doubled_normalization_matches_lens_and_orientation_controls() {
    let program = r"
        let unit = rvCycloRat[1](rvCyclo[20,0]);
        let trivial = rvCycloRat[1](rvCyclo[8,0]);
        assert rvTLSurgery[0,2,[1,1],[0,0]] == trivial;
        assert rvTLSurgeryTV[0,1,[],[0]] == trivial;
        let sphere = rvTLSurgeryTV[3,0,[],[]];
        assert sphere == rvCycloQGalois[11](rvLensTV[2,3,1,1]);
        assert rvTLSurgeryTV[3,1,[],[0]] == unit;
        assert rvTLSurgeryTV[3,1,[],[2]] == zero[RVCycloQ(20)];
        let l3 = rvTLSurgery[3,1,[],[3]];
        let mirror = rvTLSurgery[3,1,[],[-3]];
        assert mirror == rvCycloQGalois[-1](l3);
        let doubled = rvTLSurgeryTV[3,1,[],[3]];
        assert doubled == rvCycloQGalois[11](rvLensTV[2,3,3,1]);
        assert doubled == rvCycloQMul(rvCycloQMul(l3,mirror),sphere);
        assert rvTLSurgeryTV[3,2,[],[3,-1]] == doubled;
        assert rvTLSurgeryTV[3,1,[],[-3]] == doubled;
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn singular_forms_and_large_coloring_fibers_are_diagnosed_before_trace() {
    for (program, fragment) in [
        (
            "show rvTLSurgery[2,2,[1,1],[1,1]];",
            "requires a nonsingular linking matrix",
        ),
        (
            "show rvTLSurgery[2,1,[],[0]];",
            "requires a nonsingular linking matrix",
        ),
        (
            "show rvTLSurgeryTV[2,2,[1,1],[1]];",
            "requires one framing per component",
        ),
        (
            "show rvTLSurgeryTV[50,2,[1,1],[1,1]];",
            "limit of 4096 elementary crossings",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(fragment), "{program}\n{error}");
    }
    let braid: Vec<_> = (1..30).step_by(2).flat_map(|i| [i, i]).collect();
    let large = format!("show rvTLSurgeryTV[2,30,{braid:?},{:?}];", vec![1; 30]);
    let error = diagnostic(interpret(&large));
    assert!(
        error.contains("statement exceeds the computation limit")
            || error.contains("colorings exceed the computation limit"),
        "{error}"
    );
}
