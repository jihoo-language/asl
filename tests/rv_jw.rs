mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn colored_links_use_the_projected_cable_and_fusion_channels() {
    assert_eq!(
        stdout(example("rv-tl-colored")),
        concat!(
            "RVTL(5, 3, 1){[0, 1, 0, 1]: [-1, 0, 0, 0, 1, 0, 1, 0, -1, 0, -1] / 1, ",
            "[0, 1, 2, 1]: [1, 0, 0, 0, -1, 0, -1, 0, 1, 0, 1] / 1}\n",
            "RVCycloQ(28){[-2, 0, 2, 0, -1, 0, 1] / 1}\n",
            "RVCycloQ(28){[2, 0, 0, 0, -1, 0, 3, 0, -3, 0, 1] / 1}\n",
        )
    );
}

#[test]
fn block_projectors_are_idempotent_annihilate_and_truncate() {
    let program = r"
        let a = rvTLPath[4,[0,1,0,1,2,1]];
        let b = rvTLPath[4,[0,1,2,1,2,1]];
        let p = rvTLProject[2,3](add(a,b));
        assert p == add(rvTLProject[2,3](a),rvTLProject[2,3](b));
        assert rvTLProject[2,3](p) == p;
        assert rvTLE[2](p) == zero[RVTL(4,5,1)];
        assert rvTLE[3](p) == zero[RVTL(4,5,1)];
        assert rvTLProject[2,3](rvTLE[2](a)) == zero[RVTL(4,5,1)];
        assert rvTLProject[2,3](rvTLE[3](b)) == zero[RVTL(4,5,1)];
        assert rvTLProject[6,0](a) == a;
        assert rvTLProject[5,1](a) == a;
        let distant = rvTLPath[5,[0,1,0,1,0,1,0,1]];
        assert rvTLProject[2,2](rvTLE[5](distant)) == rvTLE[5](rvTLProject[2,2](distant));
        let v = rvTLPath[2,[0,1,2,1]];
        let delta = rvCycloRat[1](neg(add(rvCyclo[16,2],rvCyclo[16,-2])));
        assert rvTLProject[1,2](v) == v;
        let bent = rvTLProject[1,2](rvTLE[2](v));
        assert rvTLScale(delta,bent) == v;
        assert rvTLProject[1,3](v) == add(v,neg(rvTLScale(delta,bent)));
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn colored_unknots_mirrors_and_level_cutoff() {
    let program = r"
        let one = rvCycloRat[1](rvCyclo[28,0]);
        assert rvTLColoredJones[5,1,2,[]] == one;
        assert rvTLColoredJones[5,2,2,[1]] == one;
        assert rvTLColoredJones[5,2,2,[-1]] == one;
        assert rvTLColoredJones[5,2,2,[1,-1]] == rvTLColoredJones[5,2,2,[]];
        assert rvTLColoredJones[5,2,0,[1,1,-1]] == one;
        let t = rvTLColoredJones[5,2,2,[1,1,1]];
        assert rvTLColoredJones[5,2,2,[-1,-1,-1]] == rvCycloQGalois[-1](t);
        assert rvTLColoredJones[5,3,1,[1,-2,1,-2]] == rvTLJones[5,3,[1,-2,1,-2]];
        assert rvTLColoredJones[5,3,2,[1,-2,1,-2]] == rvTLColoredJones[5,3,2,[-2,1,-2,1]];
        assert rvTLColoredJones[2,2,2,[1,1,1]] == rvCycloRat[1](rvCyclo[16,0]);
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn invalid_blocks_colors_and_cabling_are_diagnosed() {
    for (program, expected) in [
        (
            "show rvTLProject[0,1](rvTLPath[2,[0,1]]);",
            "one-based strand index",
        ),
        (
            "show rvTLProject[2,1](rvTLPath[2,[0,1]]);",
            "outside the path strands",
        ),
        (
            "show rvTLProject[1,4](rvTLPath[2,[0,1,2,1,0]]);",
            "exceeds the level plus one",
        ),
        ("show rvTLProject[1,1](S.hi);", "requires a Temperley"),
        (
            "show rvTLColoredJones[5,2,6,[1]];",
            "color must be at most the level",
        ),
        ("show rvTLColoredJones[5,2,0,[0]];", "1 <= i < strands"),
        (
            "show rvTLColoredJones[65,2,65,[1]];",
            "4096 elementary crossings",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
