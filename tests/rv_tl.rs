mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn level_three_braid_has_exact_trace_and_galois_conjugate() {
    assert_eq!(
        stdout(example("rv-tl-jones")),
        concat!(
            "RVCycloQ(20){[1, 0, 0, 0, -1, 0, 1] / 1}\n",
            "RVCycloQ(20){[2, 0, 0, 0, 1, 0, -1] / 1}\n",
            "RVCycloQ(20){[-1, 0, 1, 0, -2] / 1}\n",
            "RVCycloQ(20){[0, 0, 0, 0, -2, 0, 2] / 1}\n",
        )
    );
}

#[test]
fn relations_hold_on_the_actual_bounded_path_vectors() {
    let program = r"
        let a = rvTLPath[3,[0,1,0,1]];
        let b = rvTLPath[3,[0,1,2,1]];
        let c = rvTLPath[3,[0,1,2,3,2]];
        let d = rvCycloRat[1](neg(add(rvCyclo[20,2],rvCyclo[20,-2])));
        assert rvTLE[2](rvTLE[2](a)) == rvTLScale(d,rvTLE[2](a));
        assert rvTLE[2](rvTLE[2](b)) == rvTLScale(d,rvTLE[2](b));
        assert rvTLE[1](rvTLE[2](rvTLE[1](a))) == rvTLE[1](a);
        assert rvTLE[2](rvTLE[1](rvTLE[2](b))) == rvTLE[2](b);
        assert rvTLE[1](b) == zero[RVTL(3,3,1)];
        assert rvTLE[1](rvTLE[3](c)) == rvTLE[3](rvTLE[1](c));
        assert rvTLE[2](rvTLE[2](c)) == rvTLScale(d,rvTLE[2](c));
        assert rvTLBraid[[1,2,1]](a) == rvTLBraid[[2,1,2]](a);
        assert rvTLBraid[[1,2,1]](b) == rvTLBraid[[2,1,2]](b);
        assert rvTLBraid[[2,-2,1,-1]](a) == a;
        assert rvTLBraid[[-1,1,-2,2]](b) == b;
        assert rvTLBraid[[1,-2]](a) == rvTLBraid[[1]](rvTLBraid[[-2]](a));
        assert rvTLBraid[[1,-2]](add(a,b)) == add(rvTLBraid[[1,-2]](a),rvTLBraid[[1,-2]](b));
        assert add(neg(a),a) == zero[RVTL(3,3,1)];
        assert rvCycloQMul(rvTLCoord[[0,1,0,1]](rvTLE[2](b)),d) == rvCycloRat[1](rvCyclo[20,0]);
        let top = rvTLPath[2,[0,1,2,1,2]];
        let delta = rvCycloRat[1](neg(add(rvCyclo[16,2],rvCyclo[16,-2])));
        assert rvTLE[3](top) == rvTLScale(delta,top);
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn closure_normalization_and_mirror_work_at_multiple_levels() {
    for level in [1, 2, 3, 5] {
        let order = 4 * (level + 2);
        let program = format!(
            "
            let one = rvCycloRat[1](rvCyclo[{order},0]);
            let delta = rvCycloRat[1](neg(add(rvCyclo[{order},2],rvCyclo[{order},-2])));
            assert rvTLTrace[{level},0,0,[]] == one;
            assert rvTLTrace[{level},0,1,[]] == zero[RVCycloQ({order})];
            assert rvCycloQMul(rvTLJones[{level},0,[]],delta) == one;
            assert rvTLJones[{level},1,[]] == one;
            assert rvTLJones[{level},2,[]] == delta;
            assert rvTLJones[{level},2,[1]] == one;
            assert rvTLJones[{level},2,[-1]] == one;
            let trefoil = rvCycloRat[1](add(add(rvCyclo[{order},-4],rvCyclo[{order},-12]),neg(rvCyclo[{order},-16])));
            assert rvTLJones[{level},2,[1,1,1]] == trefoil;
            assert rvTLJones[{level},2,[-1,-1,-1]] == rvCycloQGalois[-1](trefoil);
            assert rvTLJones[{level},3,[1,1,1,-2]] == trefoil;
            let figure = rvCycloRat[1](add(add(rvCyclo[{order},-8],rvCyclo[{order},8]),add(rvCyclo[{order},0],neg(add(rvCyclo[{order},-4],rvCyclo[{order},4])))));
            assert rvTLJones[{level},3,[1,-2,1,-2]] == figure;
            assert rvTLJones[{level},4,[1,-2,1,-2,3]] == figure;
            "
        );
        assert_eq!(stdout(interpret(&program)), "");
    }
}

#[test]
fn rational_multiplication_is_reduced_and_galois_is_multiplicative() {
    let program = r"
        let a = rvCycloRat[10](rvCycloScale[6](add(rvCyclo[20,3],rvCyclo[20,-4])));
        let b = rvCycloRat[3](rvCycloScale[5](rvCyclo[20,1]));
        let c = rvCycloQMul(a,b);
        assert c == rvCycloRat[1](add(rvCyclo[20,4],rvCyclo[20,-3]));
        assert rvCycloQGalois[7](c) == rvCycloQMul(rvCycloQGalois[7](a),rvCycloQGalois[7](b));
        assert rvCycloQGalois[3](rvCycloQGalois[7](c)) == c;
        assert rvCycloQMul(zero[RVCycloQ(20)],b) == zero[RVCycloQ(20)];
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn malformed_paths_words_fields_and_sizes_are_diagnosed() {
    let oversized = std::iter::repeat_n("1", 4097).collect::<Vec<_>>().join(",");
    let oversized = format!("show rvTLTrace[1,2,0,[{oversized}]];");
    for (program, expected) in [
        ("show zero[RVTL(0,0,0)];", "level must be at least 1"),
        ("show zero[RVTL(1,2,2)];", "height exceeds the level"),
        ("show rvTLPath[3,[]];", "must start at zero"),
        ("show rvTLPath[3,[1,2]];", "starting at zero"),
        ("show rvTLPath[3,[0,2]];", "unit steps"),
        ("show rvTLPath[1,[0,1,2,1]];", "level bounds"),
        (
            "show rvTLCoord[[0,1,0]](rvTLPath[3,[0,1]]);",
            "indicated length",
        ),
        ("show rvTLTrace[3,2,0,[0]];", "1 <= i < strands"),
        ("show rvTLJones[3,2,[-2]];", "1 <= i < strands"),
        ("show rvTLE[1](rvTLPath[3,[0]]);", "1 <= i < strands"),
        (
            "show rvTLScale(rvCycloRat[1](rvCyclo[16,0]),rvTLPath[3,[0]]);",
            "RVCycloQ(20)",
        ),
        (
            "show rvCycloQGalois[5](rvCycloRat[1](rvCyclo[20,0]));",
            "must be coprime",
        ),
        (r"show \x: RVTL(3,3,1) => x;", "cannot enumerate"),
        (oversized.as_str(), "4096 letters"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
    assert_eq!(
        stdout(interpret(
            r"card \x: RVTL(3,3,2) => S.hi; show zero[RVTL(3,3,2)];"
        )),
        "card = 1 / 1\nRVTL(3, 3, 2){}\n"
    );
}
