mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn split_place_cancellation_changes_the_initial_matroid_at_real_valuation_depth() {
    let output = stdout(example("rv-tropical-flag"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines[0], "rvQpQuadResidue = 0");
    let next: u128 = lines[1].rsplit_once('=').unwrap().1.trim().parse().unwrap();
    assert!(next != 0 && next.is_multiple_of(5u128.pow(19)));
    assert_eq!(
        &lines[2..],
        [
            "rvValPlucker = 19",
            "rvValInitial = (0, [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3]])",
            "rvValInitial = (-19, [[0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])",
            "rvValExchange = (2, 19, 0)",
            "rvValTie = (0, [2, 3])",
            "rvValInitial = (19, [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3]])",
            "rvValSupport = [[]]",
            "rvValSupport = [[0, 1, 2, 3]]",
        ],
    );
}

#[test]
fn exact_rank_three_agrees_with_a_rational_filtration_and_row_changes() {
    let source = r"
        let a = rvValMat[3,5]{[1,0,0,1,1],[0,1,0,1,10],[0,0,1,1,28]};
        let b = rvValMat[3,5]{[3,3,0,6,33],[0,1,0,1,10],[0,0,1,1,28]};
        assert a == b;
        rvphodge D[3] {
            {(1,1),(2,1),(3,1),(6,1),(9,1)}, {},
            {(2,[1,0,0,1,1]), (1,[0,1,0,1,10]), (1,[0,0,1,1,28]),
             (0,[0,0,0,1,0]), (0,[0,0,0,0,1])}
        };
        assert rvValFil[D,1] == a;
        rvValPlucker[[4,3,2]](a);
        rvValPlucker[[1,3,4]](a);
        rvValPlucker[[0,3,4]](a);
        let d = rvValDual(a);
        rvValPlucker[[0,1]](d);
        rvValTie[[0,3],[1,2,3,4]](a);
        assert rvValDual(d) == a;

        let select = \x : S => (a, d);
        assert fst(select(S.hi)) == a;
        assert snd(select(S.lo)) == d;
    ";
    let output = stdout(interpret(source));
    assert_eq!(output,
        "rvValPlucker = 2\nrvValPlucker = 3\nrvValPlucker = 2\nrvValPlucker = 2\nrvValTie = (2, [1, 4])\n");
}

#[test]
fn loops_parallel_elements_and_all_infinite_ties() {
    let output = stdout(interpret(
        r"
        let M = rvValMat[3,4]{[1,0,1,0],[0,1,0,0]};
        rvValSupport(M);
        rvValPlucker[[1,3]](M);
        rvValTie[[3],[0,1,2]](M);
        rvValSupport(rvValDual(M));
        rvValExchange[[0,1],[1,2],0](M);
    ",
    ));
    assert_eq!(output,
        "rvValSupport = [[0, 1], [1, 2]]\nrvValPlucker = infinity\nrvValTie = (infinity, [0, 1, 2])\nrvValSupport = [[0, 3], [2, 3]]\nrvValExchange = (2, 0, 0)\n");
}

#[test]
fn boundary_ranks_are_singletons_and_weights_use_widened_costs() {
    let output = stdout(interpret(
        r"
        let M = rvValMat[2,2]{};
        let D = rvValDual(M);
        assert D == rvValMat[2,2]{[2,0],[0,2]};
        assert rvValDual(rvValMat[2,0]{}) == rvValMat[2,0]{};
        rvValInitial[[-9223372036854775808,-9223372036854775808]](D);
        rvValPlucker[[]](M);
        card (\v : RVValMat(2,2,0) => S.hi);
        card (\v : RVValMat(2,2,2) => S.hi);
    ",
    ));
    assert_eq!(output,
        "rvValInitial = (18446744073709551616, [[0, 1]])\nrvValPlucker = 0\ncard = 1 / 1\ncard = 1 / 1\n");
}

#[test]
fn invalid_shapes_and_tropical_queries_fail_with_distinct_diagnostics() {
    let setup = "let M = rvValMat[3,4]{[1,0,1,0],[0,1,0,0]};";
    for (code, reason) in [
        ("show rvValMat[9,2]{[1,1]};", "not prime"),
        ("show rvValMat[3,65]{};", "ground-set limit of 64"),
        ("show rvValMat[3,1]{[1],[1]};", "rank exceeds"),
        ("show rvValMat[3,2]{[1]};", "rows of length 2"),
        ("show rvValMat[3,2]{[1,1],[2,2]};", "no finite basis"),
        ("rvValPlucker[[0,0]](M);", "repeated element"),
        ("rvValPlucker[[0,4]](M);", "outside its ground set"),
        ("rvValPlucker[[1]](M);", "subset of size 2"),
        ("rvValInitial[[0]](M);", "needs 4 weights"),
        ("rvValExchange[[0,2],[1,2],0](M);", "finite left basis"),
        ("rvValExchange[[0,1],[1,2],1](M);", "outside the right"),
        ("rvValTie[[0],[1]](M);", "subset of size 3"),
        ("show rvValDual(S.lo);", "requires an RVValMat"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{code}")));
        assert!(error.contains(reason), "{code}\n{error}");
    }
    let row = format!("[{}]", vec!["0"; 64].join(","));
    let source = format!("show rvValMat[3,64]{{{}}};", vec![row; 10].join(","));
    let error = diagnostic(interpret(&source));
    assert!(error.contains("Plücker-minor work limit"), "{error}");
}
