mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn disjoint_residue_charts_give_shells_balls_and_adlv_counts() {
    assert_eq!(
        stdout(example("rv-tree")),
        concat!(
            "rvTreeShell[3, 1, 3] = {[-3, 3]: 32, [-2, 2]: 24, [-1, 1]: 12, [0, 0]: 12}\n",
            "rvTreeShell[3, 1, 0] = {[0, 1]: 1}\n",
            "rvTreeShell[3, 1, 3] = {[-3, 4]: 64, [-2, 3]: 16}\n",
            "rvTreeBall[3, 1, 2] = {[-2, 3]: 16, [-1, 2]: 8, [0, 1]: 2}\n",
            "rvTreeCount[3, 1, 2, [-1, 2]] = 8\n",
            "rvTreeClosureCount[3, 1, 2, [-1, 2]] = 10\n",
            "rvTreeShell[3, 1, 3] = {[-3, 4]: 64, [-2, 3]: 16}\n",
            "rvTreeShell[3, 1, 3] = {[-5, 2]: 64, [-4, 1]: 16}\n",
            "rvTreeShell[3, 1, 2] = {[-2, 2]: 48, [-1, 1]: 18, [0, 0]: 6}\n",
        ),
    );
}

#[test]
fn tree_shells_may_fit_a_bound_that_the_corresponding_ball_exceeds() {
    let prefix = "qmodule X[0] { x0:0,x1:0 } {}; let a = rvU[0](rvId[X,[0],3,[1]]);";
    assert_eq!(
        stdout(interpret(&format!("{prefix} rvTreeShell[3,1,11](a);"))),
        "rvTreeShell[3, 1, 11] = {[0, 0]: 3072}\n",
    );
    for request in [
        "rvTreeBall[3,1,11](a);".into(),
        format!("rvTreeShell[3,1,{}](a);", usize::MAX),
    ] {
        let error = diagnostic(interpret(&format!("{prefix} {request}")));
        assert!(
            error.contains("exceeds the limit of 4096 vertices"),
            "{error}"
        );
    }
}

#[test]
fn zero_radius_and_closure_keep_the_true_determinant_valuation() {
    let program = r"
        qmodule X[0] { x0:0,x1:0 } {};
        let all = add(add(rvgen[X,X,[0],3,[1],0],rvgen[X,X,[0],3,[1],1]),
                      add(rvgen[X,X,[0],3,[1],2],rvgen[X,X,[0],3,[1],3]));
        let a = add(rvU[0](all),rvU[1](rvgen[X,X,[0],3,[1],3]));
        rvTreeBall[3,1,0](a);
        rvTreeCount[3,1,0,[0,1]](a);
        rvTreeClosureCount[3,1,0,[-2,3]](a);
        rvTreeClosureCount[3,1,0,[-2,2]](a);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvTreeBall[3, 1, 0] = {[0, 1]: 1}\n",
            "rvTreeCount[3, 1, 0, [0, 1]] = 1\n",
            "rvTreeClosureCount[3, 1, 0, [-2, 3]] = 1\n",
            "rvTreeClosureCount[3, 1, 0, [-2, 2]] = 0\n",
        ),
    );
}

#[test]
fn tree_queries_require_rank_two_a_rational_point_and_invertibility() {
    let prefix = "qmodule X[0] { x0:0,x1:0 } {};";
    for (request, expected) in [
        (
            "rvTreeShell[3,1,0](rvU[0](rvId[One,[0],3,[1]]));",
            "requires a rank-two Laurent Margolis endomorphism",
        ),
        (
            "rvTreeBall[3,1,0](rvU[0](rvgen[One,X,[0],3,[1],0]));",
            "requires a Laurent Margolis endomorphism",
        ),
        (
            "rvTreeShell[3,1,0](rvU[0](rvId[X,[0],7,[2]]));",
            "affine point in the specified subfield",
        ),
        (
            "rvTreeShell[3,1,0](rvU[0](rvgen[X,X,[0],7,[1],0]));",
            "requires an invertible rank-two Laurent endomorphism",
        ),
        (
            "rvTreeCount[3,1,0,[0]](rvU[0](rvId[X,[0],7,[1]]));",
            "one entry per Margolis dimension",
        ),
        (
            "rvTreeClosureCount[3,1,0,[1,-1]](rvU[0](rvId[X,[0],7,[1]]));",
            "coweight must be weakly increasing",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{prefix}\n{request}")));
        assert!(error.contains(expected), "{request}\n{error}");
    }
}
