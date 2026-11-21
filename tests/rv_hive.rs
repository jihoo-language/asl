mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn adjoint_multiplicity_stretches_and_translates() {
    assert_eq!(
        stdout(example("rv-hive")),
        concat!(
            "rvHiveCount = 2\n",
            "RVHive([0, 1, 2], [0, 1, 2], [1, 2, 3]){4}\n",
            "rvHiveAt[1, 1] = 5\n",
            "card = 2 / 2\n",
            "rvHiveCount = 8\n",
            "rvHiveAt[1, 1] = 0\n",
            "rvHiveCount = 2\n",
            "rvHorn = admissible\n",
            "rvHorn = [1], [2], [2] gives 0 < 1\n",
            "rvHorn = total mismatch (2 != 3)\n",
            "card = 0 / 0\n",
        ),
    );
}

#[test]
fn hives_are_discrete_and_low_rank_types_can_be_empty() {
    let program = r"
        let a = rvHiveWitness[[0,1,2],[0,1,2],[1,2,3]];
        let b = rvHive[[0,1,2],[0,1,2],[1,2,3]]{5};
        let id = \x : RVHive([0,1,2],[0,1,2],[1,2,3]) => x;
        assert id(a) == a;
        assert id(b) == b;
        show \x : RVHive([0],[0],[1]) => x;
        show \x : RVHive([0,0],[0,2],[1,1]) => x;
        rvHiveCount[[],[],[]];
        show rvHive[[],[],[]]{};
        rvHiveCount[[-7],[2],[-5]];
        rvHiveCount[[0],[0],[1]];
        show rvHive[[0,0],[0,2],[0,2]]{};
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "{}\n{}\n",
            "rvHiveCount = 1\n",
            "RVHive([], [], []){}\n",
            "rvHiveCount = 1\nrvHiveCount = 0\n",
            "RVHive([0, 0], [0, 2], [0, 2]){}\n",
        ),
    );
}

#[test]
fn horn_uses_positive_higher_subset_triples() {
    let program = r"
        rvHorn[[0,0,0,1],[0,0,1,3],[0,1,2,2]];
        rvHiveCount[[0,0,0,1],[0,0,1,3],[0,1,2,2]];
        rvHorn[[-2,-2,-2,-1],[3,3,4,6],[1,2,3,3]];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvHorn = [1, 2], [2, 3], [2, 3] gives 2 < 3\n",
            "rvHiveCount = 0\n",
            "rvHorn = [1, 2], [2, 3], [2, 3] gives 4 < 5\n",
        ),
    );
}

#[test]
fn reports_and_types_respect_validation_and_computation_limits() {
    assert_eq!(
        stdout(interpret(
            "rvHiveCount[[0,4100,8200],[0,4100,8200],[4100,8200,12300]];",
        )),
        "rvHiveCount = 4101\n",
    );
    for (program, expected) in [
        ("rvHiveCount[[0],[0,0],[0]];", "must have the same rank"),
        (
            "rvHiveCount[[1,0],[0,0],[0,1]];",
            "must be weakly increasing",
        ),
        (
            "show rvHive[[0],[0],[0]]{1};",
            "needs 0 interior coordinates",
        ),
        (
            "show rvHive[[0,1,2],[0,1,2],[1,2,3]]{3};",
            "rhombus inequality fails",
        ),
        (
            "show rvHiveWitness[[0,0],[0,2],[1,1]];",
            "has no integral witness",
        ),
        ("show rvHive[[0],[0],[1]]{};", "totals must satisfy"),
        (
            "rvHiveAt[1,2](rvHiveWitness[[0],[0],[0]]);",
            "outside the triangle",
        ),
        (
            r"card \x : RVHive([0,4100,8200],[0,4100,8200],[4100,8200,12300]) => S.hi;",
            "enumeration limit of 4096 elements",
        ),
        (
            "rvHiveCount[[170141183460469231731687303715884105727],[1],[1]];",
            "integer arithmetic is out of range",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
