mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn rational_primal_and_dual_check_optimality_on_the_figure_eight() {
    let output = stdout(example("annular-rasmussen"));
    for line in [
        "leeS[Figure] = 0",
        "leeDtProfile[Figure] = [((0, 1), (1, 1), -1, -1), ((1, 1), (2, 1), 1, 1)]",
        "leeDt[Figure, 0] = -1",
        "leeDt[Figure, 1/2] = -1/2",
        "leeDt[Figure, 1] = 0",
        "leeDt[Figure, 3/2] = -1/2",
        "leeDt[Figure, 2] = -1",
        "leeDtWeight[1/2] = -3/2",
        "leeS[Trefoil] = 2",
        "leeS[Mirror] = -2",
        "leeDtProfile[Trefoil] = [((0, 1), (1, 1), 1, -2), ((1, 1), (2, 1), 5, 2)]",
        "leeDtProfile[Mirror] = [((0, 1), (2, 1), -3, 0)]",
        "leeS[A] = 2",
        concat!(
            "leeDtProfile[A] = [((0, 1), (1, 2), 1, -2), ((1, 2), (2, 3), 3, 2), ",
            "((2, 3), (1, 1), -1, -4), ((1, 1), (4, 3), 7, 4), ",
            "((4, 3), (3, 2), -1, -2), ((3, 2), (2, 1), 5, 2)]",
        ),
        "leeDt[A, 3/5] = 9/5",
        "leeDt[A, 2/3] = 5/3",
    ] {
        assert!(
            output.lines().any(|actual| actual == line),
            "{line}\n{output}"
        );
    }
    assert_eq!(
        output
            .lines()
            .filter(|line| *line == "leeDtWeight[1/2] = -1/2")
            .count(),
        2,
        "{output}"
    );
    assert_eq!(
        output
            .lines()
            .filter(|line| *line == "leeDtWeight[2/3] = 5/3")
            .count(),
        2,
        "{output}"
    );
}

#[test]
fn annular_conjugation_and_braid_stabilizations_are_computed() {
    let program = r"
        khbraid F[3,[-2,1,-2,1]];
        khbraid U[1,[]];
        khbraid P[2,[1]];
        khbraid N[2,[-1]];
        leeDtProfile[F];
        leeDt[U,0,1];
        leeDt[P,0,1];
        leeDt[N,0,1];
        leeDt[U,1,1];
        leeDt[P,1,1];
        leeDt[N,1,1];
        leeDt[P,2,4];
        let f = leeDtDual[F,2,4];
        assert add(f,neg(f)) == zero[LeeDual(F,0)];
        assert leeScale[2,2](f) == f;
        leeDtWeight[0,1](zero[LeeDual(F,0)]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "leeDtProfile[F] = [((0, 1), (1, 1), -1, -1), ((1, 1), (2, 1), 1, 1)]\n",
            "leeDt[U, 0] = -1\n",
            "leeDt[P, 0] = -1\n",
            "leeDt[N, 0] = -1\n",
            "leeDt[U, 1] = 0\n",
            "leeDt[P, 1] = 1\n",
            "leeDt[N, 1] = -1\n",
            "leeDt[P, 1/2] = 0\n",
            "leeDtWeight[0] = zero\n",
        )
    );
}

#[test]
fn annular_rasmussen_checks_parameters_and_nominal_dual_types() {
    for (source, message) in [
        ("khbraid U[1,[]]; leeDt[U,1,0];", "0 <= p/q <= 2"),
        ("khbraid U[1,[]]; leeDt[U,3,1];", "0 <= p/q <= 2"),
        (
            "khbraid U[1,[]]; show leeDtDualD(leeOrientation[U]);",
            "a Lee dual",
        ),
        (
            "khbraid U[1,[]]; show leeDtDualD(zero[LeeDual(U,-9223372036854775808)]);",
            "Lee dual homological degree is out of range",
        ),
        (
            "khbraid U[1,[]]; khbraid V[1,[]]; show leeDtPair(leeDtDual[U,0,1],leeOrientation[V]);",
            "requires the same braid and homological degree",
        ),
        (
            "khbraid U[1,[]]; show leeDtPair(leeDtDual[U,0,1],leeDtNull[U,0,1]);",
            "requires the same braid and homological degree",
        ),
        (
            "khbraid U[90,[]]; leeDtProfile[U];",
            "orientation cycle exceeds the computation limit",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(message), "{source}\n{error}");
    }
}
