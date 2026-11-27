mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn actual_theta_pairing_and_mixed_genus_bridge() {
    assert_eq!(
        stdout(example("rv-hyperelliptic-graph")),
        concat!(
            "rvHyperGraph = vertices=[(0, 0, 2, 0), (1, 0, 1, 2), (1, 0, 1, 1), (2, 0, 2, 3), (3, 0, 2, 4)], ",
            "edges=[(1, 3, 0, (2, 1), 1, 1), (2, 3, 0, (2, 1), 1, 0), (0, 1, 1, (5, 1), 1, 3), (0, 2, 1, (5, 1), 1, 2), (1, 4, 2, (3, 1), 1, 5), (2, 4, 2, (3, 1), 1, 4)], ",
            "branches=[(0, 3, 2), (1, 3, 2), (2, 4, 2), (3, 4, 2), (4, 0, 2), (5, 0, 2)]\n",
            "rvHyperHarmonic = true\n",
            "rvHyperGenus = 2, toric=2, components=[0, 0, 0, 0, 0]\n",
            "rvHyperCycles = [[(0, -1), (1, 1), (2, -1), (3, 1)], [(2, -1), (3, 1), (4, -1), (5, 1)]]\n",
            "rvHyperMonodromy = [[14, 10], [10, 16]]\n",
            "rvHyperComponent = factors=[2, 62], exponent=62, order=124\n",
            "rvHyperComponent = factors=[4, 124], exponent=124, order=496\n",
            "rvHyperGraph = vertices=[(0, 1, 2, 0), (1, 1, 2, 1), (2, 0, 2, 2)], ",
            "edges=[(0, 1, 0, (1, 2), 2, 0), (0, 2, 1, (2, 1), 1, 2), (0, 2, 1, (2, 1), 1, 1)], ",
            "branches=[(0, 1, 2), (1, 1, 2), (2, 1, 2), (3, 2, 2), (4, 2, 2), (5, 0, 2), (6, 0, 2), (7, 0, 2)]\n",
            "rvHyperHarmonic = true\n",
            "rvHyperGenus = 3, toric=1, components=[1, 1, 0]\n",
            "rvHyperCycles = [[(1, -1), (2, 1)]]\n",
            "rvHyperMetric = [1, 4, 4]\n",
            "rvHyperMonodromy = [[8]]\n",
            "rvHyperComponent = factors=[8], exponent=8, order=8\n",
        )
    );
}

#[test]
fn smooth_and_purely_separating_fibers_have_trivial_component_groups() {
    let source = r"
        let smooth = rvHyperCover(rvStableLine[7]{[0,1],[1,1],[2,1],[3,1],[4,1],[5,1]});
        rvHyperGenus(smooth);
        rvHyperGraph(smooth);
        rvHyperCycles(smooth);
        rvHyperComponent(smooth);

        let bridge = rvHyperCover(rvStableLine[7]{[0,1],[7,1],[14,1],[1,1],[2,1],[3,1]});
        rvHyperGenus(bridge);
        rvHyperGraph(bridge);
        rvHyperMetric[2](bridge);
        rvHyperMonodromy[2](bridge);
        rvHyperComponent[2](bridge);

        // When the base branch length is already even the odd bridge is integral at multiplier one.
        let even = rvHyperCover(rvStableLine[7]{[0,1],[49,1],[98,1],[1,1],[2,1],[3,1]});
        rvHyperMetric(even);
        rvHyperComponent(even);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvHyperGenus = 2, toric=0, components=[2]\n",
            "rvHyperGraph = vertices=[(0, 2, 2, 0)], edges=[], branches=[(0, 0, 2), (1, 0, 2), (2, 0, 2), (3, 0, 2), (4, 0, 2), (5, 0, 2)]\n",
            "rvHyperCycles = []\n",
            "rvHyperComponent = factors=[], exponent=1, order=1\n",
            "rvHyperGenus = 2, toric=0, components=[1, 1]\n",
            "rvHyperGraph = vertices=[(0, 1, 2, 0), (1, 1, 2, 1)], edges=[(0, 1, 0, (1, 2), 2, 0)], branches=[(0, 1, 2), (1, 1, 2), (2, 1, 2), (3, 0, 2), (4, 0, 2), (5, 0, 2)]\n",
            "rvHyperMetric = [1]\n",
            "rvHyperMonodromy = []\n",
            "rvHyperComponent = factors=[], exponent=1, order=1\n",
            "rvHyperMetric = [1]\n",
            "rvHyperComponent = factors=[], exponent=1, order=1\n",
        )
    );
}

#[test]
fn even_split_at_one_quadratic_embedding_and_smooth_at_the_other() {
    let source = r"
        rvqpadic K[7,2,3];
        let zeroK = zero[RVQpQuad(K)];
        let one = rvQpQuad[K,1,0,1];
        let beta = rvQpQuad[K,3,-1,1];
        let cover = rvHyperCover(rvStableLineQuad[K]{[zeroK,one],[beta,one],[one,one],[rvQpQuad[K,2,0,1],one]});
        rvHyperGenus(cover);
        rvHyperCycles(cover);
        rvHyperMonodromy(cover);
        rvHyperComponent(cover);

        rvqpadic J[7,2,4];
        let zeroJ = zero[RVQpQuad(J)];
        let oneJ = rvQpQuad[J,1,0,1];
        let betaJ = rvQpQuad[J,3,-1,1];
        let other = rvHyperCover(rvStableLineQuad[J]{[zeroJ,oneJ],[betaJ,oneJ],[oneJ,oneJ],[rvQpQuad[J,2,0,1],oneJ]});
        rvHyperGenus(other);
        rvHyperComponent(other);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvHyperGenus = 1, toric=1, components=[0, 0]\n",
            "rvHyperCycles = [[(0, -1), (1, 1)]]\n",
            "rvHyperMonodromy = [[2]]\n",
            "rvHyperComponent = factors=[2], exponent=2, order=2\n",
            "rvHyperGenus = 1, toric=0, components=[1]\n",
            "rvHyperComponent = factors=[], exponent=1, order=1\n",
        )
    );
}

#[test]
fn relabeling_preserves_the_integral_component_group() {
    let source = r"
        let original = rvHyperCover(rvStableLine[7]{[0,1],[49,1],[1,1],[344,1],[2,1],[16809,1]});
        let reshuffled = rvHyperCover(rvStableLine[7]{[2,1],[344,1],[0,1],[16809,1],[49,1],[1,1]});
        rvHyperComponent(original);
        rvHyperComponent(reshuffled);
        let inverted = rvHyperCover(rvStableLine[7]{[1,0],[1,49],[1,1],[1,344],[1,2],[1,16809]});
        assert inverted == original;
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvHyperComponent = factors=[2, 62], exponent=62, order=124\n".repeat(2)
    );
}

#[test]
fn branch_inputs_and_metric_integrality_are_checked() {
    for (source, expected) in [
        ("show rvHyperCover(rvStableLine[2]{[0,1],[1,1],[2,1],[3,1]});", "requires odd residue characteristic"),
        ("show rvHyperCover(rvStableLine[7]{[0,1],[1,1],[2,1],[3,1],[4,1]});", "even number of at least four branch marks"),
        ("show rvHyperCover(rvStableTree(rvStableLine[7]{[0,1],[1,1],[2,1],[3,1]}));", "an RVStableMarks"),
        ("rvHyperGraph(rvStableLine[7]{[0,1],[1,1],[2,1],[3,1]});", "requires an RVHyperCover"),
        ("rvHyperMetric[0](rvHyperCover(rvStableLine[7]{[0,1],[1,1],[2,1],[3,1]}));", "ramification multiplier must be positive"),
        ("rvHyperComponent(rvHyperCover(rvStableLine[7]{[0,1],[7,1],[14,1],[1,1],[2,1],[3,1]}));", "edge 0 has half-integral length at multiplier 1"),
        ("rvHyperMonodromy[3](rvHyperCover(rvStableLine[7]{[0,1],[7,1],[14,1],[1,1],[2,1],[3,1]}));", "edge 0 has half-integral length at multiplier 3"),
        ("rvHyperComponent[18446744073709551615](rvHyperCover(rvStableLine[7]{[0,1],[49,1],[1,1],[344,1],[2,1],[16809,1]}));", "exceeds"),
    ] {
        let output = diagnostic(interpret(source));
        assert!(output.contains(expected), "{source}\n{output}");
    }
}
