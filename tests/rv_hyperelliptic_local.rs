mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn component_and_graph_factors_and_direct_nodal_counts() {
    let output = stdout(example("rv-hyperelliptic-local"));
    for line in [
        "rvHyperLocal = graph=[1, -1], abelian=[1, 4, 14, 28, 49], polynomial=[1, 3, 10, 14, 21, -49], conductor=1, root_number=-1",
        "rvHyperLocal = graph=[1, 1], abelian=[1, -4, 14, -28, 49], polynomial=[1, -3, 10, -14, 21, 49], conductor=1, root_number=1",
        "rvHyperNodalPoints[1] = 25, fixed_vertices=3, fixed_edges=3, trace=-3",
        "rvHyperNodalPoints[1] = 19, fixed_vertices=3, fixed_edges=1, trace=3",
        "rvHyperFixed = factors=[2, 6], order=12",
        "rvHyperFixed = factors=[2, 10], order=20",
    ] {
        assert!(output.contains(line), "{line}\n{output}");
    }
    for (degree, count) in [25, 159, 1009, 7199, 50185].into_iter().enumerate() {
        assert!(output.contains(&format!("rvHyperNodalPoints[{}] = {count}, ", degree + 1)));
    }
    for (degree, count) in [19, 159, 1051, 7199, 50659].into_iter().enumerate() {
        assert!(output.contains(&format!("rvHyperNodalPoints[{}] = {count}, ", degree + 1)));
    }
    assert_eq!(
        output
            .matches(
                "graph=[1, 0, -1], abelian=[1], polynomial=[1, 0, -1], conductor=2, root_number=-1"
            )
            .count(),
        2,
        "{output}"
    );
}

#[test]
fn odd_degree_nonmonic_component_at_three() {
    let input = r"
        let a = rvHyperArithmetic[3,1]{0,9,18,1,2,4};
        let b = rvHyperArithmetic[3,2]{0,9,18,1,2,4};
        rvHyperLocal(a);
        rvHyperLocalFactors(a);
        rvHyperLocal(b);
        rvHyperLocalFactors(b);
        rvHyperNodalPoints[1](a);
        rvHyperNodalPoints[1](b);
    ";
    let output = stdout(interpret(input));
    for line in [
        "graph=[1, 1], abelian=[1, 0, 3], polynomial=[1, 1, 3, 3], conductor=1, root_number=1",
        "graph=[1, -1], abelian=[1, 0, 3], polynomial=[1, -1, 3, -3], conductor=1, root_number=-1",
        "Some(1), [0, 2, 0, 1], [0], [1, 0, 3]",
        "Some(1), [0, 1, 0, 2], [0], [1, 0, 3]",
    ] {
        assert!(output.contains(line), "{line}\n{output}");
    }
}

#[test]
fn smooth_even_degree_infinity_counts_distinguish_equal_graphs() {
    let input = r"
        let a = rvHyperArithmetic[7,1]{0,1,2,3};
        let b = rvHyperArithmetic[7,1]{0,1,2,4};
        let twist = rvHyperArithmetic[7,3]{0,1,2,4};
        assert rvHyperGeometric(a)==rvHyperGeometric(b);
        assert rvHyperGeometric(b)==rvHyperGeometric(twist);
        rvHyperLocal(a);
        rvHyperLocal(b);
        rvHyperLocal(twist);
        rvHyperNodalPoints[1](a);
        rvHyperNodalPoints[1](b);
        rvHyperNodalPoints[1](twist);
    ";
    let output = stdout(interpret(input));
    for line in [
        "polynomial=[1, 0, 7], conductor=0, root_number=1",
        "polynomial=[1, 4, 7], conductor=0, root_number=1",
        "polynomial=[1, -4, 7], conductor=0, root_number=1",
        "rvHyperNodalPoints[1] = 8, fixed_vertices=1, fixed_edges=0, trace=0",
        "rvHyperNodalPoints[1] = 12, fixed_vertices=1, fixed_edges=0, trace=-4",
        "rvHyperNodalPoints[1] = 4, fixed_vertices=1, fixed_edges=0, trace=4",
    ] {
        assert!(output.contains(line), "{line}\n{output}");
    }
}

#[test]
fn local_queries_enforce_the_arithmetic_domain_and_count_degree() {
    for (input, expected) in [
        (
            "let A=rvHyperArithmetic[7,1]{0,1,2,3}; rvHyperLocal(rvHyperGeometric(A));",
            "requires an RVHyperArithmetic",
        ),
        (
            "let A=rvHyperArithmetic[7,1]{0,1,2,3}; rvHyperNodalPoints[0](A);",
            "extension degree must be positive",
        ),
    ] {
        let output = diagnostic(interpret(input));
        assert!(output.contains(expected), "{input}\n{output}");
    }
}
