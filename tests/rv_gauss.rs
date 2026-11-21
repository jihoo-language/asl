mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cusp_conic_and_fermat_compute_the_exact_projective_dual_ideals() {
    assert_eq!(
        stdout(example("rank-duality")),
        concat!(
            "rvgroebner[CuspGraph] = [x3, x2*x5 + x1*x4, x1^2*x5 + x0^2*x4, x1^3 + x0^2*x2]\n",
            "rvgroebner[CuspDual] = [x3]\n",
            "rvgroebner[CuspBidual] = [x8, x7]\n",
            "rvgroebner[ConicDual] = [x4]\n",
            "rvgroebner[ConicBidual] = [x8, x6]\n",
            "rvgroebner[FermatDual] = [x5^3 + x4^3 + x3^3]\n",
        )
    );
}

#[test]
fn group_saturation_drops_empty_and_repeated_components() {
    let program = r"
        let square = add(rvmon[[0,1,2]]{2,0,0},rvmon[[0,1,2]]{0,2,0});
        rvgauss EmptyGraph[[3,4,5]](square);
        rvdual EmptyDual[[3,4,5]](square);
        assert rvGenerator[EmptyGraph,0] == rvmon[[0,1,2,3,4,5]]{0,0,0,0,0,0};
        assert rvGenerator[EmptyDual,0] == rvmon[[3,4,5]]{0,0,0};
        rvdual EmptySource[[1]](rvmon[[0]]{3});
        assert rvGenerator[EmptySource,0] == rvmon[[1]]{0};

        // Only y=0 away from x=0 is smooth. Its tangent is constant.
        rvgauss MixedGraph[[3,4,5]](rvmon[[0,1,2]]{2,1,0});
        rvdual MixedDual[[3,4,5]](rvmon[[0,1,2]]{2,1,0});
        assert rvInIdeal[MixedGraph](rvmon[[0,1,2,3,4,5]]{0,1,0,0,0,0}) == S.hi;
        assert rvInIdeal[MixedGraph](rvmon[[0,1,2,3,4,5]]{1,0,0,0,0,0}) == S.lo;
        rvgroebner[MixedGraph];
        rvgroebner[MixedDual];

        // For two distinct lines the dual consists of two reduced points.
        rvdual TwoPoints[[3,4,5]](rvmon[[0,1,2]]{1,1,0});
        rvgroebner[TwoPoints];
        rvhilbert[TwoPoints,1000];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvgroebner[MixedGraph] = [x5, x3, x1]\n",
            "rvgroebner[MixedDual] = [x5, x3]\n",
            "rvgroebner[TwoPoints] = [x5, x3*x4]\n",
            "rvhilbert[TwoPoints, 1000] = 2\n",
        )
    );
}

#[test]
fn a_linear_graph_supports_interleaved_labels_witnesses_and_tor() {
    let program = r"
        let line = add(rvmon[[1,3,5]]{1,0,0},rvmon[[1,3,5]]{0,1,0});
        rvgauss Incidence[[0,2,4]](line);
        rvdual Point[[0,2,4]](line);
        let base = add(rvmon[[0,1,2,3,4,5]]{0,1,0,0,0,0},rvmon[[0,1,2,3,4,5]]{0,0,0,1,0,0});
        let normal = add(rvmon[[0,1,2,3,4,5]]{1,0,0,0,0,0},rvmon[[0,1,2,3,4,5]]{0,0,1,0,0,0});
        assert rvInIdeal[Incidence](base) == S.hi;
        assert rvInIdeal[Incidence](normal) == S.hi;
        assert rvGenerator[Incidence,2] == normal;
        assert rvWitness[Incidence,2](normal) == rvmon[[0,1,2,3,4,5]]{0,0,0,0,0,0};
        rvbetti[Incidence,1,1];
        rvgroebner[Point];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvbetti[Incidence, 1, 1] = 3\nrvgroebner[Point] = [x4, x2 + x0]\n"
    );
}

#[test]
fn invalid_hypersurfaces_coordinate_lists_and_generator_indices_report_errors() {
    let setup = "let a = rvmon[[0,1]]{1,1}; rvideal J[[0,1]] {};";
    for (statement, reason) in [
        (
            "rvgauss G[[2,3]](zero[RVPoly([0,1],2)]);",
            "nonzero homogeneous",
        ),
        ("rvdual D[[2,3]](rvmon[[0,1]]{0,0});", "positive degree"),
        ("rvdual D[[2,3]](S.lo);", "an RVPoly polynomial"),
        ("rvdual D[[2]](a);", "one dual coordinate"),
        ("rvgauss G[[1,2]](a);", "fresh dual coordinates"),
        ("rvgauss G[[3,2]](a);", "strictly increasing"),
        ("rvgauss J[[2,3]](a);", "already defined"),
        (
            "show rvGenerator[J,0];",
            "reduced basis index 0 is out of range",
        ),
        ("show rvGenerator[Missing,0];", "unknown rank ideal"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn a_large_sparse_gradient_respects_the_symbolic_storage_limit() {
    let base: Vec<_> = (0..16).collect();
    let dual: Vec<_> = (16..32).collect();
    let mut program = String::new();
    for pair in 0..8 {
        let degree = if pair < 6 { 3 } else { 1 };
        let mut terms = (0..=degree).map(|i| {
            let mut powers = [0; 16];
            powers[pair * 2] = i;
            powers[pair * 2 + 1] = degree - i;
            format!(
                "rvmon[{base:?}]{{{}}}",
                powers
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        });
        let first = terms.next().unwrap();
        let factor = terms.fold(first, |a, b| format!("add({a},{b})"));
        let expression = if pair == 0 {
            factor
        } else {
            format!("rvMul(f{}, {factor})", pair - 1)
        };
        writeln!(program, "let f{pair} = {expression};").unwrap();
    }
    writeln!(program, "rvgauss Oversized[{dual:?}](f7);").unwrap();
    let error = diagnostic(interpret(&program));
    assert!(
        error.contains("symbolic ideal calculation exceeds the storage limit"),
        "{error}"
    );
}
