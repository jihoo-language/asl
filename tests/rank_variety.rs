mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn nonbinary_points_detect_the_cuspidal_rank_variety() {
    assert_eq!(
        stdout(example("rank-variety")),
        concat!(
            "rvdim[One, Cusp, [0, 1, 2], 7, [1, 3, 1]] = 2\n",
            "rvdim[One, Cusp, [0, 1, 2], 7, [2, 3, 1]] = 0\n",
            "rvspectrum[One, Cusp, [0, 1, 2], 7] = {0: 16, 2: 4, 4: 1}\n",
            "rvdim[One, Cusp, [0, 1, 2], 11, [3, 4, 1]] = 2\n",
            "rvdim[One, Cusp, [0, 1, 2], 11, [2, 4, 1]] = 0\n",
            "rvspectrum[One, Cusp, [0, 1, 2], 11] = {0: 64, 2: 8, 4: 1}\n",
            "rvdim[One, Cusp, [0, 1, 2], 7, [0, 0, 1]] = 4\n",
            "rvdim[One, Cusp, [0, 1, 2], 7, [1, 0, 0]] = 2\n",
        )
    );
}

#[test]
fn homology_dimensions_at_both_endpoints_multiply() {
    let declarations = include_str!("../examples/rank-variety.asl")
        .split("\n// With lambda_2")
        .next()
        .unwrap();
    let projective_plane = include_str!("../examples/projective-plane.asl")
        .split("\nbetti RP2;")
        .next()
        .unwrap();
    let program = format!(
        "{declarations}\n{projective_plane}\n{}",
        r"
        rvdim[Cusp, Cusp, [0,1,2], 7, [0,0,0]];
        rvspectrum[Cusp, Cusp, [0,1,2], 7];
        rvspectrum[One, RP2, [0,1], 7];
        rvdim[Cusp, RP2, [0,1], 7, [0,1]];
        rvdim[Cusp, RP2, [], 11, []];
        rvspectrum[Cusp, RP2, [], 11];
    "
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvdim[Cusp, Cusp, [0, 1, 2], 7, [0, 0, 0]] = 36\n",
            "rvspectrum[Cusp, Cusp, [0, 1, 2], 7] = {0: 16, 4: 4, 16: 1}\n",
            "rvspectrum[One, RP2, [0, 1], 7] = {1: 4, 3: 1}\n",
            "rvdim[Cusp, RP2, [0, 1], 7, [0, 1]] = 0\n",
            "rvdim[Cusp, RP2, [], 11, []] = 18\n",
            "rvspectrum[Cusp, RP2, [], 11] = {}\n",
        )
    );
}

#[test]
fn sparse_large_coefficient_degrees_and_projective_point_counts() {
    let max = usize::MAX;
    let program = format!(
        "qmodule T[0] {{ a:{max}, b:{} }} {{ 0(a)=>[b] }};\n\
         rvdim[T,One,[0],283,[255]];\n\
         rvspectrum[One,One,[0],283];\n\
         rvspectrum[One,One,[0,1,2],11];",
        max - 1,
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvdim[T, One, [0], 283, [255]] = 0\n",
            "rvspectrum[One, One, [0], 283] = {1: 1}\n",
            "rvspectrum[One, One, [0, 1, 2], 11] = {1: 73}\n",
        )
    );
    if usize::BITS >= 64 {
        assert_eq!(
            stdout(interpret(
                "qmodule T[0] { a:1,b:0 } { 0(a)=>[b] };\n\
                 rvdim[T,One,[0],9223372036854775811,[9223372036854775807]];\n\
                 rvspectrum[One,One,[0],9223372036854775811];"
            )),
            concat!(
                "rvdim[T, One, [0], 9223372036854775811, [9223372036854775807]] = 0\n",
                "rvspectrum[One, One, [0], 9223372036854775811] = {1: 1}\n"
            )
        );
    }
}

#[test]
fn rank_varieties_validate_fields_points_and_primitive_selections() {
    for (program, expected) in [
        (
            "rvdim[One, One, [0], 1, [1]];",
            "positive polynomial degree",
        ),
        ("rvdim[One, One, [0], 5, [1]];", "reducible over F2"),
        ("rvspectrum[One, One, [0], 9];", "reducible over F2"),
        ("rvdim[One, One, [0,1], 7, [2]];", "one field coefficient"),
        (
            "rvdim[One, One, [0,1], 7, [4,2]];",
            "field coefficient 4 is out of range",
        ),
        ("rvspectrum[One, One, [0,0], 3];", "strictly increasing"),
        (
            "rvspectrum[One, One, [0,1,2,3,4,5], 11];",
            "projective spectrum exceeds the limit",
        ),
        (
            "qmodule U[0] {a:0} {}; rvdim[U, One, [1], 7, [0]];",
            "not declared for this module",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
