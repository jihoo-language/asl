mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn twisted_cusp_has_extension_only_strata_and_exact_zeta_coefficients() {
    assert_eq!(
        stdout(example("rank-arithmetic")),
        concat!(
            "rvorbit[7, [1, 1, 2]] = degree 2, orbit [[1, 1, 2], [1, 1, 3]]\n",
            "rvdim[One, TwistedCusp, [0, 1, 2], 7, [1, 1, 2]] = 2\n",
            "rvclosed[One, TwistedCusp, [0, 1, 2], 1] = {0: [5], 4: [1], 8: [1]}\n",
            "rvclosed[One, TwistedCusp, [0, 1, 2], 6] = {0: [5, 4, 22, 57, 210, 659], ",
            "2: [0, 3, 0, 6, 0, 20], 4: [1, 0, 0, 0, 0, 0], 8: [1, 0, 0, 0, 0, 0]}\n",
            "rvzeta[One, TwistedCusp, [0, 1, 2], 6] = {0: [1, 5, 19, 77, 307, 1229, 4915], ",
            "2: [1, 0, 3, 0, 12, 0, 48], 4: [1, 1, 1, 1, 1, 1, 1], 8: [1, 1, 1, 1, 1, 1, 1]}\n"
        )
    );
}

#[test]
fn projective_normalization_and_alternate_field_presentations_preserve_orbits() {
    let source = "
        rvorbit[7, [1, 1, 2]];
        rvorbit[7, [2, 2, 3]];
        rvorbit[7, [1, 1, 3]];
        rvorbit[11, [0, 2, 4]];
        rvorbit[13, [0, 1, 3]];
        rvorbit[11, [5]];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvorbit[7, [1, 1, 2]] = degree 2, orbit [[1, 1, 2], [1, 1, 3]]\n",
            "rvorbit[7, [2, 2, 3]] = degree 2, orbit [[1, 1, 2], [1, 1, 3]]\n",
            "rvorbit[7, [1, 1, 3]] = degree 2, orbit [[1, 1, 2], [1, 1, 3]]\n",
            "rvorbit[11, [0, 2, 4]] = degree 3, orbit [[0, 1, 2], [0, 1, 4], [0, 1, 6]]\n",
            "rvorbit[13, [0, 1, 3]] = degree 3, orbit [[0, 1, 3], [0, 1, 5], [0, 1, 6]]\n",
            "rvorbit[11, [5]] = degree 1, orbit [[1]]\n"
        )
    );

    let program = format!(
        "{}\n{}",
        include_str!("../examples/rank-arithmetic.asl"),
        "rvspectrum[One, TwistedCusp, [0, 1, 2], 19];
         rvspectrum[One, TwistedCusp, [0, 1, 2], 25];"
    );
    let output = stdout(interpret(&program));
    assert!(
        output.ends_with(concat!(
            "rvspectrum[One, TwistedCusp, [0, 1, 2], 19] = {0: 241, 2: 30, 4: 1, 8: 1}\n",
            "rvspectrum[One, TwistedCusp, [0, 1, 2], 25] = {0: 241, 2: 30, 4: 1, 8: 1}\n"
        )),
        "{output}"
    );
}

#[test]
fn projective_points_lines_and_empty_selections_have_the_expected_counts() {
    let source = "
        rvclosed[One, One, [0], 5];
        rvzeta[One, One, [0], 5];
        rvclosed[One, One, [0, 1], 3];
        rvzeta[One, One, [0, 1], 3];
        rvclosed[One, One, [], 42];
        rvzeta[One, One, [], 42];
        rvclosed[One, One, [0, 1], 0];
        rvzeta[One, One, [0, 1], 0];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvclosed[One, One, [0], 5] = {1: [1, 0, 0, 0, 0]}\n",
            "rvzeta[One, One, [0], 5] = {1: [1, 1, 1, 1, 1, 1]}\n",
            "rvclosed[One, One, [0, 1], 3] = {1: [3, 1, 2]}\n",
            "rvzeta[One, One, [0, 1], 3] = {1: [1, 3, 7, 15]}\n",
            "rvclosed[One, One, [], 42] = {}\n",
            "rvzeta[One, One, [], 42] = {}\n",
            "rvclosed[One, One, [0, 1], 0] = {}\n",
            "rvzeta[One, One, [0, 1], 0] = {}\n"
        )
    );
}

#[test]
fn invalid_projective_points_fields_and_exhausted_budgets_are_diagnosed() {
    for (source, expected) in [
        ("rvorbit[7, [0, 0]];".into(), "zero vector"),
        ("rvorbit[7, []];".into(), "zero vector"),
        ("rvorbit[7, [1, 4]];".into(), "out of range"),
        ("rvorbit[5, [1, 1]];".into(), "reducible"),
        ("rvclosed[One, One, [1, 0], 2];".into(), "increasing"),
        (
            format!("rvzeta[One, One, [0], {}];", usize::BITS),
            "extension degree",
        ),
        ("rvspectrum[One, One, [0, 1, 2], 67];".into(), "4096"),
    ] {
        let error = diagnostic(interpret(&source));
        assert!(error.contains(expected), "{source}\n{error}");
    }
    let selected = (0..30)
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    let output = interpret(&format!("rvclosed[One, One, [{selected}], 1];"));
    assert!(output.stdout.is_empty(), "{output:?}");
    let error = diagnostic(output);
    assert!(error.contains("limit"), "{error}");
}
