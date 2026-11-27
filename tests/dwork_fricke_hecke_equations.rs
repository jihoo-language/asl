mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = "dwork-fricke-hecke-equations";

#[test]
fn exact_level_three_and_five_equations_have_the_audited_coefficients() {
    let out = stdout(example(EXAMPLE));
    assert!(out.contains("rvDworkHeckeEquation[3] = bidegree=(4, 4), terms=16"));
    assert!(out.contains("rvDworkHeckeEquation[5] = bidegree=(6, 6), terms=36"));
    assert!(out.contains("(0, 0): 16039090461081600"));
    assert!(out.contains("rvDworkHeckeCoeff[5, 0, 0] = 55365991016270105811049512960000"));
    assert!(out.contains("rvDworkHeckeCoeff[3, 3, 1] = 13116"));
    assert!(out.contains("rvDworkHeckeCoeff[3, 2, 2] = 3720006"));
    assert!(out.contains("rvDworkHeckeCoeff[3, 1, 1] = -540710244352"));
    assert_eq!(
        out.matches("certified symmetry, q coefficients through 8, and Kronecker congruence")
            .count(),
        2
    );
}

#[test]
fn specializes_with_rational_and_nonsplit_neighbors_and_preserves_multiplicities() {
    let out = stdout(example(EXAMPLE));
    let p41 = out
        .lines()
        .find(|line| line.starts_with("rvDworkHeckeNeighbors[3] = H=[0]"))
        .unwrap();
    assert!(p41.contains("(2, 1, [[26], [38], [1]])"));
    assert!(p41.ends_with("rational roots (value, multiplicity)=[([7], 1), ([26], 1)]"));
    let level_five = out
        .lines()
        .find(|line| line.starts_with("rvDworkHeckeNeighbors[5]"))
        .unwrap();
    assert!(level_five.contains("(4, 1, [[38], [30], [14], [34], [1]])"));
    assert!(level_five.ends_with("rational roots (value, multiplicity)=[([15], 1), ([23], 1)]"));
    let square = out
        .lines()
        .find(|line| line.starts_with("rvDworkHeckeNeighbors[3] = H=[4]"))
        .unwrap();
    assert!(square.ends_with("rational roots (value, multiplicity)=[([17], 2)]"));

    let p11 = stdout(interpret(
        "rvunram F[11,1]{0,1}; rvDworkHeckeNeighbors[3](rvUnramElem[F]{1});",
    ));
    assert!(p11.contains("(3, 1, [[3], [3], [2], [1]])"));
    assert!(p11
        .trim_end()
        .ends_with("rational roots (value, multiplicity)=[([3], 1)]"));
}

#[test]
fn level_characteristic_distinguishes_frobenius_from_inseparable_branch() {
    let out = stdout(example(EXAMPLE));
    let row = out
        .lines()
        .find(|line| line.starts_with("rvDworkHeckeNeighbors[3] = H=[0, 1, 0]"))
        .unwrap();
    assert!(row.ends_with(
        "rational roots (value, multiplicity)=[([1, 1, 0], 1), ([2, 1, 0], 3)], characteristic equals level, Frobenius correspondence"
    ));
    let prime_field = stdout(interpret(
        "rvunram F[5,1]{0,1}; rvDworkHeckeNeighbors[5](rvUnramElem[F]{2});",
    ));
    assert!(prime_field.contains("rational roots (value, multiplicity)=[([2], 6)]"));
}

#[test]
fn verifies_the_nontrivial_elliptic_isogeny_square_as_morphisms() {
    let out = stdout(example(EXAMPLE));
    assert!(out.contains("rvECMapDegree = 6"));
    assert!(out.contains("rvECMapCertify = degree=6, cleared identity=yes"));
}

#[test]
fn rejects_bad_levels_and_neighbor_coefficient_domains() {
    for (source, fragment) in [
        ("rvDworkHeckeEquation[2];", "support good levels 3 and 5"),
        ("rvDworkHeckeEquation[9];", "support good levels 3 and 5"),
        (
            "rvunram F[3,2]{0,1}; rvDworkHeckeNeighbors[3](rvUnramElem[F]{1});",
            "odd residue field at precision one",
        ),
        (
            "rvunram F[2,1]{0,1}; rvDworkHeckeNeighbors[3](rvUnramElem[F]{1});",
            "odd residue field at precision one",
        ),
        (
            "rvunram F[7,1]{0,1}; show rvDworkHauptmodulZ(rvUnramElem[F]{});",
            "requires a nonzero input",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(fragment), "{error}");
    }
}
