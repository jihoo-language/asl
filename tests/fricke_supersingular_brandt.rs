mod common;

use common::{diagnostic, example, interpret, stdout};

const DEMO: &str = "fricke-supersingular-brandt";

#[test]
fn complete_fricke_divisor_retains_special_endpoints_and_hecke_multiplicities() {
    let out = stdout(example(DEMO));
    assert!(out.contains("geometric supersingular j values (p = 7) = [[6, 0]]"));
    assert!(out.contains("j polynomial (ascending) = [1, 1]"));
    assert!(out.contains("Fricke vertices = [[1, 0], [5, 0]]; unordered endpoints = [[[6, 0], [6, 0]], [[6, 0], [6, 0]]]"));
    assert!(out.contains("Fricke polynomial (ascending) = [5, 1, 1]"));
    assert!(out.contains("B3 = [[0, 4], [2, 2]]; weights = [1, 2]"));
    assert!(out.contains("B9 = [[4, 8], [4, 8]]; weights = [1, 2]"));
    assert!(out.contains("characteristic (descending) = [1, -2, -8]; cuspidal = [1, 2]"));
    assert_eq!(out.matches("B3^2 = B9 + 4I; B3 B5 = B5 B3").count(), 3);
}

#[test]
fn spectra_and_conjugate_vertices_match_independent_data() {
    let out = stdout(example(DEMO));
    assert!(out.contains("geometric supersingular j values (p = 23) = [[0, 0], [3, 0], [19, 0]]"));
    assert!(
        out.contains("characteristic (descending) = [1, -4, -5, 20, 0]; cuspidal = [1, 0, -5, 0]")
    );
    assert!(out.contains("common positive detailed-balance weights = [4, 1, 2, 4]"));
    assert!(out.contains("Fricke vertices = [[18, 0], [20, 7], [20, 34], [27, 0], [39, 0]]"));
    assert!(out.contains("Frobenius permutation = [0, 2, 1, 3, 4]"));
    assert!(out.contains("characteristic (descending) = [1, -2, -12, 10, 28, -16]"));
    assert!(out.contains("Frobenius +1 (rank 4) = [1, -4, -4, 18, -8]; -1 (rank 1) = [1, 2]"));
}

#[test]
fn good_levels_in_characteristic_five_and_typed_polynomial_coefficients() {
    let out = stdout(interpret(
        "rvunram F[5,1]{-2,0,1}; rvfrickess B[F]; rvFrickeSSData[B];
         rvFrickeSSRelations[B];
         assert rvFrickeSSJCoeff[B,0]==rvUnramElem[F]{0};
         assert rvFrickeSSJCoeff[B,1]==rvUnramElem[F]{1};
         assert rvFrickeSSFrickeCoeff[B,0]==neg(rvFrickeSSVertex[B,0]);
         assert rvFrickeSSFrickeCoeff[B,1]==rvUnramElem[F]{1};
         assert rvFrickeSSFrickeCoeff[B,8]==rvUnramElem[F]{0};",
    ));
    assert!(out.contains("B3^2 = B9 + 4I; good levels 3 and 9 in characteristic 5"));
    assert!(out.contains("geometric supersingular j values (p = 5) = [[0, 0]]"));
}

#[test]
fn typed_pushforward_handles_wide_cancellation_and_rejects_overflow() {
    let prefix = "rvunram F[7,1]{-3,0,1}; rvfrickess B[F];";
    let max = i128::MAX;
    let no_overflow = format!(
        "{prefix} assert rvFrickeSSHecke[9](rvFrickeSSDiv[B]{{{max},-{max}}})
         ==zero[RVFrickeSSDiv(B)];"
    );
    stdout(interpret(&no_overflow));
    let overflow = format!("{prefix} show rvFrickeSSHecke[9](rvFrickeSSDiv[B]{{{max},0}});");
    assert!(diagnostic(interpret(&overflow)).contains("exceeds the signed 128-bit range"));
}

#[test]
fn nominal_and_good_characteristic_guards() {
    for (source, fragment) in [
        (
            "rvunram F[7,1]{0,1}; rvfrickess B[F];",
            "requires a named quadratic residue field",
        ),
        (
            "rvunram F[7,2]{-3,0,1}; rvfrickess B[F];",
            "requires a named quadratic residue field",
        ),
        (
            "rvunram F[5,1]{-2,0,1}; rvfrickess B[F]; rvFrickeSSMatrix[B,5];",
            "requires level prime to the characteristic",
        ),
        (
            "rvunram F[7,1]{-3,0,1}; rvfrickess B[F]; rvfrickess C[F];
             show add(rvFrickeSSBasis[B,0],rvFrickeSSBasis[C,0]);",
            "RVFrickeSSDiv",
        ),
        (
            "rvunram F[7,1]{-3,0,1}; rvfrickess B[F]; show rvFrickeSSDiv[B]{0};",
            "supersingular Fricke divisor needs 2 coefficients",
        ),
        (
            "rvunram F[7,1]{-3,0,1}; rvfrickess B[F]; show rvFrickeSSVertex[B,2];",
            "supersingular Fricke vertex index is out of range",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(fragment), "{error}");
    }
}
