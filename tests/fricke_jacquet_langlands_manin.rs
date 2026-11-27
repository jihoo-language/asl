mod common;

use common::{diagnostic, example, interpret, stdout};

const DEMO: &str = "fricke-jacquet-langlands-manin";

#[test]
fn exact_relative_geometry_and_cuspidal_quotients_match_the_two_levels() {
    let out = stdout(example(DEMO));
    assert!(out.contains("level 46: 72 Manin symbols; relative dimension 13; cuspidal dimension 10; cuspidal J+ dimension 5; cuspidal J+ W2+ dimension 3"));
    assert!(out.contains("level 82: 126 Manin symbols; relative dimension 21; cuspidal dimension 18; cuspidal J+ dimension 9; cuspidal J+ W2+ dimension 4"));
    assert!(out.contains("cusp denominators [1, 2, 23, 46]; boundary = [1, 0, 0, -1]"));
    assert!(out.contains("cusp denominators [1, 2, 23, 46]; boundary = [0, -1, 1, 0]"));
    assert!(out.contains("T3 sector 3: characteristic (descending) = [1, 0, -5, 0]"));
}

#[test]
fn comparison_checks_both_exact_directions_and_the_frobenius_sign() {
    let out = stdout(example(DEMO));
    assert!(out.contains("Jacquet–Langlands dimension 3; determinant ="));
    assert!(out.contains("Jacquet–Langlands dimension 4; determinant ="));
    assert!(out.contains("T3: characteristic (descending) = [1, 0, -5, 0]"));
    assert!(out.contains("T3: characteristic (descending) = [1, 2, -4, -6, 4]"));
    assert!(out.contains("-Wp / Frobenius: characteristic (descending) = [1, -3, 3, -1]"));
    assert!(out.contains("-Wp / Frobenius: characteristic (descending) = [1, -2, 0, 2, -1]"));
    assert_eq!(out.matches("Manin to Brandt =").count(), 2);
    assert_eq!(out.matches("Brandt to Manin =").count(), 2);
}

#[test]
fn composite_good_hecke_and_atkin_operators_act_on_open_paths() {
    stdout(interpret(
        "rvmanin M[23]; let s=rvManinPath[M]{[-43,19],[312,5]};
         assert rvManinHecke[1](s)==s;
         assert rvManinHecke[3](rvManinHecke[5](s))==rvManinHecke[15](s);
         assert rvManinHecke[3](rvManinHecke[3](s))==add(rvManinHecke[9](s),rvManinScale[3](s));
         assert rvManinJ(rvManinJ(s))==s;
         assert rvManinAtkin[2](rvManinAtkin[2](s))==s;
         assert rvManinAtkin[23](rvManinAtkin[2](s))==rvManinAtkin[2](rvManinAtkin[23](s));
         assert rvManinAtkin[2](rvManinHecke[3](s))==rvManinHecke[3](rvManinAtkin[2](s));
         assert rvManinPath[M]{[312,5],[-43,19]}==neg(s);",
    ));
}

#[test]
fn comparison_guards_the_joint_sector_and_nominal_domains() {
    let prefix = "rvmanin M[23]; rvunram F[23,1]{-5,0,1}; rvfrickess D[F]; rvfrickejl X[M,D];";
    for (test, expected) in [
        (
            "show rvFrickeJLToBrandt[X](rvManinPath[M]{[1,0],[0,1]});",
            "requires a cuspidal J+ W2+ Manin class",
        ),
        (
            "show rvFrickeJLToManin[X](rvFrickeSSBasis[D,0]);",
            "requires a Brandt divisor of degree zero",
        ),
        (
            "rvmanin Other[23]; show rvFrickeJLToBrandt[X](rvManinBasis[Other,3,0]);",
            "Manin-to-Brandt transfer",
        ),
        (
            "rvmanin Other[41]; rvfrickejl Bad[Other,D];",
            "level must be twice the Brandt characteristic",
        ),
    ] {
        let result = diagnostic(interpret(&format!("{prefix} {test}")));
        assert!(result.contains(expected), "{result}");
    }
}

#[test]
fn bad_level_bad_cusp_and_wrong_operator_diagnostics() {
    for (source, expected) in [
        ("rvmanin M[21];", "21 is not prime"),
        ("rvmanin M[47];", "require prime 5 ≤ p ≤ 43"),
        (
            "rvmanin M[23]; show rvManinSymbol[M,2,4];",
            "must be primitive modulo 2p",
        ),
        (
            "rvmanin M[23]; show rvManinPath[M]{[0,0],[1,2]};",
            "cannot be [0,0]",
        ),
        (
            "rvmanin M[23]; show rvManinPath[M]{[1000000000001,1],[1,0]};",
            "absolute value at most 10^12",
        ),
        (
            "rvmanin M[23]; show rvManinHecke[2](rvManinBasis[M,3,0]);",
            "must be ≤ 31 and coprime to the level",
        ),
        (
            "rvmanin M[23]; show rvManinAtkin[3](rvManinBasis[M,3,0]);",
            "index must be 2 or p",
        ),
        (
            "rvmanin M[23]; show rvManinBasis[M,4,0];",
            "Manin sector must be 0",
        ),
    ] {
        let result = diagnostic(interpret(source));
        assert!(result.contains(expected), "{result}");
    }
}
