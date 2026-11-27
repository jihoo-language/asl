mod common;

use common::{diagnostic, example, interpret, stdout};

const DEMO: &str = "fricke-geometric-eichler-shimura";

#[test]
fn full_cuspidal_zeta_is_reciprocal_and_has_exact_jacobian_orders() {
    let out = stdout(example(DEMO));
    assert!(out.contains("X0(46) at ℓ=5, extension degree 1, q=5, genus 5"));
    assert!(
        out.contains("zeta numerator (ascending) = [1, 0, 5, 0, 30, -64, 150, 0, 625, 0, 3125]")
    );
    assert!(out.contains("zeta denominator = (1-T)(1-5T)"));
    assert!(out.contains("Frobenius trace = 0; #X0 = 6; #Jac = 3872"));
    assert!(out.contains("Frobenius trace = -10; #X0 = 36; #Jac = 15488000"));
    assert!(out.contains("X0(82) at ℓ=5, extension degree 1, q=5, genus 9"));
    assert!(out.contains("zeta numerator (ascending) = [1, 6, 41, 160, 680, 2088, 6864, 17536, 47570, 103204, 237850, 438400, 858000, 1305000, 2125000, 2500000, 3203125, 2343750, 1953125]"));
    assert!(out.contains("Frobenius trace = -6; #X0 = 12; #Jac = 15142400"));
    assert!(out.contains("Frobenius trace = -46; #X0 = 72; #Jac = 26076908748800"));
    assert!(out.contains("Frobenius trace = -36481726; #X0 = 95367468122352; #Jac = 652530696418020641886732909194681341324272040885850765796917909916035769493043357973424201872050000526517208906163290112000000"));
}

#[test]
fn independent_geometric_enumeration_matches_every_prime_or_square_fixture() {
    let out = stdout(example(DEMO));
    for (level, field, count) in [
        (46, 5, 6),
        (46, 25, 36),
        (46, 49, 80),
        (46, 121, 172),
        (82, 5, 12),
        (82, 25, 72),
        (82, 49, 108),
        (82, 121, 192),
        (10, 7, 8),
    ] {
        let prefix = format!("X0({level}) over F{field}: {count} points (4 cusps, {} noncuspidal); Eichler–Shimura count agrees",count-4);
        assert!(out.contains(&prefix), "{prefix}");
    }
    assert!(out.contains(
        "F25 elliptic model: trace = -10; #E = 36; rational 2-lines = 3; rational 23-lines = 24"
    ));
    assert!(out.contains(
        "F49 elliptic model: trace = -9; #E = 59; rational 2-lines = 0; rational 23-lines = 1"
    ));
    assert!(out.contains("(trace, 2-lines, 23-lines, models) ="));
    assert!(out.contains("Frobenius trace = 0; #X0 = 50; #Jac = 1"));
}

#[test]
fn exposes_typed_polynomials_and_supports_good_characteristic_three_analytically() {
    let out = stdout(interpret(
        "rvmanin M[23];
         let F=rvManinFrobPoly[M,5];
         let b=add(rvOre[0,0],add(rvOre[1,0,2],add(rvOre[2,0,6],add(rvOre[3,0,10],rvOre[4,0,25]))));
         let c=add(rvOre[0,0],add(rvOre[1,0,-4],rvOre[2,0,5]));
         assert F==rvOreMul(c,rvOreMul(b,b));
         rvManinZeta[M,3];",
    ));
    assert!(out.contains("X0(46) at ℓ=3, extension degree 1, q=3, genus 5"));
    assert!(out.contains("#X0 = 4;"));
    let wide = diagnostic(interpret("rvmanin M[41]; show rvManinFrobPoly[M,5,20];"));
    assert!(wide.contains("exceeds RVOre's signed coefficient range"));
}

#[test]
fn geometry_checks_bad_places_precision_and_field_size() {
    for (source, expected) in [
        ("rvmanin M[23]; rvManinZeta[M,23];", "coprime to 2p"),
        ("rvmanin M[23]; rvManinZeta[M,9];", "9 is not prime"),
        ("rvmanin M[23]; rvManinZeta[M,5,0];", "extension degree must lie in 1..=64"),
        ("rvmanin M[23]; rvunram F[23,1]{0,1}; rvManinPointCount[M,F];", "coprime to 2p"),
        ("rvmanin M[23]; rvunram F[3,1]{0,1}; rvManinPointCount[M,F];", "characteristic greater than three"),
        ("rvmanin M[23]; rvunram F[5,2]{0,1}; rvManinPointCount[M,F];", "at precision one"),
        ("rvmanin M[23]; rvunram F[5,1]{1,1,0,1}; rvManinPointCount[M,F];", "requires q ≤ 121"),
        ("rvmanin M[23]; rvunram F[5,1]{0,1}; rvManinCurveLines[M](rvUnramElem[F]{0},rvUnramElem[F]{0});", "require a nonsingular short Weierstrass model"),
        ("rvmanin M[23]; rvunram F[5,1]{0,1}; rvunram G[5,1]{0,1}; rvManinCurveLines[M](rvUnramElem[F]{0},rvUnramElem[G]{1});", "rvManinCurveLines"),
    ] {
        let result = diagnostic(interpret(source));
        assert!(result.contains(expected), "{result}");
    }
}

#[test]
fn nonzero_discriminant_divisible_by_p_squared_uses_the_actual_torsion_scheme() {
    let out = stdout(example(DEMO));
    assert!(out.contains(
        "F31 elliptic model: trace = 7; #E = 25; rational 2-lines = 0; rational 5-lines = 6"
    ));
    assert!(out.contains(
        "F31 elliptic model: trace = 7; #E = 25; rational 2-lines = 0; rational 5-lines = 1"
    ));
    assert!(out.contains(
        "X0(10) over F31: 32 points (4 cusps, 28 noncuspidal); Eichler–Shimura count agrees"
    ));
    assert!(out.contains("(7, 0, 1, 30), (7, 0, 6, 5)"));
    let alternative = stdout(interpret(
        "rvmanin M[5]; rvunram F[29,1]{0,1}; rvManinPointCount[M,F];",
    ));
    assert!(alternative.contains("X0(10) over F29: 30 points"));
}
