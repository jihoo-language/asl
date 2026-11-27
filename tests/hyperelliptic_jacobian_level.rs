mod common;

use common::{diagnostic, example, interpret, stdout};

const THREE_TORSION: &str = r"
    rvunram F[5,1]{0,1};
    rvunram W[5,1]{3,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic E[F]{o,z,z,o};
    rvjacbase L[E,W];
    let zW=rvUnramElem[W]{};
    let oW=rvUnramElem[W]{1};
    let alpha=rvUnramElem[W]{0,1};
    let zeta=rvUnramElem[W]{2,2};
    let p=rvJacPoint[L](zW,oW);
    let q=rvJacPoint[L](oW,alpha);
    let ramified=rvJacPoint[L](neg(oW),zW);
";

#[test]
fn example_certifies_genus_two_and_three_torsion_frobenius() {
    assert_eq!(
        stdout(example("hyperelliptic-jacobian-level")),
        concat!(
            "rvJacTorsionFrob[B2, C, 0] = multiplier=1, matrix=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], characteristic=[1, 0, 0, 0, 1]\n",
            "rvJacTorsionFrob[B2, C, 1] = multiplier=1, matrix=[[1, 0, 0, 0], [0, 1, 0, 1], [0, 0, 1, 0], [0, 0, 0, 1]], characteristic=[1, 0, 0, 0, 1]\n",
            "rvJacTorsionFrob[B3, E, 1] = multiplier=2, matrix=[[1, 0], [0, 2]], characteristic=[1, 0, 2]\n",
        )
    );
}

#[test]
fn coordinates_can_use_both_sides_of_the_pairing_and_large_powers() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{THREE_TORSION}
            rvjactorsionbasis B[L,3](zeta){{p,q}};
            let t = add(neg(p),q);
            assert rvJacTorsionCoord[B,0](t) == rvPrimeElt[3,2];
            assert rvJacTorsionCoord[B,1](t) == rvPrimeElt[3,1];
            assert rvJacTorsionCoord[B,0](zero[RVJac(L)]) == zero[PrimeScalar(3)];
            rvJacTorsionFrob[B,E,1000001];
            rvJacTorsionFrob[B,L,1000001];
        "
        ))),
        concat!(
            "rvJacTorsionFrob[B, E, 1000001] = multiplier=2, matrix=[[1, 0], [0, 2]], characteristic=[1, 0, 2]\n",
            "rvJacTorsionFrob[B, L, 1000001] = multiplier=1, matrix=[[1, 0], [0, 1]], characteristic=[1, 1, 1]\n",
        )
    );
}

#[test]
fn rejects_invalid_certificates() {
    for (statement, expected) in [
        ("rvjactorsionbasis B[L,1](zeta){p,q};", "must be prime"),
        ("rvjactorsionbasis B[L,9](zeta){p,q};", "9 is not prime"),
        (
            "rvjactorsionbasis B[L,5](zeta){p,q};",
            "differ from the residue characteristic",
        ),
        (
            "rvjactorsionbasis B[L,3](zeta){p};",
            "requires exactly 2 generators",
        ),
        (
            "rvjactorsionbasis B[L,3](zeta){p,q,p};",
            "requires exactly 2 generators",
        ),
        (
            "rvjactorsionbasis B[L,3](oW){p,q};",
            "primitive root of unity of order 3",
        ),
        (
            "rvjactorsionbasis B[L,3](alpha){p,q};",
            "primitive root of unity of order 3",
        ),
        (
            "rvjactorsionbasis B[L,3](zeta){p,ramified};",
            "generator 1 must be killed by 3",
        ),
        (
            "rvjactorsionbasis B[L,3](zeta){p,p};",
            "incorrect symplectic pairing",
        ),
        (
            "rvjactorsionbasis B[L,3](zeta){q,p};",
            "incorrect symplectic pairing",
        ),
        (
            "rvunram U[5,1]{3,0,1}; rvjactorsionbasis B[L,3](rvUnramElem[U]{2,2}){p,q};",
            "expected type RVUnram(W)",
        ),
        (
            "rvjacbase N[E,W]; rvjactorsionbasis B[L,3](zeta){p,zero[RVJac(N)]};",
            "expected type RVJac(L)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{THREE_TORSION}{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn basis_references_remain_nominal_and_frobenius_checks_the_extension() {
    for (statement, expected) in [
        ("rvjactorsionbasis B[L,3](zeta){p,q};", "already defined"),
        (
            "show rvJacTorsionCoord[Unknown,0](p);",
            "unknown Jacobian torsion basis",
        ),
        (
            "show rvJacTorsionCoord[B,2](p);",
            "coordinate index is out of range",
        ),
        (
            "show rvJacTorsionCoord[B,0](ramified);",
            "requires a class killed by 3",
        ),
        (
            "rvjacbase N[E,W]; show rvJacTorsionCoord[B,0](zero[RVJac(N)]);",
            "expected type RVJac(L)",
        ),
        (
            "assert rvJacTorsionCoord[B,0](p) == rvPrimeElt[2,1];",
            "expected type PrimeScalar(3)",
        ),
        (
            "rvhyperelliptic N[F]{o,o,z,o}; rvJacTorsionFrob[B,N,0];",
            "not the declared base extension",
        ),
    ] {
        let error = diagnostic(interpret(&format!(
            "{THREE_TORSION}rvjactorsionbasis B[L,3](zeta){{p,q}};{statement}"
        )));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn genus_zero_has_an_empty_torsion_basis_and_matrix() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            rvhyperelliptic C[F]{z,o};
            rvjactorsionbasis B[C,2](neg(o)){};
            rvJacTorsionFrob[B,C,1];
        "
        )),
        "rvJacTorsionFrob[B, C, 1] = multiplier=1, matrix=[], characteristic=[1]\n"
    );
}
