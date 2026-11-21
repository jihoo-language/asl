mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn spherical_products_normalize_iwahori_counts_and_realize_hall_polynomials() {
    assert_eq!(
        stdout(example("rv-spherical-hecke")),
        concat!(
            "rvSphCoeff = [1]\n",
            "rvSphCoeff = [1, 1]\n",
            "rvSphAt = 5\n",
            "rvSphCoeff = [1, 1, 1]\n",
            "rvSphAt = 7\n",
            "rvSphCoeff = [0, 1, 1]\n",
            "rvSphAt = 72\n",
        ),
    );
}

#[test]
fn spherical_elements_handle_repeated_coweights_negative_shifts_and_signed_sums() {
    let program = r"
        let one = rvSphId[3];
        let a = rvSphBasis[[0,0,2]];
        let shifted = rvSphBasis[[-100000000000000000000,-100000000000000000000,-100000000000000000000]];
        let result = rvSphBasis[[-100000000000000000000,-100000000000000000000,-99999999999999999998]];
        assert rvSphMul(shifted,a) == result;
        assert rvSphMul(a,shifted) == result;
        let q = rvSphQ[3];
        let mixed = add(a,neg(q));
        assert rvSphMul(mixed,one) == mixed;
        assert rvSphMul(shifted,rvSphMul(mixed,a)) ==
               rvSphMul(rvSphMul(shifted,mixed),a);
        rvHeckeCoeff[[0,1,8]](rvSphIwahori(a));
        rvHeckeCoeff[[2,0,7]](rvSphIwahori(a));
        rvHeckeCoeff[[0,1,2]](rvSphIwahori(a));
        show mixed;
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvHeckeCoeff = [1]\n",
            "rvHeckeCoeff = [1]\n",
            "rvHeckeCoeff = [0]\n",
            "RVSphHecke(3){-q*S[0, 0, 0] + S[0, 0, 2]}\n",
        ),
    );
}

#[test]
fn spherical_rank_zero_and_one_have_the_expected_identities() {
    let program = r"
        assert rvSphId[0] == rvSphBasis[[]];
        assert rvSphMul(rvSphQ[0],rvSphId[0]) == rvSphQ[0];
        assert rvSphMul(rvSphBasis[[-7]],rvSphBasis[[9]]) == rvSphBasis[[2]];
        assert rvSphIwahori(rvSphBasis[[5]]) == rvHeckeT[[5]];
        rvSphAt[-3,[]](rvSphMul(rvSphQ[0],rvSphQ[0]));
    ";
    assert_eq!(stdout(interpret(program)), "rvSphAt = 9\n");
}

#[test]
fn spherical_bounds_reject_malformed_coweights_and_large_orbits() {
    for (program, expected) in [
        (
            "show rvSphBasis[[1,0]];",
            "the spherical coweight must be weakly increasing",
        ),
        (
            "rvSphCoeff[[0]](rvSphId[2]);",
            "coweight rank must match the algebra",
        ),
        ("show rvSphMul(rvSphId[2],rvSphId[3]);", "rvSphMul"),
        (
            "show rvSphIwahori(rvSphId[8]);",
            "spherical Hecke calculation exceeds the storage limit",
        ),
        (
            "show rvSphIwahori(rvHeckeId[2]);",
            "requires a spherical Hecke element",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
