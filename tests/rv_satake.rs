mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn normalized_hall_littlewood_images_and_satake_multiplication() {
    assert_eq!(
        stdout(example("rv-satake")),
        concat!(
            "rvSatakeCoeff = {1: 1}\n",
            "rvSatakeCoeff = {1: 1}\n",
            "rvSatakeCoeff = {0: -1, 2: 1}\n",
            "rvSatakeCoeff = {2: 1}\n",
            "rvSatakeCoeff = {0: 1}\n",
            "rvSatakeAt = 5/108\n",
            "rvSatakeAt = 5/2\n",
            "rvSatakeAt = 4960\n",
            "rvSatakeCoeff = {2: -1, 4: 1}\n",
        ),
    );
}

#[test]
fn satake_handles_signed_laurent_specialization_and_linearity() {
    let program = r"
        let a = rvSatake(rvSphBasis[[-1,-1]]);
        let b = rvSatake(rvSphBasis[[0,1]]);
        let lhs = rvSatake(add(rvSphBasis[[-1,-1]],neg(rvSphBasis[[0,1]])));
        assert lhs == add(a,neg(b));
        assert add(lhs,neg(lhs)) == zero[RVSatake(2)];
        rvSatakeAt[-2,[-2,3]](rvSatakeMul(a,rvSatV[2,-1]));
        rvSatakeCoeff[[9,9]](lhs);
        show b;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvSatakeAt = 1/12\nrvSatakeCoeff = {}\nRVSatake(2){v*X2 + v*X1}\n",
    );
}

#[test]
fn huge_pure_determinant_shifts_and_rank_zero_are_exact() {
    let program = r"
        let central = rvSatake(rvSphBasis[[10000000000000000000000000,10000000000000000000000000]]);
        rvSatakeCoeff[[10000000000000000000000000,10000000000000000000000000]](central);
        rvSatakeAt[2,[1,-1]](central);
        assert rvSatake(rvSphId[0]) == rvSatV[0,0];
        rvSatakeAt[-3,[]](rvSatake(rvSphQ[0]));
        rvSatakeAt[1,[1]](rvSatake(rvSphBasis[[-170141183460469231731687303715884105728]]));
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvSatakeCoeff = {0: 1}\nrvSatakeAt = 1\nrvSatakeAt = 9\nrvSatakeAt = 1\n",
    );
}

#[test]
fn satake_checks_rank_exponents_and_storage() {
    for (program, expected) in [
        (
            "show rvSatake(rvHeckeId[2]);",
            "rvSatake requires a spherical Hecke element",
        ),
        (
            "rvSatakeCoeff[[0]](rvSatV[2,0]);",
            "coordinates must match the polynomial rank",
        ),
        ("show rvSatakeMul(rvSatV[2,0],rvSatV[3,0]);", "rvSatakeMul"),
        (
            "rvSatakeAt[0,[1,1]](rvSatV[2,-1]);",
            "cannot specialize a negative Laurent power at zero",
        ),
        (
            "rvSatakeAt[2,[0,1]](rvSatake(rvSphBasis[[-1,-1]]));",
            "cannot specialize a negative Laurent power at zero",
        ),
        (
            "show rvSatake(rvSphBasis[[0,1000000]]);",
            "Satake Laurent calculation exceeds the storage limit",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
