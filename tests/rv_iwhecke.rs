mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn generic_hecke_relations_and_specialized_complete_flag_counts() {
    assert_eq!(
        stdout(example("rv-hecke")),
        concat!(
            "rvHeckeCoeff = [0, 1]\n",
            "rvHeckeCoeff = [-1, 1]\n",
            "rvHeckeAt = 2\n",
            "rvHeckeAt = 1\n",
            "rvHeckeAt = 4\n",
            "rvHeckeAt = 3\n",
            "rvHeckeCoeff = [0, 1, -1, 1]\n",
            "rvHeckeAt = 6\n",
            "rvHeckeAt = 52\n",
        ),
    );
}

#[test]
fn hecke_kl_coefficients_reuse_singular_and_affine_schubert_polynomials() {
    let program = r"
        let finite = rvHeckeKL[[2,3,0,1]];
        rvHeckeCoeff[[0,1,2,3]](finite);
        rvHeckeCoeff[[2,3,0,1]](finite);
        rvHeckeCoeff[[3,2,1,0]](finite);
        let affine = rvHeckeKL[[-4,3,4]];
        rvHeckeCoeff[[0,1,2]](affine);
        rvHeckeCoeff[[-4,3,4]](affine);
        assert add(affine, neg(affine)) == zero[RVIwHecke(3)];
        let translate = rvHeckeT[[1,2,3]];
        assert rvHeckeKL[[-3,4,5]] == rvHeckeMul(translate,affine);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvHeckeCoeff = [1, 1]\n",
            "rvHeckeCoeff = [1]\n",
            "rvHeckeCoeff = [0]\n",
            "rvHeckeCoeff = [1, 1]\n",
            "rvHeckeCoeff = [1]\n",
        ),
    );
}

#[test]
fn hecke_handles_small_ranks_signed_scalars_and_long_rotations() {
    let program = r"
        let q0 = rvHeckeQ[0];
        rvHeckeCoeff[[]](rvHeckeMul(q0,q0));
        rvHeckeAt[-3,[]](neg(rvHeckeMul(q0,q0)));
        assert rvHeckeId[0] == rvHeckeT[[]];
        let t = rvHeckeT[[7]];
        assert rvHeckeMul(t,rvHeckeT[[-9]]) == rvHeckeT[[-2]];
        let s = rvHeckeT[[1,0,2]];
        let omega = rvHeckeT[[100000000000000000000,100000000000000000001,100000000000000000002]];
        assert rvHeckeMul(s,omega) == rvHeckeT[[99999999999999999999,100000000000000000001,100000000000000000003]];
        show neg(s);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvHeckeCoeff = [0, 0, 1]\nrvHeckeAt = -9\nRVIwHecke(3){-T[1, 0, 2]}\n",
    );
}

#[test]
fn affine_hecke_checks_types_windows_and_fixed_width_integer_bounds() {
    for (program, message) in [
        ("show rvHeckeT[[0,2]];", "window residues must form a permutation"),
        ("rvHeckeCoeff[[0]](rvHeckeId[2]);", "window rank must match the algebra"),
        ("show rvHeckeMul(rvHeckeId[2],rvHeckeId[3]);", "rvHeckeMul"),
        ("show neg(One.pt);", "neg requires an affine Hecke element"),
        ("let q=rvHeckeQ[0]; rvHeckeAt[170141183460469231731687303715884105727,[]](rvHeckeMul(q,q));", "affine Hecke arithmetic is out of range"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(message), "{program}\n{error}");
    }
    let error = diagnostic(interpret(&format!("show rvHeckeId[{}];", usize::MAX)));
    assert!(
        error.contains("affine Hecke calculation exceeds the storage limit"),
        "{error}"
    );
}
