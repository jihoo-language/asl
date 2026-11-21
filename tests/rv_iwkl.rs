mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn finite_and_affine_schubert_cells_have_nonconstant_kl_polynomials() {
    assert_eq!(
        stdout(example("rv-iwkl")),
        concat!(
            "rvIwKL = [1, 1]\n",
            "rvIwKL = [1, 1]\n",
            "rvIwKL = [0]\n",
            "rvIw = [2, 3, 0, 1]\n",
            "rvIwLength = 4\n",
            "rvIwSchubertKL = [1, 1]\n",
            "rvIwSchubertKL = [1]\n",
            "rvIwSchubertKL = [1]\n",
        ),
    );
}

#[test]
fn kl_normalizes_length_zero_translates_and_distinguishes_components() {
    let shift = 100_000_000_000_000_000_000_i128;
    let x = [shift, shift + 1, shift + 2, shift + 3];
    let w = [shift + 2, shift + 3, shift, shift + 1];
    let lhs = [shift - 1, shift, shift + 1, shift + 2];
    let program =
        format!("rvIwKL[{x:?},{w:?}];\nrvIwKL[{lhs:?},{w:?}];\nrvIwKL[[],[]];\nrvIwKL[[7],[8]];");
    assert_eq!(
        stdout(interpret(&program)),
        "rvIwKL = [1, 1]\nrvIwKL = [0]\nrvIwKL = [1]\nrvIwKL = [0]\n",
    );
}

#[test]
fn kl_checks_windows_flag_invertibility_and_computation_limits() {
    for (program, expected) in [
        ("rvIwKL[[0],[0,1]];", "windows must have the same rank"),
        (
            "rvIwKL[[0,2],[0,1]];",
            "window residues must form a permutation",
        ),
        (
            "rvIwSchubertKL[[0,1]](rvU[0](rvId[One,[0],3,[1]]));",
            "one entry per Margolis dimension",
        ),
        (
            "rvIwSchubertKL[[0]](zero[RVLaurent(One,One,[0],3,[1])]);",
            "requires invertible Laurent endomorphisms",
        ),
        (
            "rvIwKL[[0,1],[-100000000,100000001]];",
            "affine Kazhdan–Lusztig calculation exceeds the storage limit",
        ),
    ] {
        let failure = diagnostic(interpret(program));
        assert!(failure.contains(expected), "{program}\n{failure}");
    }
}
