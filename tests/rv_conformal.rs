mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn reduced_vacua_and_insertion_cancellations() {
    assert_eq!(stdout(example("rv-conformal-series")),concat!(
        "rvConformalSeries = [1, -9] / [1, -15, 45]\n",
        "rvConformalSeries = [1, -39, 180] / [1, -49, 504, -1296]\n",
        "rvConformalMod = 11974497\nrvConformalMod = 278843085\n",
        "rvConformalNum = [1, -27]\nrvConformalDen = [1, -45, 324]\n",
        "rvConformalSeries = [] / [1]\nrvConformalMod = 0\n",
        "rvConformalExact = 10021138658838911696159984319316872623069878625683050878067576471599586795946379425004124641418457031250\n",
    ));
}

#[test]
fn large_genus_and_composite_moduli_need_no_inverse() {
    let maximum = u128::MAX;
    let program = format!(
        r"
        rvConformalMod[3,2,{maximum},{maximum},[]];
        rvConformalMod[3,2,{maximum},1,[]];
        rvConformalMod[3,0,{maximum},{maximum},[[2,2,2]]];
        rvConformalExact[3,0,{maximum},[]];
        rvConformalMod[2,2,1000000000000000000,1000000000,[[0,1],[0,1]]];
        rvConformalSeries[2,2,[[0,1],[0,1]]];
        rvConformalMod[3,3,0,25,[[0,0,1],[0,1,1]]];
        rvConformalExact[3,3,0,[[0,0,1],[0,1,1]]];
        rvConformalDim[3,2,5,[]];
        rvConformalExact[3,2,5,[]];
    "
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvConformalMod = 3599435871650236134940823810588855745\n",
            "rvConformalMod = 0\nrvConformalMod = 1\nrvConformalExact = 1\n",
            "rvConformalMod = 787109376\nrvConformalSeries = [1] / [1, -4]\n",
            "rvConformalMod = 1\nrvConformalExact = 1\n",
            "rvConformalDim = 42525\nrvConformalExact = 42525\n",
        )
    );
}

#[test]
fn modular_insertions_are_reduced_during_multiplication() {
    // The vacuum coefficient of adjoint^500 is F_499 at SU(3) level two.
    let labels = std::iter::repeat_n("[0,1,2]", 500)
        .collect::<Vec<_>>()
        .join(",");
    let program = format!("rvConformalMod[3,2,0,1000000000,[{labels}]];");
    assert_eq!(stdout(interpret(&program)), "rvConformalMod = 447204501\n");
}

#[test]
fn conformal_reports_validate_bounds_and_labels() {
    for (program, expected) in [
        ("rvConformalMod[3,2,7,0,[]];", "modulus must be positive"),
        (
            "rvConformalMod[3,2,340282366920938463463374607431768211456,1,[]];",
            "genus exceeds the u128 range",
        ),
        (
            "rvConformalMod[3,2,0,1,[[0,1,3]]];",
            "outside its level alcove",
        ),
        (
            "rvConformalNum[3,2,[[1,0,2]]];",
            "must be weakly increasing",
        ),
        (
            "rvConformalExact[3,2,1000000000000000000,[]];",
            "exceeds the decimal output limit",
        ),
        (
            "rvConformalDen[3,100,[]];",
            "limit of 4096 integrable weights",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
    assert_eq!(
        stdout(interpret(
            "rvConformalExact[3,2,1000000000000000000,[[0,0,1]]];"
        )),
        "rvConformalExact = 0\n"
    );
}
