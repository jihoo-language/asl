mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/hyperelliptic-richelot-matrix.asl");

const CHARTED: &str = concat!(
    "rvJacRichelotMatrix[B3, H3] = level=3, root exponent=2, multiplier=1, rank=4, ",
    "column-action matrix=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], ",
    "column-action adjoint=[[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]\n",
    "rvJacRichelotMatrix[B3, K3] = level=3, root exponent=1, multiplier=2, rank=4, ",
    "column-action matrix=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]], ",
    "column-action adjoint=[[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]\n",
    "rvJacRichelotMatrix[B2, H2] = level=2, root exponent=1, multiplier=0, rank=2, ",
    "column-action matrix=[[0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0]], ",
    "column-action adjoint=[[0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0]]\n",
);

#[test]
fn charted_matrices_track_both_primitive_root_choices_and_the_two_torsion_kernel() {
    assert_eq!(stdout(example("hyperelliptic-richelot-matrix")), CHARTED);
}

#[test]
fn adjoint_has_the_correct_off_diagonal_sign() {
    let out = stdout(interpret(&format!(
        r"{EXAMPLE}
        rvjactorsionbasis P3[N,3](rvUnramPow[2](zeta)){{im0,im1,add(im2,im0),im3}};
        rvJacRichelotMatrix[B3,P3];
        "
    )));
    assert_eq!(
        out.strip_prefix(CHARTED).unwrap(),
        concat!(
            "rvJacRichelotMatrix[B3, P3] = level=3, root exponent=2, multiplier=1, rank=4, ",
            "column-action matrix=[[1, 0, 2, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], ",
            "column-action adjoint=[[2, 0, 2, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]\n",
        )
    );
}

#[test]
fn a_monic_quintic_target_resolves_without_a_chart() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram W[3,1]{1,0,1};
            let z=rvUnramElem[W]{};
            let o=rvUnramElem[W]{1};
            let m=neg(o);
            let i=rvUnramElem[W]{0,1};
            rvhyperelliptic C[W]{z,m,z,z,z,o};
            rvrichelot D[C]{[-1],[],[1]}{[1],[],[1]}{[],[1]};
            let t0=rvJacPoint[C](z,z);
            let t1=rvJacPoint[C](o,z);
            let t2=rvJacPoint[C](m,z);
            let t3=rvJacPoint[C](i,z);
            let e2=add(add(t0,t1),t2);
            let f2=add(add(t0,t1),t3);
            rvjactorsionbasis B[C,2](m){t0,e2,t1,f2};
            rvjactorsionbasis H[D,2](m){rvJacExtend[D](t0),rvJacExtend[D](e2),rvJacExtend[D](t1),rvJacExtend[D](f2)};
            rvJacRichelotMatrix[B,H];
        "
        )),
        concat!(
            "rvJacRichelotMatrix[B, H] = level=2, root exponent=1, multiplier=0, rank=2, ",
            "column-action matrix=[[0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0]], ",
            "column-action adjoint=[[0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0]]\n",
        )
    );
}

#[test]
fn rejects_wrong_models_levels_and_nominal_endpoints() {
    for (statement, expected) in [
        (
            "rvJacRichelotMatrix[B3,B3];",
            "has no stored Richelot splitting",
        ),
        (
            "rvJacRichelotMatrix[H3,H3];",
            "requires a basis of its named source curve 'C'",
        ),
        (
            "rvJacRichelotMatrix[B2,H3];",
            "must have the same prime level",
        ),
        (
            "rvhyperelliptic E[L]{z,one}; rvjactorsionbasis V[E,2](neg(one)){}; rvJacRichelotMatrix[V,H2];",
            "require genus-two bases",
        ),
        (
            "rvJacRichelotMatrix[B2,Unknown];",
            "unknown Jacobian torsion basis",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{EXAMPLE}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}
