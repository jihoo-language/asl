mod common;

use common::{diagnostic, example, interpret, stdout};

const NONSPLIT_F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,o,o,o,z,o};
";

#[test]
fn nonsplit_quartic_sign_choices_and_ramified_point_in_characteristic_three() {
    assert_eq!(
        stdout(example("hyperelliptic-jacobian-halving")),
        concat!(
            "RVJac(C){u=[[2], [0], [1]], v=[[2]]}\n",
            "RVJac(C){u=[[2], [2], [1]], v=[[1], [1]]}\n",
            "RVJac(C){u=[[2], [0], [1]], v=[[0], [2]]}\n",
        ),
    );
}

#[test]
fn verifies_affine_point_and_witness_over_an_extension_field() {
    // Translating the F3 quintic by alpha gives genuine extension coefficients
    // in the curve, the witness, and the output divisor.
    assert_eq!(
        stdout(interpret(
            r"
        rvunram F[3,1]{1,0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let a=rvUnramElem[F]{0,1};
        let b=add(a,neg(o));
        let c=neg(b);
        rvhyperelliptic C[F]{b,neg(a),c,z,neg(a),o};
        let h=rvJacHalfPoint[C](c,o){[-1,1],[-1],[],[1]};
        assert h==rvjac[C]{[1],[0,-1],[1]}{[2]};
        assert add(h,h)==rvJacPoint[C](c,o);
        ",
        )),
        ""
    );
}

#[test]
fn genus_zero_reduces_the_witness_and_accepts_the_zero_square_root() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{NONSPLIT_F3}
        rvhyperelliptic L[F]{{o,o}};
        assert rvJacHalfPoint[L](z,t){{[2],[1]}}==zero[RVJac(L)];
        assert rvJacHalfPoint[L](t,z){{}}==zero[RVJac(L)];
        "
        ))),
        ""
    );
}

#[test]
fn rejects_invalid_points_coordinates_and_witnesses() {
    for (source, expected) in [
        (
            "show rvJacHalfPoint[C](o,z){};",
            "an affine point on the curve",
        ),
        ("show rvJacHalfPoint[C](o,o){};", "witness must square to"),
        (
            "show rvJacHalfPoint[C](o,o){[-2],[-2],[],[-1]};",
            "witness norm must equal",
        ),
        (
            "show rvJacHalfPoint[C](o,o){[1,0]};",
            "at most 1 coordinates",
        ),
        (
            "rvunram G[3,1]{0,1}; show rvJacHalfPoint[C](rvUnramElem[G]{1},o){};",
            "rvJacHalfPoint coordinate",
        ),
        (
            "rvhyperelliptic L[F]{o,t}; show rvJacHalfPoint[L](z,o){};",
            "monic odd-degree",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{NONSPLIT_F3}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
