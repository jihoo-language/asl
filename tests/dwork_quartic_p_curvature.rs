mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn generic_rational_matrices_are_nonzero_nilpotent_and_certified() {
    let out = stdout(interpret(
        "rvDworkPCurvature[3,2]; rvDworkPCurvature[3,3];
         rvDworkPCurvature[11,2]; rvDworkPCurvature[29,3];",
    ));
    let rows: Vec<_> = out.lines().collect();
    assert!(rows[0].contains("denominator=(1-256z)^3, numerator=[[[], [1, 1, 1]], [[], []]], generic power ranks=[2, 1, 0]"));
    assert!(rows[1].contains("numerator=[[[], [1, 1, 1], []], [[], [0, 1, 2], [1, 1, 1]], [[], [0, 0, 2], [0, 2, 1]]], generic power ranks=[3, 2, 1, 0]"));
    assert!(rows[2].contains("Hasse=[1, 2, 1]"));
    assert!(rows.iter().all(|row| row.ends_with(
        "certified nilpotence, horizontality, skew-adjointness, symmetric square and Hasse section"
    )));
}

#[test]
fn prime_fiber_kernels_and_hasse_contacts_vary_with_the_actual_point() {
    let out = stdout(example("dwork-quartic-p-curvature"));
    let at_one = out
        .lines()
        .find(|row| row.starts_with("rvDworkPCurvatureFiber[11, 3, 1]"))
        .unwrap();
    assert!(at_one.contains(
        "matrix=[[7, 4, 0], [2, 3, 4], [1, 3, 1]], power ranks=[3, 2, 1, 0], nilpotency index=3"
    ));
    assert!(at_one.contains("kernel bases of Psi^1 through Psi^3=[[[8, 8, 1]],"));
    assert!(at_one.ends_with("Hasse jet=[4, 4, 6], Hasse order=0, Hasse line=ker(Psi), coordinate contact=nonzero first coordinate"));
    let double = out
        .lines()
        .find(|row| row.starts_with("rvDworkPCurvatureFiber[11, 3, -1]"))
        .unwrap();
    assert!(double.contains("power ranks=[3, 2, 1, 0]"));
    assert!(double
        .contains("kernel bases of Psi^1 through Psi^3=[[[0, 0, 1]], [[0, 1, 0], [0, 0, 1]],"));
    assert!(double.ends_with("Hasse jet=[0, 0, 2], Hasse order=2, Hasse line=ker(Psi), coordinate contact=third coordinate line"));
    let other = out
        .lines()
        .find(|row| row.starts_with("rvDworkPCurvatureFiber[13, 3, 9]"))
        .unwrap();
    assert!(other.contains("matrix=[[0, 0, 0], [9, 0, 0], [3, 9, 0]], power ranks=[3, 2, 1, 0]"));
    assert!(other.ends_with("Hasse line=ker(Psi), coordinate contact=third coordinate line"));
}

#[test]
fn prime_field_values_interoperate_with_the_hasse_witt_matrix_algebra() {
    stdout(interpret(
        "let C=rvDworkPCurvatureAt[11,2,1];
         assert C == rvHWMat[11,2,1]{9,4,10,2};
         assert rvHWCompose(C,C) == rvDworkPCurvatureAt[11,2,1,2];
         assert rvDworkPCurvatureAt[11,2,1,0] == rvHWMat[11,2,1]{1,0,0,1};
         assert rvDworkPCurvatureAt[11,2,10] == rvHWMat[11,2,1]{0,0,10,0};
         assert rvDworkPCurvatureAt[3,3,2] == rvHWMat[3,3,1]{0,2,0,0,2,2,0,1,1};",
    ));
}

#[test]
fn unramified_residue_specialization_uses_the_p_th_power_denominator() {
    let out = stdout(example("dwork-quartic-p-curvature"));
    let row = out
        .lines()
        .find(|row| row.starts_with("rvDworkPCurvatureClosedFiber"))
        .unwrap();
    assert!(row.contains("p=3, z=[0, 1], matrix=[[[0, 0], [2, 2], [0, 0]], [[0, 0], [1, 0], [2, 2]], [[0, 0], [2, 1], [2, 0]]]"));
    assert!(row.contains("power ranks=[3, 2, 1, 0], nilpotency index=3"));
    assert!(row.ends_with("Hasse line=ker(Psi), coordinate contact=nonzero first coordinate"));

    stdout(interpret(
        "rvunram K[3,1]{1,0,1};
         let v=rvDworkPCurvatureClosed[2](rvUnramElem[K]{0,1});
         assert rvHWWEntry[1,2](v) == rvUnramElem[K]{2,2};
         assert rvHWWEntry[2,2](v) == rvUnramElem[K]{};
         let vv=rvDworkPCurvatureClosed[3](rvUnramElem[K]{2});
         assert rvHWWDescend(vv) == rvDworkPCurvatureAt[3,3,2];",
    ));
}

#[test]
fn excludes_singular_coordinates_and_wrong_ring_or_rank() {
    for (source, message) in [
        ("rvDworkPCurvature[2,3];", "requires an odd prime"),
        ("rvDworkPCurvature[15,3];", "15 is not prime"),
        ("rvDworkPCurvature[263,3];", "p at most 257"),
        ("rvDworkPCurvature[11,4];", "rank must be 2 or 3"),
        ("show rvDworkPCurvatureAt[11,3,0];", "z nonzero"),
        ("rvDworkPCurvatureFiber[11,3,4];", "1-256z=0"),
        (
            "rvunram K[11,2]{0,1}; show rvDworkPCurvatureClosed[3](rvUnramElem[K]{1});",
            "requires a residue field of precision one",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(message), "{error}");
    }
}
