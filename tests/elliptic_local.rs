mod common;

use common::{diagnostic, example, interpret, stdout};

const E61: &str = "rvunram W[61,1]{0,1}; rvec E[W](rvUnramElem[W]{0},rvUnramElem[W]{2});";

#[test]
fn the_actual_nine_torsion_certificate_detects_the_sign_and_the_velu_edge_directions() {
    let output = stdout(example("ordinary-tate-orders"));
    let mut lines = output.lines();
    let certificate = lines.next().unwrap();
    assert!(
        certificate.starts_with(
            "rvEllFrobCert[E, 3, 2] = m=9, a=4, representative=4, deg(psi)=40, x=true, y=false"
        ),
        "{certificate}"
    );
    assert!(certificate.contains("frob_x=["));
    assert!(certificate.contains("frob_y=["));
    assert!(
        certificate.contains("residual_x=[], residual_y=[("),
        "{certificate}"
    );
    assert_eq!(lines.collect::<Vec<_>>(), [
        "rvEllLocalOrder[E, 3] = r=2, bound=2, m=9, a=5, trace=-1, norm=1, discriminant=-3",
        "rvEllLocalOrder[Lower, 3] = r=0, bound=1, m=1, a=0, trace=-5, norm=13, discriminant=-27",
        "rvEllLocalOrder[Upper, 3] = r=1, bound=1, m=3, a=2, trace=-3, norm=3, discriminant=-3",
        "rvEllEdge[Up] = ascending, r_source=0, r_target=1",
        "rvEllEdge[Flat] = horizontal, r_source=1, r_target=1",
        "rvEllEdge[DownA] = descending, r_source=1, r_target=0",
        "rvEllEdge[DownB] = descending, r_source=1, r_target=0",
        "rvEllEdge[DownC] = descending, r_source=1, r_target=0",
    ]);
}

#[test]
fn coefficients_in_a_proper_extension_and_a_second_nine_torsion_curve_are_certified() {
    let program = r"
        rvunram W[5,1]{2,0,1};
        rvec E[W](rvUnramElem[W]{0,1},rvUnramElem[W]{0});
        assert rvEllFrobCoordinates[E,3,0] == (S.hi,S.hi);
        assert rvEllFrobCoordinates[E,3,1] == (S.hi,S.hi);
        assert rvEllFrobCoordinates[E,3,1,-1] == (S.hi,S.lo);
        rvEllLocalOrder[E,3];
        rvunram V[73,1]{0,1};
        rvec F[V](rvUnramElem[V]{0},rvUnramElem[V]{2});
        assert rvEllFrobCoordinates[F,3,2] == (S.hi,S.hi);
        assert rvEllFrobCoordinates[F,3,2,-1] == (S.hi,S.lo);
        rvEllLocalOrder[F,3];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvEllLocalOrder[E, 3] = r=1, bound=1, m=3, a=1, trace=2, norm=2, discriminant=-4\n",
            "rvEllLocalOrder[F, 3] = r=2, bound=2, m=9, a=1, trace=-1, norm=1, discriminant=-3\n",
        )
    );
    let lift = stdout(interpret(
        r"
        rvunram W[5,1]{2,0,1};
        rvec E[W](rvUnramElem[W]{0,1},rvUnramElem[W]{0});
        rvEllFrobCert[E,11,1,5];
    ",
    ));
    assert!(
        lift.contains("m=11, a=5, representative=-6, deg(psi)=60"),
        "{lift}"
    );
}

#[test]
fn the_full_twenty_seven_torsion_scheme_is_tested_at_its_existing_budget() {
    let program = r"
        rvunram W[547,1]{0,1};
        rvec E[W](rvUnramElem[W]{0},rvUnramElem[W]{2});
        rvEllLocalOrder[E,3];
        assert rvEllFrobCoordinates[E,3,3,14] == (S.hi,S.hi);
        assert rvEllFrobCoordinates[E,3,3,-14] == (S.hi,S.lo);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvEllLocalOrder[E, 3] = r=3, bound=3, m=27, a=14, trace=-1, norm=1, discriminant=-3\n"
    );
}

#[test]
fn torsion_and_local_order_require_the_right_prime_scalar_and_cm_field() {
    for (source,expected) in [
        ("rvEllLocalOrder[E,2];","must be an odd prime"),
        ("rvEllFrobCert[E,61,1];","equals the residue characteristic"),
        ("rvEllFrobCert[E,3,2,3];","scalar candidate must be prime to the torsion order"),
        ("rvEllFrobCert[E,3,90,1];","torsion order ell^k exceeds u64"),
        ("rvunram V[5,1]{0,1}; rvec F[V](rvUnramElem[V]{0},rvUnramElem[V]{1}); rvEllLocalOrder[F,3];", "requires an ordinary curve"),
        ("rvunram V[73,1]{0,1}; rvec F[V](rvUnramElem[V]{0},rvUnramElem[V]{2}); rvellipticcm K[E]; show rvEllLocalGen[K,F,3];", "requires the curve's Frobenius field"),
    ] {
        let error=diagnostic(interpret(&format!("{E61} {source}")));
        assert!(error.contains(expected),"{source}\n{error}");
    }
}
