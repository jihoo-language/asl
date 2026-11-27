mod common;

use common::{diagnostic, example, interpret, stdout};

const F41: &str = include_str!("../examples/elliptic-cm-torelli.asl");

#[test]
fn all_twelve_principal_quotients_have_the_independent_torelli_profile() {
    let output = stdout(example("elliptic-cm-torelli"));
    for expected in [
        "{(1, 20, Some([1, 1])): 2, (2, 80, Some([1, 2])): 4, (3, 180, None): 2, (6, 720, None): 4}",
        "rank=4, signature=(1,3), determinant=-80",
        "rank=2, signature=(2,0), determinant=80",
        "rank=4, signature=(1,3), determinant=-720",
        "rank=2, signature=(2,0), determinant=720",
        "factor_polarization_degrees=[1, 1]",
        "splitting=indecomposable, geometric_complex_type=genus_two_jacobian",
        "quotient_definition_degree=2, fixed_factors=[0, 1]",
        "multiplier=1681",
        "rvCMQuotientWedgePair = 0",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn algebraic_riemann_form_is_identity_in_the_integral_rosati_order() {
    let output = stdout(interpret(&format!(
        "{F41}
        let e=rvCMThetaElem[K6]{{0,0,1,0}};
        rvcmquotient Full(rvCMThetaSubgroup[K6]{{e}});
        let h=rvCMQuotientWedge[Full]{{1,-1,-1,-1,-3,3}};
        assert rvCMQuotientToRosati(h)==rvCMQuotientEndIdentity[Full];
        rvCMQuotientWedgePair(h,h);
        rvCMQuotientWedgePair(h,rvCMQuotientTransBasis[Full,1]);
        rvCMQuotientTorelli[Full];"
    )));
    assert!(output.contains("rvCMQuotientWedgePair = 2"));
    assert!(output.contains("rvCMQuotientWedgePair = 0"));
    assert!(output.contains("factor_cm_conductors=[1, 1]"));
}

#[test]
fn nonprincipal_polarization_still_has_saturated_hodge_lattices() {
    let output = stdout(interpret(&format!(
        "{F41}
        rvcmquotient Original(rvCMThetaSubgroup[K6]{{}});
        rvCMQuotientNS[Original];
        rvCMQuotientTrans[Original];
        let a=rvCMQuotientNSBasis[Original,3];
        let b=rvCMQuotientTransBasis[Original,1];
        rvCMQuotientWedgePair(a,b);"
    )));
    assert!(
        output.contains("signature=(1,3), determinant=-20"),
        "{output}"
    );
    assert!(
        output.contains("signature=(2,0), determinant=20"),
        "{output}"
    );
    assert!(output.contains("rvCMQuotientWedgePair = 0"), "{output}");
}

#[test]
fn wedge_pairing_uses_full_wide_products_and_checks_nominality() {
    let output = stdout(interpret(&format!(
        "{F41}
        let big=rvCMQuotientWedge[Split]{{1000000000000000000000000,1000000000000000000000000,0,0,1000000000000000000000000,1000000000000000000000000}};
        assert add(big,neg(big))==zero[RVCMQuotientWedge(Split)];
        rvCMQuotientWedgePair(big,big);
        rvCMQuotientWedgePair(rvCMQuotientWedge[Split]{{1,0,0,0,0,0}},rvCMQuotientWedge[Split]{{0,0,0,0,0,1}});"
    )));
    assert!(
        output.ends_with(
            "rvCMQuotientWedgePair = 0
rvCMQuotientWedgePair = 1
"
        ),
        "{output}"
    );
}

#[test]
fn requires_principal_and_algebraic_input_for_torelli_and_rosati() {
    for (query, expected) in [
        (
            "rvcmquotient Original(rvCMThetaSubgroup[K6]{}); rvCMQuotientTorelli[Original];",
            "requires a principal polarization",
        ),
        (
            "show rvCMQuotientToRosati(tr);",
            "requires an algebraic wedge form",
        ),
        (
            "show rvCMQuotientProjector[Jacobian,0];",
            "splitting projector does not exist",
        ),
        (
            "show rvCMQuotientProjector[Split,2];",
            "splitting projector does not exist",
        ),
        (
            "rvCMQuotientFactorFrob[Jacobian];",
            "factor Frobenius requires a polarized elliptic product",
        ),
        (
            "show rvCMQuotientTransBasis[Split,2];",
            "wedge basis index is out of range",
        ),
        (
            "rvCMQuotientWedgePair(ns,rvCMQuotientNSBasis[Jacobian,0]);",
            "CM quotient wedge pairing argument",
        ),
        (
            "show rvCMQuotientToRosati(point);",
            "expected a CM quotient wedge form",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F41}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn odd_fundamental_discriminant_uses_the_actual_nonzero_cm_trace() {
    let output = stdout(interpret(
        r"
        rvunram k[7,1]{0,1};
        rvec E[k](rvUnramElem[k]{0},rvUnramElem[k]{3});
        rvec F[k](rvUnramElem[k]{0},rvUnramElem[k]{3});
        let f=rvECChart[E,F](rvUnramElem[k]{1});
        let g=rvECChart[E,F](rvUnramElem[k]{2});
        let fd=rvECChart[F,E](rvUnramElem[k]{1});
        let gd=rvECChart[F,E](rvUnramElem[k]{4});
        rvcmhomframe H(f,g);
        rvcmhomdual HD[H](fd,gd);
        rvcmhomframe RE(rvECMultMap[E,1],rvECMapCompose(fd,g));
        rvcmhomframe RF(rvECMultMap[F,1],rvECMapCompose(g,fd));
        rvcmabelian A[H,HD,RE,RF];
        rvcmthetakernel K(rvCMAbelianDiv[A,1,1](zero[RVCMHomCoordinate(H)]));
        rvcmquotient Q(rvCMThetaSubgroup[K]{});
        rvCMQuotientNS[Q];
        rvCMQuotientTrans[Q];
        rvCMQuotientTorelli[Q];",
    ));
    assert!(
        output.contains("signature=(1,3), determinant=-3"),
        "{output}"
    );
    assert!(output.contains("oriented_form=[1, 1, 1]"), "{output}");
    assert!(output.contains("factor_cm_conductors=[1, 1]"), "{output}");
}
