mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn same_trace_has_different_kernel_counts_and_irreducible_kernels_are_effective() {
    let out = stdout(example("cyclic-velu"));
    for expected in [
        "rvECPointCount[C] = 19, trace=-5",
        "rvECPointCount[D] = 19, trace=-5",
        "rvECIsogenies[C, 3] = count=1, kernels=[[[2], [1]]]",
        "rvECIsogenies[D, 3] = count=4, kernels=[[[0], [1]], [[2], [1]], [[5], [1]], [[6], [1]]]",
        "RVECCurve(Q){a=[0], b=[11]}",
        "rvECCertify[I] = degree=3, cleared identity=yes",
        "rvECPointCount[E] = 5, trace=3",
        "rvECIsogenies[E, 5] = count=2, kernels=[[[1], [3], [1]], [[0], [5], [1]]]",
        "xnum=[[0], [1], [1], [1], [6], [1]], xden=[[1], [6], [4], [6], [1]]",
        "rvECCertify[J] = degree=5, cleared identity=yes",
        "rvECPointCount[A] = 28, trace=-2",
        "rvECCertify[L] = degree=3, cleared identity=yes",
    ] {
        assert!(
            out.contains(expected),
            "{expected}
{out}"
        );
    }
}

#[test]
fn an_actual_degree_seven_kernel_uses_the_orbit_and_maps_its_generator_to_infinity() {
    let result = stdout(interpret(
        r"
        rvunram F[11,1]{0,1};
        rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{1});
        let k = rveckernel[E,7]{[0],[7],[2],[1]};
        rvecisogeny V[E,B](k);
        let p=rvECPoint[E](rvUnramElem[F]{0},rvUnramElem[F]{1});
        assert rvECScale[7](p) == rvECInfinity[E];
        assert rvECApply(rvECIso[V],p) == rvECInfinity[B];
        assert rvECApply(rvECIso[V],rvECScale[-5](p)) == rvECInfinity[B];
        assert rvECXNum[7](rvECIso[V]) == rvUnramElem[F]{1};
        assert rvECXDen[6](rvECIso[V]) == rvUnramElem[F]{1};
        assert rvECKCoeff[2](rvECGetKernel(rvECIso[V])) == rvUnramElem[F]{2};
        assert rvECKCoeff[7](k) == zero[RVUnram(F)];
        rvECPointCount[E];
        rvECCertify[V];
        ",
    ));
    assert_eq!(
        result,
        concat!(
            "rvECPointCount[E] = 14, trace=-2
",
            "rvECCertify[V] = degree=7, cleared identity=yes
",
        )
    );
}

#[test]
fn a_divisor_of_psi_five_can_fail_the_generator_orbit_test() {
    let source = r"
        rvunram F[19,1]{0,1};
        rvec E[F](rvUnramElem[F]{1},zero[RVUnram(F)]);
        let first = rveckernel[E,5]{[7],[5],[1]};
        let second = rveckernel[E,5]{[7],[14],[1]};
    ";
    let good = stdout(interpret(source));
    assert!(good.is_empty());
    let failure = diagnostic(interpret(&format!(
        "{source} let mixed = rveckernel[E,5]{{[12],[4],[1]}};"
    )));
    assert!(
        failure.contains("not stable under the cyclic generator"),
        "{failure}"
    );
}

#[test]
fn supplied_kernels_do_not_require_field_or_kernel_enumeration() {
    let source = r"
        rvunram F[4294967291,1]{0,1};
        rvec E[F](zero[RVUnram(F)],rvUnramElem[F]{1});
        let k = rveckernel[E,3]{[],[1]};
        rvecisogeny V[E,B](k);
        assert rvECB[B] == rvUnramElem[F]{-27};
        assert rvECXNum[0](rvECIso[V]) == rvUnramElem[F]{4};
        rvECCertify[V];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvECCertify[V] = degree=3, cleared identity=yes
"
    );
    let failed = diagnostic(interpret(&format!("{source} rvECIsogenies[E,3];")));
    assert!(failed.contains("computation limit"), "{failed}");
}

#[test]
fn empty_kernel_types_and_actual_rational_point_domains_tabulate() {
    assert_eq!(
        stdout(interpret(
            r"
        rvunram W[5,1]{0,1};
        rvec E[W](rvUnramElem[W]{1},zero[RVUnram(W)]);
        card (\p : RVECPoint(E) => S.hi);
        card (\k : RVECKernel(E,3) => S.hi);
        rvECIsogenies[E,3];
        card (\e : RVECCurve(E) => S.hi);
    "
        )),
        concat!(
            "card = 4 / 4
",
            "card = 0 / 0
",
            "rvECIsogenies[E, 3] = count=0, kernels=[]
",
            "card = 1 / 1
",
        )
    );
}

#[test]
fn singular_input_inseparability_nondivisibility_and_point_errors_are_explicit() {
    let start = r"
        rvunram F[7,1]{0,1};
        rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{1});
    ";
    for (rest, expected) in [
        (
            "rvec Bad[F](zero[RVUnram(F)],zero[RVUnram(F)]);",
            "singular",
        ),
        ("rvECIsogenies[E,7];", "equals the residue characteristic"),
        ("rvECIsogenies[E,2];", "odd prime"),
        ("rvECIsogenies[E,9];", "not prime"),
        ("let a=rveckernel[E,5]{[1],[1]};", "monic of degree"),
        ("let a=rveckernel[E,5]{[0],[0],[1]};", "repeated roots"),
        ("let a=rveckernel[E,3]{[0],[1]};", "does not divide"),
        (
            "show rvECPoint[E](zero[RVUnram(F)],zero[RVUnram(F)]);",
            "does not satisfy",
        ),
        ("show rvECX(rvECInfinity[E]);", "no affine coordinate"),
        (
            "rvunram R[7,2]{0,1}; rvec Bad[R](rvUnramElem[R]{1},rvUnramElem[R]{1});",
            "precision one",
        ),
    ] {
        let message = diagnostic(interpret(&format!("{start}{rest}")));
        assert!(message.contains(expected), "{rest} -> {message}");
    }
}

#[test]
fn degree_eleven_larger_than_the_characteristic_preserves_nominal_endpoints() {
    let source = r"
        rvunram F[7,1]{0,1};
        rvec E[F](rvUnramElem[F]{1},rvUnramElem[F]{6});
        let k=rveckernel[E,11]{[3],[2],[6],[4],[5],[1]};
        rvecisogeny V[E,B](k);
        let p=rvECPoint[E](rvUnramElem[F]{1},rvUnramElem[F]{1});
        assert rvECScale[11](p) == rvECInfinity[E];
        assert rvECApply(rvECIso[V],p) == rvECInfinity[B];
        assert rvECA[B] == rvECA[E];
        assert rvECB[B] == rvECB[E];
        assert rvECXNum[11](rvECIso[V]) == rvUnramElem[F]{1};
        rvECCertify[V];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvECCertify[V] = degree=11, cleared identity=yes\n"
    );
    for (extra, expected) in [
        ("show rvECApply(rvECIso[V],rvECInfinity[B]);", "rvECApply"),
        (
            "rvecisogeny Wrong[B,Q](k);",
            "certified kernel on the source curve",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{source}{extra}")));
        assert!(error.contains(expected), "{error}");
    }
}
