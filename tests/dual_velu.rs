mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn coprime_diamonds_match_both_rational_coordinates_at_common_endpoints() {
    assert_eq!(
        stdout(example("dual-velu")),
        concat!(
            "rvECMapCertify = degree=3, cleared identity=yes\n",
            "rvECMapCertify = degree=5, cleared identity=yes\n",
            "RVECCurve(Q35){a=[3], b=[6]}\n",
            "rvECMapDegree = 15\n",
            "rvECMapCertify = degree=15, cleared identity=yes\n",
            "rvECMapDegree = 15\n",
            "RVECCurve(P35){a=[1], b=[7]}\n",
        )
    );
}

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let one=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    rvec E[F](one,one);
    let k=rveckernel[E,5]{[1],[3],[1]};
    rvecisogeny V[E,B](k);
    let iso=rvECIso[V];
    let f=rvECMap(iso);
    let d=rvECDual(f);
";

#[test]
fn irreducible_kernel_duals_and_their_twists_have_the_right_scale() {
    let source = format!(
        r"{F7}
        assert rvECMapKCoeff[0](d) == zero[RVUnram(F)];
        assert rvECMapKCoeff[1](d) == rvUnramElem[F]{{6}};
        assert rvECMapKCoeff[2](d) == one;
        assert rvECMapXNum[0](d) == one;
        assert rvECMapYNum[0](d) == rvUnramElem[F]{{6}};
        assert rvECMapCompose(d,f) == rvECMultMap[E,5];
        assert rvECMapCompose(f,d) == rvECMultMap[B,5];
        assert rvECDual(d) == f;
        let u=rvECChart[E,B](two);
        let v=rvECChart[B,E](rvUnramElem[F]{{4}});
        assert rvECMapCompose(v,u) == rvECId[E];
        let p=rvECPoint[E](two,two);
        assert rvECMapApply(f,p) == rvECApply(iso,p);
        assert rvECMapApply(rvECMapCompose(d,f),p) == rvECInfinity[E];
        rvECMapCertify(d);

        rvec Twist[F](two,rvUnramElem[F]{{6}});
        let kt=rveckernel[Twist,5]{{[2],[2],[1]}};
        rvecisogeny P[Twist,T](kt);
        let twist=rvECMap(rvECIso[P]);
        let dt=rvECDual(twist);
        assert rvECMapKCoeff[1](dt) == rvUnramElem[F]{{4}};
        assert rvECMapCompose(dt,twist) == rvECMultMap[Twist,5];
        assert rvECMapCompose(twist,dt) == rvECMultMap[T,5];
        assert rvECDual(dt) == twist;
        rvECMapCertify(dt);
        "
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvECMapCertify = degree=5, cleared identity=yes\n".repeat(2)
    );
}

#[test]
fn signed_multiplication_maps_include_the_inseparable_case() {
    let source = format!(
        r"{F7}
        let p=rvECPoint[E](two,two);
        let a=rvECMultMap[E,2];
        let b=rvECMultMap[E,3];
        assert rvECMapCompose(a,b) == rvECMultMap[E,6];
        assert rvECMapCompose(a,rvECMultMap[E,-1]) == rvECMultMap[E,-2];
        assert rvECMapApply(rvECMultMap[E,-4],p) == rvECScale[-4](p);
        assert zero[RVECMap(E,E)] == rvECMultMap[E,0];
        assert rvECMapCompose(a,rvECMultMap[E,0]) == rvECMultMap[E,0];
        assert rvECMapCompose(rvECMultMap[E,0],a) == rvECMultMap[E,0];
        let inseparable=rvECMultMap[E,7];
        assert rvECMapApply(inseparable,p) == rvECScale[7](p);
        rvECMapDegree(inseparable);
        rvECMapCertify(inseparable);
        "
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvECMapDegree = 49\nrvECMapCertify = degree=49, cleared identity=yes\n"
    );
}

#[test]
fn unramified_residue_coefficients_survive_exact_duality() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram W[5,1]{2,0,1};
            let z=rvUnramElem[W]{0,1};
            rvec E[W](z,rvUnramElem[W]{1});
            let k=rveckernel[E,3]{[1,1],[1]};
            rvecisogeny I[E,B](k);
            let f=rvECMap(rvECIso[I]);
            let d=rvECDual(f);
            assert rvECMapCompose(d,f) == rvECMultMap[E,3];
            assert rvECMapCompose(f,d) == rvECMultMap[B,3];
            assert rvECDual(d) == f;
            rvECMapCertify(d);
            "
        )),
        "rvECMapCertify = degree=3, cleared identity=yes\n"
    );
}

#[test]
fn endpoints_degrees_and_poles_are_checked() {
    for (extra, expect) in [
        ("show rvECMapCompose(f,f);", "mismatched endpoints"),
        ("show rvECMapApply(f,rvECInfinity[B]);", "rvECMapApply"),
        ("show rvECChart[E,B](one);", "does not give the target"),
        ("show rvECPush(iso,k);", "coprime prime"),
        ("show rvECDual(rvECId[E]);", "certified prime kernel"),
        (
            "show rvECMapXNum[0](rvECMultMap[E,0]);",
            "no affine rational coordinates",
        ),
        (
            "show rvECMultMap[E,18446744073709551616];",
            "degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F7}{extra}")));
        assert!(error.contains(expect), "{extra} -> {error}");
    }
}
