mod common;

use common::{diagnostic, example, interpret, stdout};

const F41: &str = include_str!("../examples/elliptic-cm-hom-ideal-composition.asl");

#[test]
fn actual_f41_dual_products_generate_the_full_principal_order() {
    let output = stdout(example("elliptic-cm-hom-ideal-composition"));
    for expected in [
        "rational_rank=2, product_degrees=[4, 6, 6, 9]",
        "coordinates=[([2, 0], 1), ([1, 1], 1), ([1, -1], 1), ([3, 0], 1)]",
        "common_denominator=1, hermite_numerators=[[1, 0], [0, 1]], smith_numerators=[1, 1]",
        "integral_index=1, target_intersection_index=1, product_intersection_index=1",
        "orientation_outer_inner=[-1, 1], full_ideal_product_certified=true",
        "rvCMHomCoordinates = numerator=[0, 1], denominator=1, in_integral_span=true, degree=5",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn thin_product_and_fractional_target_have_different_index_certificates() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvcmhom HT(f,rvCMHomMap[H,0,2]);
         rvcmhomcompose Thin[HD,HT,R];
         rvCMHomProductData[Thin];
         assert rvCMHomCoordMap(rvCMHomProductBasis[Thin,0])==rvECMultMap[E,2];
         let tw=rvCMHomMap[R,0,2];
         rvcmhomframe RT(id,tw);
         let wc=rvCMHomCoords[RT](w);
         assert rvCMHomCoordMap(wc)==w;
         rvCMHomCoordinates(wc);
         assert add(wc,neg(wc))==zero[RVCMHomCoordinate(RT)];
         assert add(wc,wc)==rvCMHomCoord[RT]{{0,1}};
         rvcmhomcompose Fractional[HD,H,RT];
         rvCMHomProductData[Fractional];
         assert rvCMHomCoordMap(rvCMHomProductBasis[Fractional,0])==id;
         assert rvCMHomCoordMap(rvCMHomProductBasis[Fractional,1])==w;"
    )));
    for expected in [
        "smith_numerators=[1, 2], integral_index=2, target_intersection_index=2, product_intersection_index=1",
        "rvCMHomCoordinates = numerator=[0, 1], denominator=2, in_integral_span=false, degree=5",
        "coordinates=[([2, 0], 1), ([2, 1], 2), ([2, -1], 2), ([3, 0], 1)]",
        "common_denominator=2, hermite_numerators=[[2, 0], [0, 1]], smith_numerators=[1, 2], integral_index=nonintegral, target_intersection_index=1, product_intersection_index=2",
    ] {
        assert!(output.contains(expected), "{output}");
    }
    assert_eq!(
        output.matches("full_ideal_product_certified=false").count(),
        2
    );
}

#[test]
fn separately_oriented_frames_agree_and_actual_composition_is_associative() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvcmhomframe Rev(id,rvECNegate(w));
         rvcmhomcompose Opposite[HD,H,Rev];
         rvCMHomProductData[Opposite];
         assert rvCMHomCoordMap(rvCMHomProductBasis[Opposite,1])==rvECNegate(w);
         rvcmhom HShear(f,rvCMHomMap[H,1,1]);
         rvcmhomcompose Shear[HD,HShear,R];
         rvCMHomProductData[Shear];
         assert rvCMHomProductColumn[Shear,1]==rvCMHomCoord[R]{{3,1}};
         assert rvCMHomProductColumn[Shear,3]==rvCMHomCoord[R]{{4,-1}};
         let idF=rvECMultMap[F,1];
         let wF=rvECSum(rvECMapCompose(g,fd),rvECNegate(idF));
         assert rvECMapCompose(wF,wF)==rvECMultMap[F,-5];
         rvcmhomframe RF(idF,wF);
         rvcmhomcompose ReverseProduct[H,HD,RF];
         rvCMHomProductData[ReverseProduct];
         let lhs=rvCMHomCoordCompose[H](rvCMHomCoordCompose[RF](fc,gdc),fc);
         let rhs=rvCMHomCoordCompose[H](fc,rvCMHomCoordCompose[R](gdc,fc));
         assert lhs==rhs;
         assert rvCMHomCoordMap(lhs)==rvECMapCompose(rvECMapCompose(f,gd),f);"
    )));
    assert!(output.contains("orientation_outer_inner=[1, -1], full_ideal_product_certified=true"));
    assert_eq!(
        output.matches("full_ideal_product_certified=true").count(),
        4
    );
}

#[test]
fn actual_inverse_has_reversed_named_endpoints_and_canonical_coordinates() {
    stdout(interpret(&format!(
        "{F41}
         rvec S[K](rvUnramElem[K]{{16}},rvUnramElem[K]{{22}});
         let scale=rvECChart[E,S](rvUnramElem[K]{{2}});
         let back=rvECChart[S,E](rvUnramElem[K]{{21}});
         rvcmhomframe Iso(scale,rvECMapCompose(scale,w));
         rvcmhomframe Back(back,rvECMapCompose(w,back));
         let iso=rvCMHomCoord[Iso]{{1,0}};
         let inverse=rvCMHomCoordInverse[Back](iso);
         assert rvCMHomCoordMap(inverse)==back;
         assert rvCMHomCoordCompose[R](inverse,iso)==rvCMHomCoord[R]{{1,0}};
         assert rvCMHomCoords[R](rvECMultMap[E,0])==rvCMHomCoord[R]{{0,0}};"
    )));
}

#[test]
fn exact_spanning_after_base_change_still_needs_the_endpoint_order_proof() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvunram L[41,1]{{38,0,1}};
         rvec EL[L](rvUnramElem[L]{{1}},rvUnramElem[L]{{17}});
         rvec FL[L](rvUnramElem[L]{{30}},rvUnramElem[L]{{3}});
         rvcmhom HL(rvECExtend[EL,FL](f),rvECExtend[EL,FL](g));
         rvcmhomdual DL[HL](rvECExtend[FL,EL](fd),rvECExtend[FL,EL](gd));
         rvcmhomframe RL(rvECMultMap[EL,1],rvECExtend[EL,EL](w));
         rvcmhomcompose ExactSpan[DL,HL,RL];
         rvCMHomProductData[ExactSpan];"
    )));
    let summary = output
        .lines()
        .find(|line| line.starts_with("rvCMHomProductData[ExactSpan]"))
        .expect("base-changed product report");
    assert!(summary.contains("integral_index=1"), "{summary}");
    assert!(
        summary.contains("full_ideal_product_certified=false"),
        "{summary}"
    );
}

#[test]
fn rejects_unverified_duals_bad_endpoints_and_nonunit_morphism_inversion() {
    for (query, expected) in [
        (
            "rvcmhomdual Bad[H](rvECNegate(fd),gd);",
            "exact Rosati sandwich",
        ),
        ("rvcmhomdual Bad[H](gd,fd);", "same degree"),
        ("rvcmhomdual Bad[H](f,g);", "CM Hom dual map"),
        ("rvcmhomcompose Bad[H,HD,R];", "composable named endpoints"),
        ("show rvCMHomCoords[R](f);", "coordinate map endpoints"),
        ("show rvCMHomCoordInverse[HD](fc);", "only a degree-one"),
        (
            "show rvCMHomCoordCompose[R](fc,fc);",
            "composable named endpoints",
        ),
        (
            "show rvCMHomProductColumn[P,4];",
            "column index is out of range",
        ),
        (
            "show rvCMHomProductBasis[P,2];",
            "basis index is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F41}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}
