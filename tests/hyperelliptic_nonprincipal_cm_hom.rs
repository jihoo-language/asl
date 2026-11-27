mod common;

use common::{diagnostic, example, interpret, stdout};

const F41: &str = include_str!("../examples/hyperelliptic-nonprincipal-cm-hom.asl");

#[test]
fn exact_nonprincipal_norm_form_and_actual_geometric_graphs() {
    let out = stdout(example("hyperelliptic-nonprincipal-cm-hom"));
    for expected in [
        "j=[[40], [12]], trace=12, frobenius_conductor=1, fundamental_discriminant=-20, maximal_order_certified=true",
        "degrees=[2, 3], degree_sum=7, degree_difference=3, norm_form=[2, 2, 3], primitive_reduced=[2, 2, 3]",
        "span_minimum=2, full_minimum=2",
        "graph_gram=[[0, 1, 2, 3], [1, 0, 1, 1], [2, 1, 0, 3], [3, 1, 3, 0]], graph_discriminant=-20",
        "image_basis_1_omega=[[2, 0], [1, 1]], span_index=1, full_lattice=true, nonprincipal=true",
        "rvECMapDegree = 7",
        "rvECMapDegree = 3",
        "polarization_units=[1, 1], map_columns=[[[1, 0], [0, 2]], [[0, 0], [0, 2]]]",
        "cycle_lifts=[[0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 1, 2, 1, 2, 0], [0, 1, 2, 0, 0, 0]]",
        "twisted_frobenius=[[2, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 2]]",
        "abelian_tate_rank=4, locally_saturated=true",
        "rvCMHomCyclePair = product=3, normalized_kummer=3/2",
        "rvCMHomWedgePair = residue=1, modulus=3",
    ] {
        assert!(out.contains(expected), "{out}");
    }
}

#[test]
fn graph_generators_change_and_proper_index_can_still_be_a_local_unit() {
    let out = stdout(interpret(&format!(
        "{F41}
         rvcmhom Minus(phi,rvECNegate(psi));
         rvCMHomData[Minus];
         rvcmhom Swap(psi,phi);
         rvCMHomData[Swap];
         rvcmhom Thin(phi,rvCMHomMap[H,0,2]);
         rvCMHomData[Thin];
         rvcmhommark MinusM[Minus,C,B,BE,BF];
         rvcmhommark SwapM[Swap,C,B,BE,BF];
         rvcmhommark ThinM[Thin,C,B,BE,BF];
         let swapped=rvCMHomLift[SwapM](rvCMHomCycle[Swap]{{0,0,0,1}});
         assert swapped==rvCMHomChange[SwapM](g2);
         let negative=rvCMHomLift[MinusM](rvCMHomCycle[Minus]{{0,0,0,1}});
         let twiceaxes=rvCMHomLift[M](rvCMHomCycle[H]{{2,6,0,0}});
         assert add(rvCMHomChange[M](negative),g3)==twiceaxes;
         let x=add(trans,g2);
         assert rvCMHomChange[ThinM](rvCMHomAlgebraic(x))==rvCMHomAlgebraic(rvCMHomChange[ThinM](x));
         assert add(graph3,neg(graph3))==zero[RVCMHomCycle(H)];"
    )));
    assert!(out.contains("degree_sum=3, degree_difference=7, norm_form=[2, -2, 3]"));
    assert!(out.contains("norm_form=[3, 2, 2], primitive_reduced=[2, 2, 3]"));
    assert!(out.contains("norm_form=[2, 4, 12], primitive_reduced=[1, 0, 5], content=2"));
    assert!(out.contains("graph_discriminant=-80, image_basis_1_omega=[[2, 0], [2, 2]], span_index=2, full_lattice=false"));
}

#[test]
fn changing_the_weil_root_intertwines_projectors_and_base_frobenius() {
    let out = stdout(interpret(&format!(
        r"{F41}
         rvjactorsionbasis BT[SL,3](rvUnramElem[L]{{20,16}}){{e1,e2,neg(add(f1,e1)),f2}};
         rvcmhommark MT[H,C,BT,BE,BF];
         let gt=rvCMHomLift[MT](graph3);
         assert gt==rvCMHomChange[MT](g3);
         assert rvCMHomChange[M](gt)==g3;
         assert (\x:RVCMHomWedge(M)=>rvCMHomChange[MT](rvCMHomAlgebraic(x)))
            ==(\x:RVCMHomWedge(M)=>rvCMHomAlgebraic(rvCMHomChange[MT](x)));
         assert (\x:RVCMHomWedge(M)=>rvCMHomChange[MT](rvCMHomFrob[1](x)))
            ==(\x:RVCMHomWedge(M)=>rvCMHomFrob[1](rvCMHomChange[MT](x)));
         rvCMHomMarkData[MT];"
    )));
    assert!(out.contains("rvCMHomMarkData[MT] = level=3, polarization_units=[2, 2]"));
}

#[test]
fn nonfundamental_frobenius_alone_does_not_certify_the_maximal_orders() {
    let setup = format!(
        "{F41}
         rvec EL[L](rvUnramElem[L]{{1}},rvUnramElem[L]{{17}});
         rvec FL[L](rvUnramElem[L]{{30}},rvUnramElem[L]{{3}});
         rvcmhom Ext(rvECExtend[EL,FL](phi),rvECExtend[EL,FL](psi));
         rvCMHomData[Ext];
         rvecendtorsion EO[E];
         rvecendtorsion FO[F];
         rvcmhom Certified[EO,FO](phi,psi);
         assert rvCMHomClass[O,Certified]==rvCMHomClass[O,H];"
    );
    let out = stdout(interpret(&setup));
    let data = out
        .lines()
        .find(|line| line.starts_with("rvCMHomData[Ext]"))
        .unwrap();
    assert!(data.contains("trace=62, frobenius_conductor=12, fundamental_discriminant=-20, maximal_order_certified=false"), "{data}");
    assert!(
        data.contains("span_index=unproved, full_lattice=unproved, nonprincipal=unproved"),
        "{data}"
    );
    let error = diagnostic(interpret(&format!("{setup}\nshow rvCMHomClass[O,Ext];")));
    assert!(
        error.contains("requires a full certified maximal-order span"),
        "{error}"
    );
    let error = diagnostic(interpret(&format!(
        "{setup}\nrvcmhom Wrong[FO,EO](phi,psi);"
    )));
    assert!(error.contains("actual map endpoints"), "{error}");
}

#[test]
fn rejects_dependent_maps_and_wrong_geometry_without_fabricating_projectors() {
    for (request, expected) in [
        ("rvcmhom Bad(phi,rawpsi);", "same named source and target"),
        ("rvcmhom Bad(phi,rvCMHomMap[H,0,0]);", "positive degree"),
        ("rvcmhom Bad(phi,rvCMHomMap[H,2,0]);", "must be independent"),
        ("rvcmhom Bad(rvECMultMap[E,1],rvECMultMap[E,2]);", "distinct elliptic j-invariants"),
        ("rvcmhommark Bad[H,C,B,BF,BE];", "actual product factors"),
        ("rvcmhommark Bad[H,C,BE,BE,BF];", "same odd prime-power level"),
        ("rvcmorder U[-4]; show rvCMHomClass[U,H];", "matching order"),
        ("rvcmhom Lower(phi,rvCMHomMap[H,0,3]); rvcmhommark LM[Lower,C,B,BE,BF]; show rvCMHomAlgebraic(zero[RVCMHomWedge(LM)]);", "Gram is not a local unit"),
    ] {
        let error = diagnostic(interpret(&format!("{F41}\n{request}")));
        assert!(error.contains(expected), "{error}");
    }
}
