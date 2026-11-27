mod common;

use common::{diagnostic, example, interpret, stdout};

const F13: &str = "rvunram F[13,1]{0,1};
    rvec G[F](rvUnramElem[F]{4},rvUnramElem[F]{0});
    rvec V[F](rvUnramElem[F]{2},rvUnramElem[F]{1});";

#[test]
fn normalized_two_velu_edges_and_local_schematic_signs() {
    let output = stdout(example("two-adic-volcano"));
    for expected in [
        "rvEllLocalOrder[G, 2] = r=1, bound=1, m=2, a=1, trace=2, norm=2, discriminant=-4",
        "rvEllLocalOrder[V, 2] = r=0, bound=1, m=1, a=0, trace=6, norm=13, discriminant=-16",
        "rvEllEdge[Flat] = horizontal, r_source=1, r_target=1",
        "rvEllEdge[DownA] = descending, r_source=1, r_target=0",
        "rvEllEdge[DownB] = descending, r_source=1, r_target=0",
        "rvEllEdge[Up] = ascending, r_source=0, r_target=1",
        "rvECMapDegree = 2",
        "rvEllLocalOrder[S4, 2] = r=2, bound=2, m=4, a=1, trace=0, norm=1, discriminant=-4",
        "rvEllLocalOrder[N4, 2] = r=1, bound=2, m=2, a=1, trace=0, norm=4, discriminant=-16",
        "rvEllNeighborhood[G, 2] = r_source=1, ascending=0, horizontal=1, descending=2",
        "rvEllNeighborhood[V, 2] = r_source=0, ascending=1, horizontal=0, descending=0",
        "(kernel=[[0], [1]], a=[10], b=[0], j=[12], horizontal, r=1)",
        "(kernel=[[10], [1]], a=[5], b=[12], j=[1], descending, r=0)",
        "(kernel=[[3], [1]], a=[5], b=[1], j=[1], descending, r=0)",
        "(kernel=[[11], [1]], a=[10], b=[0], j=[12], ascending, r=1)",
    ] {
        assert!(output.contains(expected), "{expected}\n{output}");
    }
}

#[test]
fn the_linear_kernel_and_both_dual_composites_use_the_true_maps() {
    let source = format!(
        "{F13}
         let k=rveckernel[G,2]{{[10],[1]}};
         rvecisogeny D[G,L](k);
         let i=rvECIso[D]; let f=rvECMap(i); let d=rvECDual(f);
         assert rvECGetKernel(i)==k;
         assert rvECXNum[0](i)==rvUnramElem[F]{{5}};
         assert rvECXNum[1](i)==rvUnramElem[F]{{10}};
         assert rvECXNum[2](i)==rvUnramElem[F]{{1}};
         assert rvECXDen[0](i)==rvUnramElem[F]{{10}};
         assert rvECXDen[1](i)==rvUnramElem[F]{{1}};
         assert rvECXDen[2](i)==rvUnramElem[F]{{0}};
         assert rvECDNum[0](i)==rvUnramElem[F]{{4}};
         assert rvECDNum[1](i)==rvUnramElem[F]{{7}};
         assert rvECDDen[0](i)==rvUnramElem[F]{{9}};
         assert rvECKCoeff[0](k)==rvUnramElem[F]{{10}};"
    );
    assert_eq!(
        stdout(interpret(&format!(
            "{source}
         assert rvECMapKCoeff[0](d)==rvUnramElem[F]{{6}};
         assert rvECMapCompose(d,f)==rvECMultMap[G,2];
         assert rvECMapCompose(f,d)==rvECMultMap[L,2];
         assert rvECDual(d)==f;
         let p=rvECPoint[G](rvUnramElem[F]{{2}},rvUnramElem[F]{{4}});
         assert rvECApply(i,p)==rvECMapApply(f,p);
         let kpoint=rvECPoint[G](rvUnramElem[F]{{3}},rvUnramElem[F]{{0}});
         assert rvECApply(i,kpoint)==rvECInfinity[L];
         assert rvECApply(i,rvECInfinity[G])==rvECInfinity[L];
         rvECCertify[D]; rvECMapCertify(d);"
        ))),
        concat!(
            "rvECCertify[D] = degree=2, cleared identity=yes\n",
            "rvECMapCertify = degree=2, cleared identity=yes\n",
        )
    );
}

#[test]
fn degree_two_duality_handles_nonprime_coefficients_and_coprime_push() {
    assert_eq!(
        stdout(interpret(
            "rvunram F[13,1]{11,0,1};
         rvec V[F](rvUnramElem[F]{2},rvUnramElem[F]{1});
         rvecisogeny I[V,L](rveckernel[V,2]{[1,11],[1]});
         let f=rvECMap(rvECIso[I]); let d=rvECDual(f);
         assert rvECMapKCoeff[0](d)==rvUnramElem[F]{11,4};
         assert rvECMapCompose(d,f)==rvECMultMap[V,2];
         assert rvECMapCompose(f,d)==rvECMultMap[L,2];
         assert rvECDual(d)==f;
         rvunram G[37,1]{0,1}; rvec E[G](rvUnramElem[G]{1},rvUnramElem[G]{0});
         let k2=rveckernel[E,2]{[0],[1]};
         let k3=rveckernel[E,3]{[34],[1]};
         rvecisogeny V2[E,P](k2); rvecisogeny V3[E,Q](k3);
         rvecisogeny U3[P,S](rvECPush(rvECIso[V2],k3));
         rvecisogeny U2[Q,T](rvECPush(rvECIso[V3],k2));
         let a=rvECMapCompose(rvECMap(rvECIso[U3]),rvECMap(rvECIso[V2]));
         let b=rvECMapCompose(rvECMap(rvECIso[U2]),rvECMap(rvECIso[V3]));
         assert a==rvECMapCompose(rvECChart[T,S](rvUnramElem[G]{1}),b);
         rvECMapDegree(a);",
        )),
        "rvECMapDegree = 6\n"
    );
}

#[test]
fn full_even_eliminant_separates_branch_and_nonbranch_coordinates() {
    let source = format!("{F13} rvEllFrobCert[V,2,1]; rvEllFrobCert[G,2,1];");
    let output = stdout(interpret(&source));
    for expected in [
        "rvEllFrobCert[V, 2, 1] = m=2, a=1, representative=1, deg(psi)=3, x=false, y=true",
        "rvEllFrobCert[G, 2, 1] = m=2, a=1, representative=1, deg(psi)=3, x=true, y=true",
    ] {
        assert!(output.contains(expected), "{output}");
    }
    assert_eq!(output.matches("frob_y=[]").count(), 2);
    let source = "rvunram W[17,1]{0,1}; rvec E[W](rvUnramElem[W]{1},rvUnramElem[W]{0});
        rvEllFrobCert[E,2,2,-1]; rvEllFrobCert[E,2,3,1];";
    let output = stdout(interpret(source));
    assert!(
        output.contains("m=4, a=3, representative=-1, deg(psi)=9, x=true, y=false"),
        "{output}"
    );
    assert!(
        output.contains("m=8, a=1, representative=1, deg(psi)=33"),
        "{output}"
    );
}

#[test]
fn the_default_uses_the_second_half_trace_candidate_when_required() {
    let output = stdout(interpret(
        "rvunram F[17,1]{0,1}; rvec B[F](rvUnramElem[F]{9},rvUnramElem[F]{0});
         assert rvEllFrobCoordinates[B,2,2]==(S.hi,S.hi);
         assert rvEllFrobCoordinates[B,2,2,1]==(S.hi,S.lo);
         rvEllLocalOrder[B,2];
         rvunram G[11,1]{0,1}; rvec Z[G](rvUnramElem[G]{1},rvUnramElem[G]{8});
         rvEllLocalOrder[Z,2];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvEllLocalOrder[B, 2] = r=2, bound=2, m=4, a=3, trace=-2, norm=2, discriminant=-4\n",
            "rvEllLocalOrder[Z, 2] = r=0, bound=0, m=1, a=0, trace=6, norm=11, discriminant=-8\n",
        )
    );
}

#[test]
fn inert_and_split_neighbors_keep_every_actual_edge() {
    let output = stdout(interpret(
        "rvunram F[37,1]{0,1}; rvec A[F](rvUnramElem[F]{1},rvUnramElem[F]{0}); rvEllNeighborhood[A,3];
         rvunram G[41,1]{0,1}; rvec B[G](rvUnramElem[G]{6},rvUnramElem[G]{0}); rvEllNeighborhood[B,5];",
    ));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(
        lines[0].contains("ascending=0, horizontal=0, descending=4"),
        "{output}"
    );
    for model in [
        "a=[17], b=[19], j=[15]",
        "a=[25], b=[28], j=[30]",
        "a=[25], b=[9], j=[30]",
        "a=[17], b=[18], j=[15]",
    ] {
        assert!(lines[0].contains(model), "{output}");
    }
    assert!(
        lines[1].contains("ascending=0, horizontal=2, descending=4"),
        "{output}"
    );
    for model in [
        "a=[15], b=[0], j=[6]",
        "a=[24], b=[0], j=[6]",
        "a=[35], b=[23], j=[19]",
        "a=[11], b=[21], j=[34]",
        "a=[35], b=[18], j=[19]",
        "a=[11], b=[20], j=[34]",
    ] {
        assert!(lines[1].contains(model), "{output}");
    }
    assert_eq!(lines[1].matches(", horizontal, r=1").count(), 2);
}

#[test]
fn rejects_wrong_two_kernels_and_keeps_the_other_domain_guards() {
    for (query, expected) in [
        (
            "show rveckernel[G,2]{[0],[0],[1]};",
            "must be monic and linear",
        ),
        (
            "show rveckernel[G,2]{[12],[1]};",
            "not a root of the Weierstrass cubic",
        ),
        (
            "rvEllFrobCert[G,2,2,2];",
            "scalar candidate must be prime to the torsion order",
        ),
        ("rvEllLocalOrder[G,4];", "must be 2 or an odd prime"),
        (
            "rvEllNeighborhood[G,13];",
            "equals the residue characteristic",
        ),
        (
            "rvunram K[3,1]{0,1}; rvec X[K](rvUnramElem[K]{1},rvUnramElem[K]{0});",
            "p>=5",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F13} {query}")));
        assert!(error.contains(expected), "{error}");
    }
    let empty = stdout(interpret(
        "rvunram F[13,1]{0,1}; rvec E[F](rvUnramElem[F]{0},rvUnramElem[F]{2}); rvECIsogenies[E,2];",
    ));
    assert!(empty.ends_with("[]\n"), "{empty}");
}
