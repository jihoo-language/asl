mod common;

use common::{diagnostic, example, interpret, stdout};

const F41: &str = include_str!("../examples/elliptic-cm-fourier-mukai.asl");

#[test]
fn exact_fourier_mukai_pull_push_and_tensor_on_the_witnessed_product() {
    let output = stdout(example("elliptic-cm-fourier-mukai"));
    for expected in [
        "radicand=-5, adjoint_weight=2, matrix=[[[1, 0, 1], [2, 0, 1]], [[1, 1, 2], [1, 0, 1]]], determinant=[0, -1, 1], degree=5",
        "rvCMAbelianDet = 6",
        "rank=2, matrix=[[[3, 0, 1], [1, -1, 1]], [[1, 1, 2], [3, 0, 1]]], last=-1, square=16",
        "rank=2, matrix=[[[18, 0, 1], [16, -4, 1]], [[8, 2, 1], [11, 0, 1]]], last=-5, square=80",
        "rank=6, matrix=[[[-3, 0, 1], [1, -1, 1]], [[1, 1, 2], [-3, 0, 1]]], last=1, square=0",
    ] {
        assert!(output.contains(expected), "{output}");
    }
    let pairing: Vec<_> = output
        .lines()
        .filter_map(|line| line.strip_prefix("rvCMAbelianPair = "))
        .collect();
    assert_eq!(pairing, ["9", "9", "9", "67", "67", "45"]);
}

#[test]
fn changing_actual_anchor_and_endomorphism_bases_preserves_the_correspondence() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvcmhom HS(g,f);
         rvcmhomdual DS[HS](gd,fd);
         rvcmhomframe ES(rvECSum(id,w),w);
         rvcmhomframe FS(rvECNegate(wF),idF);
         rvcmabelian B[HS,DS,ES,FS];
         let d=rvCMAbelianDiv[B,3,3](rvCMHomCoords[HS](g));
         let x=rvCMAbelianMukai[2,-1](d);
         let m=rvCMAbelianMap[B](id,fd,g,idF);
         let md=rvCMAbelianDagger(m);
         assert rvCMAbelianBlock[1](md)==gd;
         assert rvCMAbelianBlock[2](md)==f;
         let detm=rvCMAbelianCompose(m,rvCMAbelianAdjugate(m));
         assert rvCMAbelianBlock[0](detm)==rvECNegate(w);
         assert rvCMAbelianBlock[3](detm)==rvECNegate(wF);
         rvCMAbelianMapData(m);
         rvCMAbelianData(rvCMAbelianPull(m,x));
         assert rvCMAbelianGraph(d)==rvCMHomCycle[HS]{{4,6,-1,0}};"
    )));
    assert!(output.contains("adjoint_weight=3, matrix=[[[1, 0, 1], [1, -1, 1]], [[1, 0, 1], [1, 0, 1]]], determinant=[0, 1, 1], degree=5"));
    assert!(output.contains(
        "rank=2, matrix=[[[18, 0, 1], [18, -6, 1]], [[6, 2, 1], [11, 0, 1]]], last=-5, square=80"
    ));
}

#[test]
fn graph_intersections_match_the_actual_marked_kummer_wedges() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvunram L[41,1]{{38,0,1}};
         rvhyperelliptic E[K]{{rvUnramElem[K]{{17}},rvUnramElem[K]{{1}},rvUnramElem[K]{{0}},rvUnramElem[K]{{1}}}};
         rvhyperelliptic F[K]{{rvUnramElem[K]{{3}},rvUnramElem[K]{{30}},rvUnramElem[K]{{0}},rvUnramElem[K]{{1}}}};
         rvellglue C[E,F](rvUnramElem[K]{{28}},rvUnramElem[K]{{2}},rvUnramElem[K]{{1}},rvUnramElem[K]{{40}});
         rvjacbase EL[E,L];
         rvjacbase FL[F,L];
         rvellglue CL[EL,FL](rvUnramElem[L]{{28}},rvUnramElem[L]{{2}},rvUnramElem[L]{{1}},rvUnramElem[L]{{40}});
         rvjacchart SL[CL](rvUnramElem[L]{{0,4}});
         let P=rvJacPoint[EL](rvUnramElem[L]{{34}},rvUnramElem[L]{{6}});
         let Q=rvJacPoint[EL](rvUnramElem[L]{{23}},rvUnramElem[L]{{0,16}});
         let Rr=rvJacPoint[FL](rvUnramElem[L]{{7}},rvUnramElem[L]{{8}});
         let T=neg(rvJacPoint[FL](rvUnramElem[L]{{31}},rvUnramElem[L]{{0,13}}));
         let e1=rvEllGlueMap[SL](P,zero[RVJac(FL)]);
         let e2=rvEllGlueMap[SL](zero[RVJac(EL)],Rr);
         let f1=rvEllGlueMap[SL](Q,zero[RVJac(FL)]);
         let f2=rvEllGlueMap[SL](zero[RVJac(EL)],T);
         rvjactorsionbasis B[SL,3](rvUnramElem[L]{{20,25}}){{e1,e2,f1,f2}};
         rvjactorsionbasis BE[EL,3](rvJacWeil[3](P,Q)){{P,Q}};
         rvjactorsionbasis BF[FL,3](rvJacWeil[3](Rr,T)){{Rr,T}};
         rvcmhommark G[H,C,B,BE,BF];
         rvCMAbelianPair(D,J);
         rvCMHomCyclePair(rvCMAbelianGraph(D),rvCMAbelianGraph(J));
         let graph_d=rvCMHomLift[G](rvCMAbelianGraph(D));
         let graph_j=rvCMHomLift[G](rvCMAbelianGraph(J));
         rvCMHomWedgePair(graph_d,graph_j);
         assert rvCMHomFrob[1](graph_d)==graph_d;"
    )));
    assert!(output.contains("rvCMAbelianPair = 16"));
    assert!(output.contains("rvCMHomCyclePair = product=16, normalized_kummer=8"));
    assert!(output.contains("rvCMHomWedgePair = residue=2, modulus=3"));
}

#[test]
fn rational_scaling_cancelled_wide_products_and_singular_pullback() {
    let output = stdout(interpret(&format!(
        "{F41}
         let thirds=rvCMAbelianScale[1,3](v);
         rvCMAbelianPair(thirds,thirds);
         assert rvCMAbelianFourier(thirds)==rvCMAbelianScale[1,3](rvCMAbelianFourier(v));
         assert rvCMAbelianPull(M,thirds)==rvCMAbelianScale[1,3](rvCMAbelianPull(M,v));
         assert add(thirds,neg(thirds))==zero[RVCMAbelianMukai(A)];
         let bound=rvCMAbelianScale[1267650600228229401496703205376,1](rvCMAbelianFromGraph[A](rvCMHomCycle[H]{{0,0,1,0}}));
         rvCMAbelianDet(bound);
         rvCMAbelianPair(bound,bound);
         let t=rvCMAbelianScale[1267650600228229401496703205376,1](rvCMAbelianLine(I));
         rvCMAbelianPair(t,t);
         let proj=rvCMAbelianMap[A](id,rvCMHomMap[HD,0,0],rvCMHomMap[H,0,0],rvECMultMap[F,0]);
         let pull=rvCMAbelianPull(proj,v);
         assert pull==rvCMAbelianMukai[2,0](rvCMAbelianDiv[A,3,0](zero[RVCMHomCoordinate(H)]));
         rvCMAbelianMapData(proj);"
    )));
    assert!(output.contains("rvCMAbelianPair = 16/9"));
    assert!(output.contains("rvCMAbelianDet = 0"));
    assert_eq!(output.matches("rvCMAbelianPair = 0\n").count(), 2);
    assert!(output.contains("determinant=[0, 0, 1], degree=0"));
}

#[test]
fn odd_fundamental_cm_discriminant_and_geometric_unit_blocks() {
    let output = stdout(interpret(
        r"
        rvunram k[7,1]{0,1};
        rvec E[k](rvUnramElem[k]{0},rvUnramElem[k]{3});
        rvec F[k](rvUnramElem[k]{0},rvUnramElem[k]{3});
        let f=rvECChart[E,F](rvUnramElem[k]{1});
        let g=rvECChart[E,F](rvUnramElem[k]{3});
        let fd=rvECChart[F,E](rvUnramElem[k]{1});
        let gd=rvECChart[F,E](rvUnramElem[k]{5});
        rvcmhomframe H(f,g);
        rvcmhomdual HD[H](fd,gd);
        let e0=rvECMultMap[E,1];
        let f0=rvECMultMap[F,1];
        rvcmhomframe RE(e0,rvECMapCompose(fd,g));
        rvcmhomframe RF(f0,rvECMapCompose(g,fd));
        rvcmabelian A[H,HD,RE,RF];
        let m=rvCMAbelianMap[A](e0,fd,g,f0);
        let h=rvCMAbelianDiv[A,2,2](rvCMHomCoord[H]{0,1});
        let v=rvCMAbelianMukai[1,1](h);
        assert rvCMAbelianPush(m,rvCMAbelianPull(m,v))==v;
        assert rvCMAbelianBlock[1](rvCMAbelianDagger(m))==gd;
        rvCMAbelianDet(h);
        rvCMAbelianMapData(m);",
    ));
    assert!(output.contains("rvCMAbelianDet = 3"));
    assert!(output.contains("radicand=-3, adjoint_weight=1"));
    assert!(output.contains("determinant=[1, -1, 2], degree=1"));
}

#[test]
fn rejects_singular_push_wrong_endpoints_and_unproved_markings() {
    for (query, expected) in [
        ("rvcmabelian Bad[H,HD,R,R];", "two endomorphism frames on the named factors"),
        ("rvcmhom Thin(f,rvCMHomMap[H,0,2]); rvcmabelian Bad[Thin,HD,R,RF];", "four full certified maximal-order"),
        ("rvcmhomframe Dbad(fd,rvECNegate(gd)); rvcmabelian Bad[H,Dbad,R,RF];", "Rosati sandwich"),
        ("show rvCMAbelianMap[A](id,g,fd,idF);", "product map block endpoints"),
        ("show rvCMAbelianDiv[A,1,1](rvCMHomCoord[R]{1,0});", "off-diagonal Hom coordinate"),
        ("show rvCMAbelianBlock[4](M);", "between zero and three"),
        ("show rvCMAbelianGraph(rvCMAbelianScale[1,2](D));", "outside the certified integral lattice"),
        ("show rvCMAbelianScale[1,0](v);", "denominator must be positive"),
        ("let m=rvCMAbelianMap[A](id,rvCMHomMap[HD,0,0],rvCMHomMap[H,0,0],rvECMultMap[F,0]); show rvCMAbelianPush(m,v);", "requires a nonsingular product isogeny"),
        ("rvcmabelian X[H,HD,R,RF]; let mx=rvCMAbelianMap[X](id,fd,g,idF); show rvCMAbelianPull(mx,v);", "CM abelian second operand"),
    ] {
        let error = diagnostic(interpret(&format!("{F41}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}
