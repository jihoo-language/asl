mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const CYCLIC: &str = "
    rvgroup H[4];rvgroup E[4];
    rvcrossed C[H,E]{0,2,0,2}{{0,1,2,3},{0,3,2,1},{0,1,2,3},{0,3,2,1}};
    rvpostnikov P[C,2];let k=rvPostnikovK[P];
";
const TWISTED: &str = "
    rvgroup H[3];rvgroup E[2];
    rvcrossed C[H,E]{0,0,0}{{0,1,2},{0,2,1}};
    rvpostnikov P[C,1];let k=rvPostnikovK[P];
";

#[test]
fn rp3_cycle_detects_the_nonzero_obstruction_and_the_split_lifts() {
    let actual = stdout(example("postnikov-flat-lift"));
    let lines: Vec<_> = actual.lines().collect();
    assert_eq!(lines.len(), 8);
    assert_eq!(
        lines[0],
        "rvFlatData[Cover] = prime=2, chains=[40, 232, 384, 192, 0], betti=[1, 1, 1, 1]"
    );
    assert!(lines[1].starts_with("rvFlatObstruction[Cover] = obstructed=true, pulled_support=1"));
    assert!(lines[1].ends_with("fixed_cardinality=0"));
    assert_eq!(lines[2], "rvFlatCoeff[dp,cpdm,bpcmdp,apbmcpdm] = 1");
    assert_eq!(lines[3], "rvFlatCyclePair = 1");
    assert_eq!(
        lines[4],
        "rvFlatObstruction[Cover] = obstructed=false, pulled_support=0, reduced=[], lift_classes=2, fixed_cardinality=2"
    );
    assert_eq!(lines[5], "rvFlatLiftCheck[Cover] = true");
    assert_eq!(
        lines[6],
        "rvFlatObstruction[Trivial] = obstructed=false, pulled_support=0, reduced=[], lift_classes=2, fixed_cardinality=2"
    );
    assert_eq!(lines[7], "rvFlatLiftCheck[Trivial] = true");
}

#[test]
fn actual_obstruction_prevents_a_lift_but_a_truncated_chain_is_not_a_cycle() {
    let source = include_str!("../examples/postnikov-flat-lift.asl");
    let err = diagnostic(interpret(&format!("{source}show rvFlatLift[Cover](k);")));
    assert!(
        err.contains("flat lift is obstructed by a nonzero pulled-back Postnikov class"),
        "{err}"
    );
    let err = diagnostic(interpret(&format!(
        "{source}rvFlatCyclePair(pulled){{1,{}}};",
        vec!["0"; 191].join(",")
    )));
    assert!(
        err.contains("supplied dual local chain is not a cycle"),
        "{err}"
    );
}

#[test]
fn twisted_pullbacks_commute_with_coboundary_and_gauge_translates_checked_lifts() {
    let source = format!(
        "{TWISTED}{}",
        r"
        space X{x0<x1<x2<x3<x4};
        rvflat F[X,P]{x0<x1:1,x1<x2:1,x2<x3:1,x3<x4:1,x0<x4:0};
        rvFlatEdge[F,x0,x3]; rvFlatEdge[F,x0,x4]; rvFlatEdge[F,x1,x1];
        rvFlatData[F];
        let f=rvFlatCochain[F,0]{0,1,0,0,0};
        let df=rvFlatCoboundary(f);
        rvFlatCoeff[x0,x1](df); rvFlatCoeff[x1,x2](df);
        assert rvFlatCoboundary(df)==zero[RVFlatCochain(F,2)];
        assert rvFlatCoboundary(rvFlatBoundaryWitness(df))==df;
        let b=rvPostnikovCochain[P,2]{1};
        let changed=rvPostnikovGauge(k,b);
        let local=rvFlatPullback[F](changed);
        let b_local=rvFlatPullback[F](b);
        assert local==rvFlatCoboundary(b_local);
        assert rvFlatCoboundary(local)==zero[RVFlatCochain(F,4)];
        rvFlatCoeff[x0,x1,x2,x3](local);
        rvFlatObstruction[F](changed);
        let solved=rvFlatLift[F](changed);
        let translated=add(rvFlatLift[F](k),b_local);
        rvFlatLiftCheck[F](changed,solved);
        rvFlatLiftCheck[F](changed,translated);
        let difference=add(translated,neg(solved));
        assert rvFlatCoboundary(rvFlatBoundaryWitness(difference))==difference;
        assert rvFlatLiftDifference[F](changed,solved,translated)==zero[RVFlatCochain(F,2)];
        assert rvFlatCoboundary(rvFlatBoundaryWitness(local))==local;
        "
    );
    let actual = stdout(interpret(&source));
    let lines: Vec<_> = actual.lines().collect();
    assert_eq!(lines[0], "rvFlatEdge[F,x0,x3] = 1");
    assert_eq!(lines[1], "rvFlatEdge[F,x0,x4] = 0");
    assert_eq!(lines[2], "rvFlatEdge[F,x1,x1] = 0");
    assert_eq!(
        lines[3],
        "rvFlatData[F] = prime=3, chains=[5, 10, 10, 5, 1], betti=[1, 0, 0, 0]"
    );
    assert_eq!(lines[4], "rvFlatCoeff[x0,x1] = 2");
    assert_eq!(lines[5], "rvFlatCoeff[x1,x2] = 2");
    assert_eq!(lines[6], "rvFlatCoeff[x0,x1,x2,x3] = 1");
    assert!(lines[7].starts_with("rvFlatObstruction[F] = obstructed=false, pulled_support="));
    assert!(lines[7].ends_with("reduced=[], lift_classes=1, fixed_cardinality=3"));
    assert_eq!(
        lines[8..],
        ["rvFlatLiftCheck[F] = true", "rvFlatLiftCheck[F] = true"]
    );
}

#[test]
fn inversion_holonomy_and_rational_cardinality_depend_on_the_frozen_background() {
    let source = format!(
        "{TWISTED}{}",
        r"
        space Graph{a,b,c,x,y,z,w,a<x,a<y,a<z,a<w,b<x,b<y,c<z,c<w};
        rvflat F[Graph,P]{a<x:1,a<y:0,a<z:0,a<w:0,b<x:0,b<y:0,c<z:0,c<w:0};
        rvflat T[Graph,P]{a<x:0,a<y:0,a<z:0,a<w:0,b<x:0,b<y:0,c<z:0,c<w:0};
        rvFlatData[F];rvFlatData[T];rvFlatObstruction[F](k);rvFlatObstruction[T](k);
        let edge=rvFlatCochain[F,1]{0,0,1,0,0,0,0,0};
        let vertex=rvFlatCochain[F,0]{0,0,1,0,0,0,0};
        rvFlatCyclePair(edge){0,0,1,2,0,0,2,1};
        rvFlatCyclePair(rvFlatCoboundary(vertex)){0,0,1,2,0,0,2,1};
        space Empty{};
        rvflat Void[Empty,P]{};
        rvFlatObstruction[Void](k);
        assert (\x : RVFlatCochain(Void,2) => add(x,neg(x)))==
            (\x : RVFlatCochain(Void,2) => zero[RVFlatCochain(Void,2)]);
        "
    );
    let actual = stdout(interpret(&source));
    let lines: Vec<_> = actual.lines().collect();
    assert_eq!(
        lines[0],
        "rvFlatData[F] = prime=3, chains=[7, 8, 0, 0, 0], betti=[0, 1, 0, 0]"
    );
    assert_eq!(
        lines[1],
        "rvFlatData[T] = prime=3, chains=[7, 8, 0, 0, 0], betti=[1, 2, 0, 0]"
    );
    assert!(lines[2].ends_with("fixed_cardinality=1/3"));
    assert!(lines[3].ends_with("fixed_cardinality=1/3"));
    assert_eq!(lines[4..6], ["rvFlatCyclePair = 1", "rvFlatCyclePair = 0"]);
    assert!(lines[6].ends_with("fixed_cardinality=1"));

    let points = (0..100)
        .map(|i| format!("x{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let large = stdout(interpret(&format!(
        "{TWISTED}space X{{{points}}};rvflat F[X,P]{{}};rvFlatObstruction[F](k);"
    )));
    assert!(large
        .trim()
        .ends_with("lift_classes=1, fixed_cardinality=3^100"));
}

#[test]
fn covers_use_ordered_nonabelian_products() {
    let mul = |x: usize, y: usize| {
        let (a, b, c, d) = (x % 3, x / 3, y % 3, y / 3);
        (a + if b == 0 { c } else { 3 - c }) % 3 + 3 * ((b + d) % 2)
    };
    let mut source = String::from("rvgroup E{");
    for x in 0..6 {
        write!(
            source,
            "{{{}}},",
            (0..6)
                .map(|y| mul(x, y).to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
        .unwrap();
    }
    source.push_str(
        "};
        rvgroup H[2];rvcrossed C[H,E]{0,0}{{0,1},{0,1},{0,1},{0,1},{0,1},{0,1}};
        rvpostnikov P[C,1];
        space X{a<b<d,a<c<d};
        rvflat F[X,P]{a<b:3,b<d:1,a<c:1,c<d:4};
        rvFlatEdge[F,a,d];",
    );
    assert_eq!(stdout(interpret(&source)), "rvFlatEdge[F,a,d] = 5\n");
    let inconsistent = source.replace("c<d:4", "c<d:3");
    assert!(diagnostic(interpret(&inconsistent)).contains("flat labels disagree along"));
}

#[test]
fn flatness_and_lift_claims_are_checked_and_bounded() {
    let diamond = "space X{a<b<d,a<c<d};";
    let chain = "space X{a<b<c<d};rvflat F[X,P]{a<b:1,b<c:1,c<d:1};";
    let cases = [
        (format!("{diamond}rvflat F[X,P]{{a<b:1,b<d:0,a<c:0,c<d:0}};"), "flat labels disagree"),
        (format!("{diamond}rvflat F[X,P]{{a<b:1,b<d:1,c<d:0}};"), "missing flat label on cover"),
        (format!("{diamond}rvflat F[X,P]{{a<b:1,a<b:1}};"), "duplicate flat edge label"),
        (format!("{diamond}rvflat F[X,P]{{a<b:2}};"), "skeletal quotient object index is out of range"),
        (format!("{diamond}rvflat F[X,P]{{d<a:0}};"), "flat edge requires x < y"),
        (format!("{diamond}rvflat F[X,P]{{a<b:0,b<d:0,a<c:0,c<d:0,a<d:1}};"), "flat labels disagree"),
        (format!("{chain}rvFlatLiftCheck[F](k,zero[RVFlatCochain(F,2)]);"), "supplied lift does not have"),
        (format!("{chain}rvflat G[X,P]{{a<b:1,b<c:1,c<d:1}};rvFlatLiftCheck[F](k,zero[RVFlatCochain(G,2)]);"), "expected type RVFlatCochain(F,2)"),
        (format!("{chain}rvpostnikov Q[C,2];show rvFlatPullback[F](rvPostnikovK[Q]);"), "expected type RVPostnikovCochain(P,3)"),
        (format!("{chain}rvFlatCoeff[b,a](zero[RVFlatCochain(F,1)]);"), "coefficient needs a strict chain"),
        (format!("{chain}show rvFlatH2Basis[F,0];"), "local H^2 basis index"),
        (format!("{chain}show rvFlatCochain[F,2]{{0,0,0,2}};"), "degree-2 requires 4 residues"),
    ];
    for (source, expected) in cases {
        let error = diagnostic(interpret(&format!("{CYCLIC}{source}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
    let vertices = (0..28).map(|i| format!("x{i}")).collect::<Vec<_>>();
    let edges = vertices
        .windows(2)
        .map(|w| format!("{}<{}:0", w[0], w[1]))
        .collect::<Vec<_>>()
        .join(",");
    let error = diagnostic(interpret(&format!(
        "{CYCLIC}space X{{{}}};rvflat F[X,P]{{{edges}}};",
        vertices.join("<")
    )));
    assert!(
        error.contains("order complex exceeds the simplex limit"),
        "{error}"
    );
}
