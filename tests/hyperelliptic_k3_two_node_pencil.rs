mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F5: &str = r"
    rvunram F[5,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let three=rvUnramElem[F]{3};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    let r0=rvJacPoint[C](z,z);
    let r1=rvJacPoint[C](o,z);
    let r2=rvJacPoint[C](two,z);
    let r3=rvJacPoint[C](three,z);
    rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
";

#[test]
fn all_eight_sections_invert_the_pencil_on_their_entire_rational_tropes() {
    let mut program = format!(
        r"{F5}
        rvk3pencil P[B,0,1];
        rvk3pencil Reverse[B,1,0];
        let f=rvK3PencilDiv[P];
        assert f==rvK3PencilDiv[Reverse];
        assert f==rvKummerDiv[B]{{2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
        let central0=add(rvKummerDivScale[2](rvKummerTropeDiv[B,0]),add(add(rvKummerNodeDiv[B,4],rvKummerNodeDiv[B,7]),add(rvKummerNodeDiv[B,13],rvKummerNodeDiv[B,15])));
        let central1=add(rvKummerDivScale[2](rvKummerTropeDiv[B,1]),add(add(rvKummerNodeDiv[B,5],rvKummerNodeDiv[B,6]),add(rvKummerNodeDiv[B,12],rvKummerNodeDiv[B,14])));
        assert f==central0;
        assert f==central1;
        assert rvK3PencilParam[P](two,neg(o))==rvK3PencilParam[P](o,two);
        rvKummerDivPair(f,f);
        "
    );
    for mask in [4, 5, 6, 7, 12, 13, 14, 15] {
        writeln!(
            program,
            r"
            let section{mask}=\p:RVK3PencilBase(P)=>rvK3PencilSection[P,{mask}](p);
            let affine{mask}=\s:RVUnram(F)=>rvKummerTropeParam[B,{mask}](o,s);
            let infinity{mask}=rvKummerTropeParam[B,{mask}](z,o);
            assert (\p:RVK3PencilBase(P)=>rvK3PencilMap[P](section{mask}(p)))==(\p:RVK3PencilBase(P)=>p);
            assert (\s:RVUnram(F)=>section{mask}(rvK3PencilMap[P](affine{mask}(s))))==affine{mask};
            assert section{mask}(rvK3PencilMap[P](infinity{mask}))==infinity{mask};
            assert (\s:RVUnram(F)=>section{mask}(rvK3PencilParam[P](o,s)))==(\s:RVUnram(F)=>rvK3PencilSection[Reverse,{mask}](rvK3PencilParam[Reverse](o,s)));
            assert section{mask}(rvK3PencilParam[P](z,o))==rvK3PencilSection[Reverse,{mask}](rvK3PencilParam[Reverse](z,o));
            rvKummerDivPair(f,rvKummerTropeDiv[B,{mask}]);
            "
        )
        .unwrap();
    }
    assert_eq!(
        stdout(interpret(&program)),
        format!("rvKummerDivPair = 0\n{}", "rvKummerDivPair = 1\n".repeat(8))
    );
}

#[test]
fn entire_f5_map_matches_the_plain_secant_forms_and_its_census() {
    let output = stdout(interpret(&format!(
        r"{F5}
        rvk3pencil P[B,0,1];
        show (\x:RVKummerK3(C)=>rvK3PencilMap[P](x));
        rvK3PencilCensus[P];
        rvKummerK3Points[C,1];
        rvK3PencilFiber[P](rvK3PencilParam[P](z,o));
        rvK3PencilFiber[P](rvK3PencilParam[P](o,z));
        rvK3PencilFiber[P](rvK3PencilParam[P](o,o));
        rvK3PencilFiber[P](rvK3PencilParam[P](o,two));
        rvK3PencilFiber[P](rvK3PencilParam[P](o,three));
        rvK3PencilFiber[P](rvK3PencilParam[P](o,neg(o)));
        "
    )));
    let mut lines = output.lines();
    let graph = lines.next().unwrap();
    let mut counts = BTreeMap::new();
    let mut exceptional_on_secant = 0;
    for (point, _) in graph.match_indices("RVKummerK3(C){") {
        let point = &graph[point..];
        let end = point.find('}').unwrap();
        let coordinates = digits(&point[point.find('{').unwrap() + 1..end]);
        let base = &coordinates[..4];
        let relevant = if base == [0, 0, 0, 1] || base == [0, 1, 0, 0] {
            exceptional_on_secant += 1;
            &coordinates[4..8]
        } else {
            base
        };
        // For this secant the two reduced annihilator rows pick w and p.
        let [u, v] = [relevant[0], relevant[2]];
        let expected = if u == 0 {
            assert_ne!(v, 0);
            [0, 1]
        } else {
            let inverse = (1..5).find(|a| a * u % 5 == 1).unwrap();
            [1, inverse * v % 5]
        };
        let parameter = point[end..].split_once("RVK3PencilBase(P)[[").unwrap().1;
        let actual = digits(parameter.split_once("]]").unwrap().0);
        assert_eq!(actual, expected);
        *counts.entry(expected).or_insert(0) += 1;
    }
    assert_eq!(exceptional_on_secant, 12);
    assert_eq!(counts.len(), 6);
    assert_eq!(counts.values().sum::<usize>(), 96);
    let mut histogram = BTreeMap::new();
    for &count in counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }
    assert_eq!(
        lines.next().unwrap(),
        format!(
            "rvK3PencilCensus[P] = base points=6, total=96, fiber sizes={:?}",
            histogram.into_iter().collect::<Vec<_>>()
        )
    );
    assert_eq!(
        lines.next().unwrap(),
        "rvKummerK3Points[C, 1] = singular=16, nodes=16, resolved=96"
    );
    for (_, count) in counts {
        assert_eq!(
            lines.next().unwrap(),
            format!("rvK3PencilFiber[P] = {count}")
        );
    }
    assert!(lines.next().is_none());
}

#[test]
fn another_secant_uses_the_reduced_annihilator_of_two_finite_nodes() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            rvk3pencil P[B,4,7];
            let node=rvKummerTropeParam[B,0](o,o);
            assert rvK3PencilMap[P](node)==rvK3PencilParam[P](z,o);
            let other=rvKummerK3Translate[B,4](rvKummerTropeParam[B,0](o,z));
            assert rvKummerBlowdown(other)==rvKummer[C](o,o,z,neg(o));
            assert rvK3PencilMap[P](other)==rvK3PencilParam[P](o,three);
            let x=rvK3PencilSection[P,4](rvK3PencilParam[P](o,three));
            assert rvKummerOnTrope[B,4](x)==S.hi;
            assert rvK3PencilMap[P](x)==rvK3PencilParam[P](o,three);
            "
        ))),
        ""
    );
}

#[test]
fn nonsplit_jacobian_lifts_are_counted_on_the_resolved_surface() {
    let output = stdout(interpret(
        r"
        rvunram U[3,1]{1,0,1};
        let z=rvUnramElem[U]{};
        let o=rvUnramElem[U]{1};
        let i=rvUnramElem[U]{0,1};
        rvhyperelliptic C[U]{z,neg(o),z,z,z,o};
        let e1=rvJacPoint[C](z,z);
        let f1=rvJacPoint[C](o,z);
        let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
        let f2=add(add(e1,f1),rvJacPoint[C](i,z));
        rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};
        rvk3pencil P[B,0,1];
        let t=add(o,i);
        let twist=rvKummer[C](z,o,t,rvUnramMul(t,t));
        assert rvKummerLiftable(twist)==S.lo;
        let p=rvK3PencilParam[P](z,o);
        assert rvK3PencilMap[P](rvKummerSmooth(twist))==p;
        rvK3PencilCensus[P];
        rvKummerK3Points[C,1];
    ",
    ));
    assert!(output.contains("rvK3PencilCensus[P] = base points=10, total=248,"));
    assert!(output.contains("rvKummerK3Points[C, 1] = singular=104, nodes=16, resolved=248"));
}

#[test]
fn direct_section_inversion_works_without_enumerating_the_large_surface() {
    let program = format!(
        r"{F5}
        rvunram U[5,1]{{1,4,0,0,0,1}};
        rvjacbase L[C,U];
        let zu=rvUnramElem[U]{{}};
        let ou=rvUnramElem[U]{{1}};
        let twoU=rvUnramElem[U]{{2}};
        let threeU=rvUnramElem[U]{{3}};
        let radical=rvUnramElem[U]{{0,1}};
        let x0=rvJacPoint[L](zu,zu);
        let x1=rvJacPoint[L](ou,zu);
        let x2=rvJacPoint[L](twoU,zu);
        let x3=rvJacPoint[L](threeU,zu);
        rvjactorsionbasis BU[L,2](neg(ou)){{x0,add(add(x0,x1),x2),x1,add(add(x0,x1),x3)}};
        rvk3pencil P[BU,0,1];
        let p=rvK3PencilParam[P](ou,radical);
        let x=rvK3PencilSection[P,4](p);
        assert rvK3PencilMap[P](x)==p;
        assert rvKummerOnTrope[BU,4](x)==S.hi;
        show p;
        "
    );
    assert_eq!(
        stdout(interpret(&program)),
        "RVK3PencilBase(P)[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]]\n"
    );
    let error = diagnostic(interpret(&format!("{program} rvK3PencilCensus[P];")));
    assert!(error.contains("computation limit"), "{error}");
}

#[test]
fn rejects_degenerate_pairs_unusable_tropes_and_nominal_mismatches() {
    for (program, expected) in [
        ("rvk3pencil Q[B,1,1];", "requires two distinct nodes"),
        (
            "rvk3pencil Q[B,0,16];",
            "node mask must be between 0 and 15",
        ),
        (
            "rvk3pencil Q[Absent,0,1];",
            "unknown Jacobian torsion basis",
        ),
        ("rvhyperelliptic E[F]{z,neg(o),z,o}; rvjactorsionbasis Low[E,2](neg(o)){rvJacPoint[E](z,z),rvJacPoint[E](o,z)}; rvk3pencil Q[Low,0,1];", "certified genus-two level-2 basis"),
        ("rvk3pencil P[B,0,1];", "already defined"),
        ("show rvK3PencilParam[P](z,z);", "nonzero projective pair"),
        (
            "show rvK3PencilMap[P](zero[RVKummer(C)]);",
            "expected type RVKummerK3(C)",
        ),
        (
            "show rvK3PencilSection[P,0](rvK3PencilParam[P](o,o));",
            "must contain exactly one",
        ),
        (
            "show rvK3PencilSection[P,2](rvK3PencilParam[P](o,o));",
            "must contain exactly one",
        ),
        (
            "rvk3pencil Q[B,1,0]; rvK3PencilFiber[P](rvK3PencilParam[Q](o,o));",
            "expected type RVK3PencilBase(P)",
        ),
        (
            "rvunram G[5,1]{0,1}; show rvK3PencilParam[P](o,rvUnramElem[G]{1});",
            "expected type RVUnram(F)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F5} rvk3pencil P[B,0,1]; {program}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn runnable_example_constructs_a_fiber_and_a_section() {
    let output = stdout(example("hyperelliptic-k3-two-node-pencil"));
    assert!(output.starts_with("rvKummerDivPair = 0\nRVK3PencilBase(P)[[1], [2]]\n"));
    assert!(output.contains("rvK3PencilCensus[P] = base points=6, total=96,"));
}

fn digits(s: &str) -> Vec<u8> {
    s.bytes()
        .filter(u8::is_ascii_digit)
        .map(|digit| digit - b'0')
        .collect()
}
