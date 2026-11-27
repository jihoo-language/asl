mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,t,z,z,z,o};
    rvunram W[3,1]{1,0,1};
    rvjacbase D[C,W];
    let dw=zero[RVKummer(D)];
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let tw=rvUnramElem[W]{2};
    let alpha=rvUnramElem[W]{0,1};
    let tangent=rvKummerExceptional(dw,ow,alpha,tw,zw);
";

const F9: &str = r"
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[W]{};
    let o=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    let m=neg(o);
    let ni=neg(i);
    rvhyperelliptic C[W]{z,m,z,z,z,o};
    let t0=rvJacPoint[C](z,z);
    let t1=rvJacPoint[C](o,z);
    let tm=rvJacPoint[C](m,z);
    let ti=rvJacPoint[C](i,z);
    let e2=add(add(t0,t1),tm);
    let f2=add(add(t0,t1),ti);
    rvjactorsionbasis B[C,2](m){t0,e2,t1,f2};
    let id=zero[RVKummer(C)];
";

#[test]
fn extension_descends_both_coordinates_and_tangent_and_frobenius_is_relative() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            let r=rvKummerExceptional(zero[RVKummer(C)],o,o,o,z);
            let baseImage=rvKummerK3Extend[D](r);
            assert rvKummerK3Descend[C](baseImage)==r;
            assert rvKummerBlowdown(baseImage)==rvKummerExtend[D](rvKummerBlowdown(r));
            assert rvKummerK3Frob[C,1](baseImage)==baseImage;
            assert rvKummerK3Frob[C,1](tangent)==rvKummerExceptional(dw,ow,neg(alpha),tw,zw);
            assert rvKummerK3Frob[C,2](tangent)==tangent;
            assert rvKummerDescend[C](rvKummerBlowdown(tangent))==zero[RVKummer(C)];

            rvunram V[3,1]{{2,0,1,0,1}};
            rvjacbase E[D,V];
            let extended=rvKummerK3Extend[E](tangent);
            assert rvKummerDirection[1](extended)==rvUnramElem[V]{{2,0,1,0}};
            assert rvKummerK3Frob[D,1](extended)==extended;
            assert rvKummerK3Descend[D](extended)==tangent;

            rvhyperelliptic T[F]{{z,o,o,o,z,o}};
            rvjacbase U[T,W];
            let point=rvKummer[T](o,o,t,t);
            assert rvKummerLiftable(point)==S.lo;
            let smooth=rvKummerSmooth(point);
            assert rvKummerK3Extend[U](smooth)==rvKummerSmooth(rvKummerExtend[U](point));
            assert rvKummerK3Descend[T](rvKummerK3Extend[U](smooth))==smooth;
            "
        ))),
        ""
    );
    assert!(diagnostic(interpret(&format!(
        "{F3} show rvKummerK3Descend[C](tangent);"
    )))
    .contains("tangent coordinate lies outside the requested subfield"));
}

#[test]
fn translates_the_entire_resolution_and_preserves_the_trope_action() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F9}
            assert (\x:RVKummerK3(C)=>rvKummerK3Translate[B,5](rvKummerK3Translate[B,3](x)))
                ==(\x:RVKummerK3(C)=>rvKummerK3Translate[B,6](x));
            assert (\x:RVKummerK3(C)=>rvKummerBlowdown(rvKummerK3Translate[B,13](x)))
                ==(\x:RVKummerK3(C)=>rvKummerTranslate[B,13](rvKummerBlowdown(x)));
            assert (\x:RVKummerK3(C)=>rvKummerOnTrope[B,2](rvKummerK3Translate[B,5](x)))
                ==(\x:RVKummerK3(C)=>rvKummerOnTrope[B,7](x));
            "
        ))),
        ""
    );
}

fn rendered_points(output: &str) -> impl Iterator<Item = (&str, &str)> {
    output.match_indices("RVKummerK3(C){").map(|(start, _)| {
        let end = start + output[start..].find('}').unwrap() + 1;
        (&output[start..end], &output[end..])
    })
}

#[test]
fn strict_tropes_are_disjoint_complete_and_have_six_distinct_incident_directions() {
    let branches = [
        ("z", "o", 0),
        ("o", "z", 1),
        ("o", "o", 4),
        ("o", "m", 7),
        ("o", "i", 13),
        ("o", "ni", 15),
    ];
    let mut source = F9.to_owned();
    for a in 0..16 {
        for (u, v, branch) in branches {
            let mask = a ^ branch;
            writeln!(source, "assert rvKummerBlowdown(rvKummerTropeParam[B,{a}]({u},{v}))==rvKummerTranslate[B,{mask}](id);").unwrap();
        }
        writeln!(
            source,
            "show (\\x:RVKummerK3(C)=>rvKummerOnTrope[B,{a}](x));"
        )
        .unwrap();
        writeln!(
            source,
            "show (\\s:RVUnram(W)=>rvKummerTropeParam[B,{a}](o,s));"
        )
        .unwrap();
        writeln!(source, "show rvKummerTropeParam[B,{a}](z,o);").unwrap();
    }
    source.push_str(r"show (\x:RVKummerK3(C)=>rvKummerIsExceptional(x));");
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    let mut all_trope_points = BTreeSet::new();
    let mut incident = BTreeMap::<_, usize>::new();
    for _ in 0..16 {
        let marked: BTreeSet<_> = rendered_points(lines.next().unwrap())
            .filter_map(|(point, after)| after.starts_with(" => S.hi").then_some(point))
            .collect();
        let parameterized: BTreeSet<_> = [lines.next().unwrap(), lines.next().unwrap()]
            .into_iter()
            .flat_map(rendered_points)
            .map(|(point, _)| point)
            .collect();
        assert_eq!(parameterized.len(), 10);
        assert_eq!(marked, parameterized);
        let exceptional: Vec<_> = marked
            .iter()
            .filter_map(|point| point.split_once(", direction="))
            .collect();
        assert_eq!(exceptional.len(), 6);
        for (node, _) in exceptional {
            *incident.entry(node).or_default() += 1;
        }
        for point in marked {
            assert!(all_trope_points.insert(point), "{point}");
        }
    }
    assert_eq!(incident.len(), 16);
    assert!(incident.values().all(|&count| count == 6));
    let mut fibers = BTreeMap::<_, usize>::new();
    for (point, after) in rendered_points(lines.next().unwrap()) {
        if after.starts_with(" => S.hi") {
            let (node, _) = point.split_once(", direction=").unwrap();
            *fibers.entry(node).or_default() += 1;
        }
    }
    assert_eq!(fibers.len(), 16);
    assert_eq!(
        fibers.keys().collect::<Vec<_>>(),
        incident.keys().collect::<Vec<_>>()
    );
    assert!(fibers.values().all(|&count| count == 10));
    assert!(lines.next().is_none());
}

#[test]
fn projective_parameter_scaling_and_tangent_membership_are_exact() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-resolved-tropes")),
        "S.hi\nS.lo\n"
    );
    assert_eq!(
        stdout(interpret(&format!(
            r"{F9}
            let x=add(o,i);
            let point=rvKummerTropeParam[B,0](o,x);
            assert rvKummerTropeParam[B,0](i,rvUnramMul(i,x))==point;
            assert rvKummerIsExceptional(point)==S.lo;
            assert rvKummerLiftable(rvKummerBlowdown(point))==S.lo;
            let infinity=rvKummerTropeParam[B,0](z,i);
            assert infinity==rvKummerExceptional(id,z,z,o,z);
            assert rvKummerTropeParam[B,0](o,i)
                ==rvKummerExceptional(rvJacKummerPoint(ti),z,z,o,neg(i));
            let differentDirection=rvKummerExceptional(id,o,z,z,z);
            assert rvKummerOnTrope[B,0](infinity)==S.hi;
            assert rvKummerOnTrope[B,0](differentDirection)==S.lo;
            "
        ))),
        ""
    );
}

#[test]
fn requires_nominal_fields_and_curves_and_valid_masks_and_parameters() {
    for (source, expected) in [
        ("show rvKummerTropeParam[B,0](z,z);", "nonzero projective pair"),
        ("show rvKummerTropeParam[B,16](o,z);", "between 0 and 15"),
        ("show rvKummerOnTrope[B,0](id);", "expected type RVKummerK3(C)"),
        ("rvunram V[3,1]{1,0,1}; show rvKummerTropeParam[B,0](rvUnramElem[V]{1},z);", "expected type RVUnram(W)"),
        ("rvhyperelliptic D[W]{z,m,z,z,z,o}; let x=rvKummerExceptional(zero[RVKummer(D)],z,z,o,z); show rvKummerK3Translate[B,0](x);", "expected type RVKummerK3(C)"),
    ] {
        let error = diagnostic(interpret(&format!("{F9}{source}")));
        assert!(error.contains(expected), "{error}");
    }
    for (source, expected) in [
        ("show rvKummerK3Extend[D](zero[RVKummer(C)]);", "an RVKummerK3 point"),
        ("rvhyperelliptic E[W]{ow,tw,zw,zw,zw,ow}; let x=rvKummerExceptional(zero[RVKummer(C)],z,z,o,z); show rvKummerK3Extend[E](x);", "not the declared base extension"),
    ] {
        let error = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
