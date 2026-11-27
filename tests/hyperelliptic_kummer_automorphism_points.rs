mod common;

use std::collections::{BTreeMap, BTreeSet};

use common::{diagnostic, example, interpret, stdout};

const F5: &str = r"
    rvunram F[5,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let th=rvUnramElem[F]{3};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    let r0=rvJacPoint[C](z,z);
    let r1=rvJacPoint[C](o,z);
    let r2=rvJacPoint[C](t,z);
    let r3=rvJacPoint[C](th,z);
    rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
";

const F9: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase D[C,W];
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    let r0=rvJacPoint[D](zw,zw);
    let r1=rvJacPoint[D](ow,zw);
    let r2=rvJacPoint[D](neg(ow),zw);
    let r3=rvJacPoint[D](i,zw);
    rvjactorsionbasis B[D,2](neg(ow)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
";

#[test]
fn runnable_example_exhibits_period_nine_and_hyperbolic_finite_cycles() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-automorphism-points")),
        concat!(
            "rvKummerK3Points[C, 1] = singular=16, nodes=16, resolved=96\n",
            "rvK3AutoPeriod[Mix] = 9\n",
            "rvK3AutoCycles[Mix] = total=96, fixed=0, cycles=[(2, 6), (4, 12), (9, 4)]\n",
            "rvK3AutoCycles[Hyperbolic] = total=96, fixed=8, cycles=[(1, 8), (2, 2), (3, 18), (5, 6)]\n",
            "rvK3Entropy[Hyperbolic] = characteristic=[1, -31, 312, -1672, 5628, -12740, 19656, -19448, 8294, 8294, -19448, 19656, -12740, 5628, -1672, 312, -31, 1], cyclotomic=[(1, 14), (2, 1)], reciprocal quadratic=[1, -18, 1], dynamical degree=largest real root greater than 1 of [1, -18, 1]\n",
        ),
    );
}

#[test]
fn whole_resolution_action_and_census_match_explicit_nested_generators() {
    let output = stdout(interpret(&format!(
        r"{F9}
        rvk3autom Mix[B]{{translate(3),switch(5),project(1)}};
        rvk3autom Reverse[B]{{project(1),switch(5),translate(3)}};
        rvk3autom Idle[B]{{}};
        let center=rvKummerTranslate[B,1](zero[RVKummer(D)]);
        let explicit=\x:RVKummerK3(D)=>rvKummerNodeInvolution(center,rvKummerK3Translate[B,5](rvKummerSwitch(rvKummerK3Translate[B,3](x))));
        assert (\x:RVKummerK3(D)=>rvK3AutoPoint[Mix](x))==explicit;
        assert (\x:RVKummerK3(D)=>rvK3AutoIter[Mix,-1](x))==(\x:RVKummerK3(D)=>rvK3AutoPoint[Reverse](x));
        assert (\x:RVKummerK3(D)=>rvK3AutoIter[Mix,-1](rvK3AutoPoint[Mix](x)))==(\x:RVKummerK3(D)=>x);
        assert (\x:RVKummerK3(D)=>rvK3AutoIter[Mix,0](x))==(\x:RVKummerK3(D)=>x);
        assert (\x:RVKummerK3(D)=>rvK3AutoIter[Mix,2](x))==(\x:RVKummerK3(D)=>explicit(explicit(x)));
        rvKummerK3Points[D,1];
        rvK3AutoCycles[Mix];
        rvK3AutoCycles[Reverse];
        rvK3AutoCycles[Idle];
        show explicit;
        "
    )));
    let mut lines = output.lines();
    assert_eq!(
        lines.next().unwrap(),
        "rvKummerK3Points[D, 1] = singular=104, nodes=16, resolved=248"
    );
    assert_eq!(
        lines.next().unwrap(),
        "rvK3AutoCycles[Mix] = total=248, fixed=0, cycles=[(2, 6), (3, 20), (4, 24), (10, 8)]"
    );
    assert_eq!(
        lines.next().unwrap(),
        "rvK3AutoCycles[Reverse] = total=248, fixed=0, cycles=[(2, 6), (3, 20), (4, 24), (10, 8)]"
    );
    assert_eq!(
        lines.next().unwrap(),
        "rvK3AutoCycles[Idle] = total=248, fixed=248, cycles=[(1, 248)]"
    );
    let mapping = named_map(lines.next().unwrap());
    assert_eq!(mapping.len(), 248);
    let cycles = independently_count_cycles(&mapping);
    assert_eq!(cycles, BTreeMap::from([(2, 6), (3, 20), (4, 24), (10, 8)]));
    assert_eq!(
        cycles
            .iter()
            .map(|(length, count)| length * count)
            .sum::<usize>(),
        248
    );
    assert!(lines.next().is_none());
}

#[test]
fn huge_positive_and_negative_iterates_and_periods_use_the_first_return() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            rvk3autom Mix[B]{{translate(3),switch(5),project(1)}};
            rvk3autom Idle[B]{{}};
            let x=rvKummerExceptional(rvJacKummerPoint(r1),o,z,z,z);
            let two_cycle=rvKummerTropeParam[B,0](z,o);
            assert rvK3AutoIter[Mix,170141183460469231731687303715884105727](x)==rvK3AutoPoint[Mix](x);
            assert rvK3AutoIter[Mix,-170141183460469231731687303715884105728](x)==rvK3AutoIter[Mix,-2](x);
            assert rvK3AutoIter[Mix,18](x)==x;
            assert rvK3AutoIter[Idle,-170141183460469231731687303715884105728](x)==x;
            rvK3AutoPeriod[Mix](x);
            rvK3AutoPeriod[Mix](rvK3AutoPoint[Mix](x));
            rvK3AutoPeriod[Mix](two_cycle);
            rvK3AutoPeriod[Idle](x);
        "
        ))),
        "rvK3AutoPeriod[Mix] = 9\nrvK3AutoPeriod[Mix] = 9\nrvK3AutoPeriod[Mix] = 2\nrvK3AutoPeriod[Idle] = 1\n",
    );
}

#[test]
fn individual_orbits_do_not_enumerate_large_point_spaces() {
    let setup = format!(
        r"{F5}
        rvunram U[5,1]{{1,4,0,0,0,1}};
        rvjacbase L[C,U];
        let zu=rvUnramElem[U]{{}};
        let ou=rvUnramElem[U]{{1}};
        rvjactorsionbasis Large[L,2](neg(ou)){{rvJacExtend[L](r0),rvJacExtend[L](add(add(r0,r1),r2)),rvJacExtend[L](r1),rvJacExtend[L](add(add(r0,r1),r3))}};
        let x=rvKummerExceptional(rvJacKummerPoint(rvJacExtend[L](r1)),ou,zu,zu,zu);
        rvk3autom T[Large]{{translate(3)}};
        assert rvK3AutoPoint[T](x)==rvKummerK3Translate[Large,3](x);
        assert rvK3AutoIter[T,100000000000000000000000000000000001](x)==rvK3AutoPoint[T](x);
        rvK3AutoPeriod[T](x);
        "
    );
    assert_eq!(stdout(interpret(&setup)), "rvK3AutoPeriod[T] = 2\n");
    let problem = diagnostic(interpret(&format!("{setup} rvK3AutoCycles[T];")));
    assert!(problem.contains("computation limit"), "{problem}");
}

#[test]
fn point_api_checks_word_names_curve_types_and_iteration_range() {
    let setup = format!("{F5} rvk3autom A[B]{{project(0)}}; let x=rvKummerTropeParam[B,0](z,o);");
    for (source, expected) in [
        ("show rvK3AutoPoint[Missing](x);", "unknown K3 automorphism word"),
        ("rvK3AutoCycles[Missing];", "unknown K3 automorphism word"),
        ("show rvK3AutoPoint[A](zero[RVKummer(C)]);", "expected type RVKummerK3(C)"),
        ("rvhyperelliptic Twin[F]{z,neg(o),z,z,z,o}; let y=rvKummerExceptional(zero[RVKummer(Twin)],z,z,o,z); rvK3AutoPeriod[A](y);", "expected type RVKummerK3(C), found RVKummerK3(Twin)"),
        ("show rvK3AutoIter[A,170141183460469231731687303715884105728](x);", "signed i128 range"),
    ] {
        let problem = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}

fn named_map(output: &str) -> BTreeMap<&str, &str> {
    let mut remaining = output.strip_prefix('{').unwrap().strip_suffix('}').unwrap();
    let mut map = BTreeMap::new();
    while !remaining.is_empty() {
        let (from, after_arrow) = remaining.split_once(" => ").unwrap();
        let end = after_arrow.find('}').unwrap() + 1;
        let (to, rest) = after_arrow.split_at(end);
        assert!(map.insert(from, to).is_none());
        remaining = rest.strip_prefix(", ").unwrap_or(rest);
    }
    map
}

fn independently_count_cycles<'a>(mapping: &BTreeMap<&'a str, &'a str>) -> BTreeMap<usize, usize> {
    let mut seen = BTreeSet::new();
    let mut counts = BTreeMap::new();
    for &point in mapping.keys() {
        if seen.contains(point) {
            continue;
        }
        let mut image = point;
        let mut length = 0;
        while seen.insert(image) {
            length += 1;
            image = mapping[image];
        }
        assert_eq!(image, point);
        *counts.entry(length).or_default() += 1;
    }
    counts
}
