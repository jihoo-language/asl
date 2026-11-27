mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

const BASIS: &str = r"
    let r0=rvJacPoint[C](z,z);
    let r1=rvJacPoint[C](o,z);
    let r2=rvJacPoint[C](two,z);
    let r3=rvJacPoint[C](three,z);
    rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
    rvk3pencil P[B,0,1];
";

const F9: &str = r"
    rvunram F[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=neg(o);
    let three=rvUnramElem[F]{0,1};
    rvhyperelliptic C[F]{z,two,z,z,z,o};
";

fn split(p: i64) -> String {
    format!(
        r"
        rvunram F[{p},1]{{0,1}};
        let z=rvUnramElem[F]{{}};
        let o=rvUnramElem[F]{{1}};
        let two=rvUnramElem[F]{{2}};
        let three=rvUnramElem[F]{{3}};
        let c1=rvUnramElem[F]{{24}};
        let c2=rvUnramElem[F]{{-50}};
        let c3=rvUnramElem[F]{{35}};
        let c4=rvUnramElem[F]{{-10}};
        rvhyperelliptic C[F]{{z,c1,c2,c3,c4,o}};
        {BASIS}"
    )
}

#[test]
fn every_smooth_fiber_in_two_prime_fields_matches_the_resolved_census() {
    for prime in [13, 17] {
        let mut source = split(prime);
        source.push_str(
            "
            let inf=rvK3PencilParam[P](z,o);
            rvK3FiberQuartic[P](inf);
            rvK3FiberJacobian[P](inf);
            ",
        );
        for lambda in 0..prime {
            writeln!(
                source,
                "let p{lambda}=rvK3PencilParam[P](o,rvUnramElem[F]{{{lambda}}});
                 rvK3FiberQuartic[P](p{lambda});
                 rvK3FiberJacobian[P](p{lambda});"
            )
            .unwrap();
            if is_smooth(lambda) {
                writeln!(
                    source,
                    "show rvK3FiberJ[P](p{lambda});
                     rvK3PencilFiber[P](p{lambda});
                     rvK3FiberPoints[P,1](p{lambda});
                     rvk3fiberjac E{lambda}[P](p{lambda});
                     rvECPointCount[E{lambda}];"
                )
                .unwrap();
            }
        }
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        assert!(lines
            .next()
            .unwrap()
            .ends_with("status=singular (zero binary quartic)"));
        assert!(lines.next().unwrap().contains("discriminant=[0]"));
        for lambda in 0..prime {
            let expected = expected_quartic(prime, lambda);
            let report = lines.next().unwrap();
            let actual = report
                .split_once("F4=")
                .unwrap()
                .1
                .split_once(", status=")
                .unwrap()
                .0;
            assert_eq!(scalars(actual), expected);
            let jacobian = lines.next().unwrap();
            let [a2, a4, a6] = long_coefficients(prime, expected);
            assert!(jacobian.contains(&format!("a2=[{a2}], a4=[{a4}], a6=[{a6}]")));
            if is_smooth(lambda) {
                assert!(report.ends_with("status=smooth"));
                let discriminant = cubic_discriminant(prime, a2, a4, a6);
                assert_ne!(discriminant, 0);
                assert!(
                    jacobian.ends_with(&format!("discriminant=[{discriminant}], status=smooth"))
                );
                let c4 = (16 * (a2 * a2 - 3 * a4)).rem_euclid(prime);
                let expected_j = (c4.pow(3) * inverse(prime, discriminant)).rem_euclid(prime);
                assert_eq!(scalars(lines.next().unwrap()), [expected_j]);
                let direct = count(lines.next().unwrap());
                assert_eq!(count(lines.next().unwrap()), direct);
                assert_eq!(count(lines.next().unwrap()), direct);
            } else {
                assert!(report.contains("status=singular"));
                assert!(jacobian.contains("discriminant=[0]"));
            }
        }
        assert!(lines.next().is_none());
    }
}

#[test]
fn characteristic_three_uses_the_long_jacobian_for_all_four_smooth_fibers() {
    let mut first = format!("{F9}{BASIS}");
    for (index, value) in ["1,1", "1,2", "2,1", "2,2"].into_iter().enumerate() {
        writeln!(
            first,
            "let p{index}=rvK3PencilParam[P](o,rvUnramElem[F]{{{value}}});
             rvK3FiberJacobian[P](p{index});"
        )
        .unwrap();
    }
    let reports = stdout(interpret(&first));
    let mut second = first;
    for (index, report) in reports.lines().enumerate() {
        assert!(report.ends_with("status=smooth"), "{report}");
        let [a6, a4, a2] = ["a6=[", "a4=[", "a2=["].map(|label| {
            report
                .split_once(label)
                .unwrap()
                .1
                .split_once(']')
                .unwrap()
                .0
        });
        writeln!(
            second,
            "rvhyperelliptic E{index}[F]{{rvUnramElem[F]{{{a6}}},rvUnramElem[F]{{{a4}}},rvUnramElem[F]{{{a2}}},o}};
             rvK3FiberPoints[P,1](p{index});
             rvK3PencilFiber[P](p{index});
             rvHWWPoints[E{index},1];
             rvK3FiberPoints[P,2](p{index});
             rvHWWPoints[E{index},2];"
        )
        .unwrap();
    }
    let result = stdout(interpret(&second));
    let mut lines = result.lines().skip(4);
    for _ in 0..4 {
        let first = count(lines.next().unwrap());
        assert_eq!(count(lines.next().unwrap()), first);
        assert_eq!(count(lines.next().unwrap()), first);
        assert_eq!(count(lines.next().unwrap()), count(lines.next().unwrap()));
    }
    assert!(lines.next().is_none());

    let error = diagnostic(interpret(&format!("{second} rvk3fiberjac Short[P](p0);")));
    assert!(
        error.contains("requires characteristic at least five"),
        "{error}"
    );
}

#[test]
fn coordinate_changes_preserve_j_and_extension_counts() {
    let setup = split(13);
    let [a2, a4, a6] = long_coefficients(13, expected_quartic(13, 1));
    let output = stdout(interpret(&format!(
        r"{setup}
        rvk3pencil Reverse[B,1,0];
        rvk3pencil Translated[B,4,5];
        let p=rvK3PencilParam[P](o,o);
        let r=rvK3PencilParam[Reverse](two,two);
        let x=rvK3PencilSection[P,4](p);
        let q=rvK3PencilMap[Translated](rvKummerK3Translate[B,4](x));
        assert rvK3FiberJ[P](p)==rvK3FiberJ[Reverse](r);
        assert rvK3FiberJ[P](p)==rvK3FiberJ[Translated](q);
        rvK3FiberQuartic[P](p);
        rvK3FiberQuartic[Reverse](r);
        rvhyperelliptic Long[F]{{rvUnramElem[F]{{{a6}}},rvUnramElem[F]{{{a4}}},rvUnramElem[F]{{{a2}}},o}};
        rvK3FiberPoints[P,2](p);
        rvK3FiberPoints[Translated,2](q);
        rvHWWPoints[Long,2];
        "
    )));
    let mut lines = output.lines();
    assert_eq!(
        lines.next().unwrap().split_once(" = ").unwrap().1,
        lines.next().unwrap().split_once(" = ").unwrap().1
    );
    let value = count(lines.next().unwrap());
    assert_eq!(count(lines.next().unwrap()), value);
    assert_eq!(count(lines.next().unwrap()), value);
    assert!(lines.next().is_none());
}

#[test]
fn larger_field_counts_only_the_curve() {
    let setup = split(101);
    let source = format!(
        r"{setup}
        let p=rvK3PencilParam[P](o,o);
        rvk3fiberjac E[P](p);
        rvK3FiberPoints[P,1](p);
        rvECPointCount[E];
        assert add(rvECInfinity[E],rvECInfinity[E])==zero[RVECPoint(E)];
        "
    );
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    assert_eq!(count(lines.next().unwrap()), count(lines.next().unwrap()));
    let error = diagnostic(interpret(&format!("{source} rvK3PencilCensus[P];")));
    assert!(error.contains("limit"), "{error}");
}

#[test]
fn rejects_singular_fibers_invalid_parameters_and_overflow() {
    let setup = split(13);
    let setup = format!(
        "{setup} let p=rvK3PencilParam[P](o,o); let central=rvK3PencilParam[P](o,z); let repeated=rvK3PencilParam[P](o,two);"
    );
    for (input, message) in [
        ("show rvK3FiberJ[P](central);", "fiber is singular"),
        ("show rvK3FiberJ[P](repeated);", "fiber is singular"),
        ("rvK3FiberPoints[P,1](repeated);", "fiber is singular"),
        ("rvk3fiberjac E[P](central);", "fiber is singular"),
        (
            "rvK3FiberPoints[P,0](p);",
            "extension degree must be positive",
        ),
        (
            "rvK3FiberPoints[P,340282366920938463463374607431768211455](p);",
            "exceeds unsigned 128-bit range",
        ),
        ("rvK3FiberQuartic[P](o);", "expected type RVK3PencilBase(P)"),
        (
            "rvk3pencil Q[B,1,0]; show rvK3FiberJ[P](rvK3PencilParam[Q](o,o));",
            "expected type RVK3PencilBase(P)",
        ),
        (
            "rvk3fiberjac E[P](p); rvk3fiberjac E[P](p);",
            "elliptic curve is already defined",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {input}")));
        assert!(error.contains(message), "{error}");
    }
}

#[test]
fn example_reports_equations_counts_and_a_central_double_plane() {
    let output = stdout(example("hyperelliptic-k3-pencil-fiber-jacobian"));
    assert!(output.contains("status=smooth"));
    assert!(output.contains(concat!(
        "A2=[[1], [0], [0]], B2=[[0], [11], [0]], C2=[[0], [0], [1]], ",
        "F4=[[0], [0], [0], [0], [0]], status=singular (zero binary quartic)"
    )));
    assert!(output.contains("rvK3FiberPoints[P, 1] = 12\n"));
    assert!(output.contains("rvK3FiberPoints[P, 2] = 192\n"));
    assert!(output.contains("rvECPointCount[E] = 12, trace=2\n"));
}

fn is_smooth(lambda: i64) -> bool {
    ![0, 2, 3, 4, 6, 8, 12].contains(&lambda)
}

fn expected_quartic(p: i64, lambda: i64) -> [i64; 5] {
    let mut quartic = vec![16 * lambda];
    for alpha in 1..=4 {
        let mut next = vec![0; quartic.len() + 1];
        for (i, &coefficient) in quartic.iter().enumerate() {
            next[i] = (next[i] - alpha * coefficient).rem_euclid(p);
            next[i + 1] = (next[i + 1] + (lambda + alpha * alpha) * coefficient).rem_euclid(p);
        }
        quartic = next;
    }
    quartic.try_into().unwrap()
}

fn long_coefficients(p: i64, [a, b, c, d, e]: [i64; 5]) -> [i64; 3] {
    [c, b * d - 4 * a * e, a * d * d + b * b * e - 4 * a * c * e].map(|value| value.rem_euclid(p))
}

fn cubic_discriminant(p: i64, c: i64, a: i64, b: i64) -> i64 {
    (16 * (c * c * a * a - 4 * a.pow(3) - 4 * c.pow(3) * b - 27 * b * b + 18 * c * a * b))
        .rem_euclid(p)
}

fn inverse(p: i64, value: i64) -> i64 {
    (1..p).find(|n| n * value % p == 1).unwrap()
}

fn scalars(report: &str) -> Vec<i64> {
    report
        .split(|c: char| !c.is_ascii_digit())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}

fn count(report: &str) -> u128 {
    report
        .split_once(" = ")
        .unwrap()
        .1
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
