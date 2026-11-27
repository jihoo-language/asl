mod common;

use std::collections::BTreeMap;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,t,z,z,z,o};
    let identity=zero[RVKummer(C)];
    let boundary=rvJacKummerPoint(rvJacPoint[C](z,z));
    let affine=rvJacKummerPoint(add(rvJacPoint[C](z,z),rvJacPoint[C](o,z)));
";

#[test]
fn exceptional_points_normalize_scaling_and_translation_by_the_base_point() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-resolution")),
        concat!(
            "RVKummerK3(C){base=[[0], [0], [0], [1]], direction=[[1], [1], [1], [0]]}\n",
            "RVKummerK3(C){base=[[1], [1], [0], [2]], direction=[[0], [1], [0], [1]]}\n",
        )
    );
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            let x=rvKummerExceptional(boundary,o,z,t,o);
            assert x==rvKummerExceptional(boundary,t,o,o,t);
            assert rvKummerDirection[0](x)==o;
            assert rvKummerDirection[1](x)==z;
            assert rvKummerDirection[2](x)==t;
            assert rvKummerDirection[3](x)==o;
            assert rvKummerBlowdown(x)==boundary;
            let y=rvKummerExceptional(affine,z,o,z,o);
            assert y==rvKummerExceptional(affine,o,z,z,o);
            "
        ))),
        ""
    );
}

#[test]
fn tests_full_hasse_conics_at_identity_boundary_and_affine_nodes() {
    for (name, pivot) in [("identity", 3), ("boundary", 1), ("affine", 0)] {
        let mut accepted = 0;
        for a in 0i32..3 {
            for b in 0i32..3 {
                for c in 0i32..3 {
                    let coordinates = [a, b, c];
                    if coordinates.iter().copied().find(|&x| x != 0) != Some(1) {
                        continue;
                    }
                    // These are the ε² coefficients in the homogeneous
                    // equation of y²=x⁵−x, in the three indicated affine charts.
                    let jet = match pivot {
                        3 => b * b - a * c,
                        1 => c * c + a * b,
                        0 => (a - c) * (a - c) + 2 * b * (a + c),
                        _ => unreachable!(),
                    };
                    let mut direction = coordinates.map(|x| ["z", "o", "t"][x as usize]).to_vec();
                    direction.insert(pivot, "z");
                    let output = interpret(&format!(
                        "{F3} let x=rvKummerExceptional({name},{}); assert rvKummerBlowdown(x)=={name};",
                        direction.join(",")
                    ));
                    if jet.rem_euclid(3) == 0 {
                        assert_eq!(stdout(output), "");
                        accepted += 1;
                    } else {
                        assert!(diagnostic(output).contains("quadratic tangent cone"));
                    }
                }
            }
        }
        assert_eq!(accepted, 4);
    }
}

#[test]
fn finite_type_matches_independent_k3_counts_and_fibers_over_f3_and_f9() {
    let output = stdout(interpret(
        r"
        rvunram F[3,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let t=rvUnramElem[F]{2};
        rvhyperelliptic C[F]{z,t,z,t,o,o};
        rvunram W[3,1]{1,0,1};
        rvjacbase L[C,W];
        rvKummerK3Points[C,1];
        show (\x:RVKummerK3(C)=>rvKummerBlowdown(x));
        show (\x:RVKummerK3(C)=>rvKummerIsExceptional(x));
        rvKummerK3Points[L,1];
        show (\x:RVKummerK3(L)=>rvKummerBlowdown(x));
        show (\x:RVKummerK3(L)=>rvKummerIsExceptional(x));
        ",
    ));
    let mut lines = output.lines();
    for (curve, q, singular, resolved) in [("C", 3, 8, 20), ("L", 9, 98, 134)] {
        assert_eq!(
            lines.next().unwrap(),
            format!(
                "rvKummerK3Points[{curve}, 1] = singular={singular}, nodes=4, resolved={resolved}"
            )
        );
        let blowdown = lines.next().unwrap();
        assert_eq!(
            blowdown.matches(&format!("RVKummerK3({curve})")).count(),
            resolved
        );
        let mut fibers = BTreeMap::<_, usize>::new();
        for (index, _) in blowdown.match_indices(&format!("=> RVKummer({curve})")) {
            let start = index + "=> ".len();
            let end = start + blowdown[start..].find("]]").unwrap() + 2;
            *fibers.entry(&blowdown[start..end]).or_default() += 1;
        }
        assert_eq!(fibers.len(), singular);
        assert_eq!(fibers.values().filter(|&&size| size == q + 1).count(), 4);
        assert_eq!(
            fibers.values().filter(|&&size| size == 1).count(),
            singular - 4
        );
        let exceptional = lines.next().unwrap();
        assert_eq!(exceptional.matches("=> S.hi").count(), 4 * (q + 1));
        assert_eq!(exceptional.matches("=> S.lo").count(), singular - 4);
    }
    assert!(lines.next().is_none());
}

#[test]
fn twist_only_points_fill_smoothly_and_have_no_direction() {
    let setup = r"
        rvunram F[3,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let t=rvUnramElem[F]{2};
        rvhyperelliptic C[F]{z,o,o,o,z,o};
        let point=rvKummer[C](o,o,t,t);
        let resolved=rvKummerSmooth(point);
    ";
    assert_eq!(
        stdout(interpret(&format!(
            "{setup} assert rvKummerLiftable(point)==S.lo; assert rvKummerIsExceptional(resolved)==S.lo; assert rvKummerBlowdown(resolved)==point; show resolved;"
        ))),
        "RVKummerK3(C){base=[[1], [1], [2], [2]]}\n"
    );
    for (source, expected) in [
        (
            "show rvKummerDirection[0](resolved);",
            "requires an exceptional point",
        ),
        (
            "show rvKummerExceptional(point,o,z,z,z);",
            "requires a Kummer node",
        ),
    ] {
        assert!(diagnostic(interpret(&format!("{setup}{source}"))).contains(expected));
    }
}

#[test]
fn rejects_invalid_directions_and_nominal_type_mismatches() {
    for (source, expected) in [
        ("show rvKummerSmooth(identity);", "requires a nonsingular Kummer point"),
        ("show rvKummerExceptional(identity,z,z,z,o);", "nonzero projective tangent direction"),
        ("show rvKummerExceptional(identity,z,z,z,z);", "nonzero projective tangent direction"),
        ("show rvKummerExceptional(identity,z,o,z,z);", "quadratic tangent cone"),
        ("show rvKummerDirection[4](rvKummerExceptional(identity,o,z,z,z));", "between 0 and 3"),
        ("show rvKummerBlowdown(identity);", "an RVKummerK3 point"),
        ("show rvKummerExceptional(zero[RVJac(C)],o,z,z,z);", "an RVKummer point"),
        ("rvunram G[3,1]{0,1}; show rvKummerExceptional(identity,rvUnramElem[G]{1},z,z,z);", "expected type RVUnram(F)"),
        ("rvhyperelliptic D[F]{z,t,z,z,z,o}; assert rvKummerExceptional(identity,o,z,z,z)==rvKummerExceptional(zero[RVKummer(D)],o,z,z,z);", "expected type RVKummerK3(C)"),
    ] {
        let error = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(error.contains(expected), "{error}");
    }
    assert!(diagnostic(interpret(&format!("{F3} show zero[RVKummerK3(C)];"))).contains("zero"));
}
