mod common;

use std::collections::{BTreeMap, BTreeSet};

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,o,o,o,z,o};
";

const F5: &str = r"
    rvunram F[5,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let m=rvUnramElem[F]{4};
    rvhyperelliptic C[F]{o,o,z,z,z,o};
";

#[test]
fn projective_scaling_repeated_support_and_twist_on_the_boundary() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-surface")),
        "RVKummer(C)[[1], [0], [0], [4]]\nRVKummer(C)[[0], [1], [1], [1]]\n",
    );
}

#[test]
fn identity_and_rational_nodes_have_unique_lifts() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F5}
            rvhyperelliptic N[F]{{z,m,z,z,z,o}};
            let id=rvKummer[N](z,z,z,m);
            assert id==zero[RVKummer(N)];
            assert id==rvJacKummerPoint(zero[RVJac(N)]);
            assert rvKummerLift(id)==zero[RVJac(N)];
            let branch=rvJacPoint[N](o,z);
            let k=rvJacKummerPoint(branch);
            assert k==rvKummer[N](z,o,o,o);
            assert rvKummerLift(k)==branch;
            let node=add(rvJacPoint[N](z,z),branch);
            let image=rvKummer[N](o,o,z,m);
            assert rvJacKummerPoint(node)==image;
            assert rvKummerLiftable(image)==S.hi;
            assert rvKummerLift(image)==node;
            assert neg(node)==node;
        "
        ))),
        "",
    );
}

#[test]
fn extension_field_scaling_and_square_root_use_the_named_basis() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=rvUnramElem[F]{2};
            let alpha=rvUnramElem[F]{0,1};
            rvhyperelliptic C[F]{t,z,z,z,z,o};
            let d=rvJacPoint[C](z,alpha);
            let point=rvKummer[C](z,alpha,z,z);
            assert rvJacKummerPoint(d)==point;
            assert rvKummerCoord[1](point)==o;
            assert rvKummerLift(point)==d;
        ",
        )),
        "",
    );
}

#[test]
fn all_small_projections_and_fibers_match_finite_space_enumeration() {
    for (setup, surface_size, rational_classes, singleton_fibers, paired_fibers) in
        [(F3, 12, 20, 2, 9), (F5, 36, 36, 4, 16)]
    {
        let output = stdout(interpret(&format!(
            r"{setup}
            let project=\d:RVJac(C)=>rvJacKummerPoint(d);
            let projectNegative=\d:RVJac(C)=>rvJacKummerPoint(neg(d));
            let projectLift=\d:RVJac(C)=>rvJacKummerPoint(rvKummerLift(project(d)));
            assert project==projectNegative;
            assert project==projectLift;
            let liftable=\k:RVKummer(C)=>rvKummerLiftable(k);
            show project;
            show liftable;
        "
        )));
        let mut lines = output.lines();
        let project = lines.next().unwrap();
        let liftable = lines.next().unwrap();
        assert!(lines.next().is_none());
        assert_eq!(project.matches("RVJac(C){").count(), rational_classes);
        assert_eq!(liftable.matches("RVKummer(C)").count(), surface_size);
        let mut multiplicities = BTreeMap::<_, usize>::new();
        for image in project.match_indices("=> RVKummer(C)") {
            let start = image.0 + "=> ".len();
            let end = project[start..].find("]]").unwrap() + start + 2;
            *multiplicities.entry(&project[start..end]).or_default() += 1;
        }
        let marked: BTreeSet<_> = liftable
            .match_indices("RVKummer(C)")
            .filter_map(|(start, _)| {
                let end = liftable[start..].find("]]").unwrap() + start + 2;
                liftable[end..]
                    .starts_with(" => S.hi")
                    .then_some(&liftable[start..end])
            })
            .collect();
        assert_eq!(
            multiplicities.keys().copied().collect::<BTreeSet<_>>(),
            marked
        );
        assert_eq!(
            multiplicities.values().filter(|&&n| n == 1).count(),
            singleton_fibers
        );
        assert_eq!(
            multiplicities.values().filter(|&&n| n == 2).count(),
            paired_fibers
        );
        assert_eq!(multiplicities.len(), singleton_fibers + paired_fibers);
    }
}

#[test]
fn nonsquare_quadratic_recovery_is_a_rational_kummer_point() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            let point=rvKummer[C](o,o,t,t);
            let identity=\k:RVKummer(C)=>k;
            assert identity(point)==point;
            assert rvKummerLiftable(point)==S.lo;
            show point;
        "
        ))),
        "RVKummer(C)[[1], [1], [2], [2]]\n",
    );
}

#[test]
fn diagnoses_invalid_coordinates_lifts_and_nominal_types() {
    for (source, expected) in [
        ("show rvKummer[C](z,z,z,z);", "nonzero projective quadruple"),
        ("show rvKummer[C](z,o,o,z);", "Kummer quartic"),
        ("show rvKummer[C](o,z,z,z);", "Kummer quartic"),
        (
            "show rvKummerLift(rvKummer[C](o,o,t,t));",
            "nonsquare recovery data",
        ),
        (
            "show rvKummerCoord[4](zero[RVKummer(C)]);",
            "between 0 and 3",
        ),
        (
            "rvhyperelliptic D[F]{z,o,o,o,z,o}; assert zero[RVKummer(C)]==zero[RVKummer(D)];",
            "expected type RVKummer(C)",
        ),
        (
            "rvunram E[3,1]{0,1}; show rvKummer[C](rvUnramElem[E]{1},o,t,t);",
            "rvKummer coordinate",
        ),
        (
            "rvhyperelliptic E[F]{o,o,z,o}; show rvJacKummerPoint(zero[RVJac(E)]);",
            "smooth monic quintic",
        ),
        (
            "show rvJacKummerPoint(zero[RVKummer(C)]);",
            "RVJac divisor class",
        ),
        ("show rvKummerLift(zero[RVJac(C)]);", "RVKummer point"),
    ] {
        let error = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
