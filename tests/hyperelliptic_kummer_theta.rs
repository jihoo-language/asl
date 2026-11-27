mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let three=rvUnramElem[F]{3};
    let four=rvUnramElem[F]{4};
    rvhyperelliptic C[F]{z,three,neg(o),z,four,o};
    let t0=rvJacPoint[C](z,z);
    let t1=rvJacPoint[C](o,z);
    let t2=rvJacPoint[C](two,z);
    let t3=rvJacPoint[C](three,z);
    let e1=t0;
    let e2=add(add(t0,t1),t2);
    let f1=t1;
    let f2=add(add(t0,t1),t3);
    rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};
";

fn marked_classes() -> String {
    let mut source = format!("{F7}\nlet n0=zero[RVJac(C)];\n");
    let generators = ["e1", "e2", "f1", "f2"];
    for mask in 1_usize..16 {
        let generator = generators[mask.trailing_zeros() as usize];
        let rest = mask & (mask - 1);
        writeln!(source, "let n{mask}=add(n{rest},{generator});").unwrap();
    }
    source
}

#[test]
fn every_translation_agrees_with_cantor_on_the_full_group_and_on_a_twist() {
    let mut source = marked_classes();
    source.push_str("rvJacCount[C];\n");
    for mask in 0..16 {
        writeln!(
            source,
            r"assert (\d:RVJac(C)=>rvKummerTranslate[B,{mask}](rvJacKummerPoint(d)))
                ==(\d:RVJac(C)=>rvJacKummerPoint(add(d,n{mask})));"
        )
        .unwrap();
    }
    source.push_str(
        r"
        // T is the monic twist with x_C=-x_T. Its Kummer coordinates
        // change by (w,x,y,z) ↦ (w,-x,y,-z).
        rvhyperelliptic T[F]{z,three,o,z,three,o};
        let back=\k:RVKummer(T)=>rvKummer[C](
            rvKummerCoord[0](k),neg(rvKummerCoord[1](k)),
            rvKummerCoord[2](k),neg(rvKummerCoord[3](k)));
        let p=rvJacPoint[T](o,o);
        let k=back(rvJacKummerPoint(p));
        let t=rvJacPoint[T](neg(o),z);
        let expected=back(rvJacKummerPoint(add(p,t)));
        assert rvKummerLiftable(k)==S.lo;
        assert rvKummerLiftable(expected)==S.lo;
        assert rvKummerTranslate[B,4](k)==expected;
        assert rvKummerTranslate[B,4](expected)==k;
        ",
    );
    assert_eq!(stdout(interpret(&source)), "rvJacCount[C] = 48\n");
}

#[test]
fn normalized_tropes_and_nodes_form_the_sixteen_six_configuration() {
    let mut source = marked_classes();
    for mask in 0..16 {
        writeln!(source, "let k{mask}=rvJacKummerPoint(n{mask});").unwrap();
        for index in 0..4 {
            writeln!(source, "show rvKummerCoord[{index}](k{mask});").unwrap();
        }
    }
    for mask in 0..16 {
        for index in 0..4 {
            writeln!(source, "show rvKummerTropeCoeff[B,{mask},{index}];").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let coordinates = output
        .lines()
        .map(|line| {
            line.strip_prefix("RVUnram(F)[")
                .unwrap()
                .strip_suffix(']')
                .unwrap()
                .parse::<u64>()
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(coordinates.len(), 2 * 16 * 4);
    let (nodes, tropes) = coordinates.split_at(16 * 4);
    let incidents: Vec<Vec<bool>> = tropes
        .chunks_exact(4)
        .map(|plane| {
            assert_eq!(
                plane.iter().find(|&&coefficient| coefficient != 0),
                Some(&1)
            );
            nodes
                .chunks_exact(4)
                .map(|node| {
                    let dot: u64 = plane.iter().zip(node).map(|(a, b)| a * b).sum();
                    dot.is_multiple_of(7)
                })
                .collect()
        })
        .collect();
    for (i, plane) in incidents.iter().enumerate() {
        assert_eq!(plane.iter().filter(|&&value| value).count(), 6);
        assert_eq!(incidents.iter().filter(|row| row[i]).count(), 6);
        for (node, &on_plane) in plane.iter().enumerate() {
            assert_eq!(on_plane, incidents[0][i ^ node]);
        }
        for other in &incidents[i + 1..] {
            assert_eq!(
                plane.iter().zip(other).filter(|&(a, b)| *a && *b).count(),
                2
            );
        }
    }
}

#[test]
fn the_central_factors_obey_the_cocycle_and_the_weil_commutator() {
    let mut source = marked_classes();
    for row in 0..4 {
        for col in 0..4 {
            let expected = if row == col { "o" } else { "z" };
            writeln!(
                source,
                "assert rvKummerThetaEntry[B,0,{row},{col}]=={expected};"
            )
            .unwrap();
        }
    }
    source.push_str("assert rvJacWeil2(e1,f1)==neg(o);\n");
    for a in 0..16 {
        writeln!(source, "assert rvKummerThetaFactor[B,0,{a}]==o;").unwrap();
        for b in 0..16 {
            writeln!(
                source,
                "assert rvKummerThetaFactor[B,{a},{b}]==rvUnramMul(rvKummerThetaFactor[B,{b},{a}],rvJacWeil2(n{a},n{b}));"
            )
            .unwrap();
        }
    }
    for (a, b, c) in [(1, 4, 8), (4, 1, 4), (2, 8, 5), (3, 12, 6), (15, 13, 7)] {
        let ab = a ^ b;
        let bc = b ^ c;
        writeln!(
            source,
            "assert rvUnramMul(rvKummerThetaFactor[B,{a},{b}],rvKummerThetaFactor[B,{ab},{c}])==rvUnramMul(rvKummerThetaFactor[B,{b},{c}],rvKummerThetaFactor[B,{a},{bc}]);"
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn extension_field_example_uses_a_nonprime_translation() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-theta")),
        concat!(
            "RVKummer(C)[[0, 0], [1, 0], [0, 1], [2, 0]]\n",
            "RVUnram(W)[2, 0]\n",
            "RVUnram(W)[0, 2]\n",
        )
    );
}

#[test]
fn theta_operations_check_masks_shapes_basis_and_nominal_type() {
    for (source, expected) in [
        ("show rvKummerTranslate[B,16](zero[RVKummer(C)]);", "mask must be between 0 and 15"),
        ("show rvKummerThetaEntry[B,1,4,0];", "row must be between 0 and 3"),
        ("show rvKummerThetaEntry[B,1,0,4];", "column must be between 0 and 3"),
        ("show rvKummerThetaFactor[B,1,16];", "mask must be between 0 and 15"),
        ("show rvKummerTropeCoeff[B,1,4];", "coefficient must be between 0 and 3"),
        ("show rvKummerThetaFactor[Unknown,1,0];", "unknown Jacobian torsion basis"),
        ("show rvKummerTranslate[B,1](zero[RVJac(C)]);", "expected type RVKummer(C)"),
        (
            "rvhyperelliptic Twin[F]{z,three,neg(o),z,four,o}; show rvKummerTranslate[B,1](zero[RVKummer(Twin)]);",
            "expected type RVKummer(C), found RVKummer(Twin)",
        ),
        (
            "rvhyperelliptic L[F]{z,o}; rvjactorsionbasis Empty[L,2](neg(o)){}; show rvKummerThetaFactor[Empty,0,0];",
            "certified genus-two level-2 basis",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F7}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
