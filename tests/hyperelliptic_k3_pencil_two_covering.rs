mod common;

use std::collections::BTreeSet;

use common::{diagnostic, example, interpret, stdout};

const F23: &str = include_str!("../examples/hyperelliptic-k3-pencil-two-covering.asl");
const POINT: &str =
    "rvK3PencilSectionPoint[P]{[0],[18],[2],[5]}{[1]}{[0],[0],[18],[11],[10],[5],[6]}{[1]}";
const SECOND: &str =
    "rvK3PencilSectionPoint[P]{[0],[18],[4],[5]}{[1]}{[0],[0],[14],[7],[18],[5],[6]}{[1]}";

fn setup() -> &'static str {
    F23.split_once("\nlet p=").unwrap().0
}

fn span(rows: &[u16]) -> BTreeSet<u16> {
    (0..1 << rows.len())
        .map(|mask| {
            rows.iter()
                .enumerate()
                .filter(|&(i, _)| mask & (1 << i) != 0)
                .fold(0, |v, (_, row)| v ^ row)
        })
        .collect()
}

fn permute_places(vector: u16) -> u16 {
    // The independent computation uses t-values in increasing order. The
    // interpreter sorts by ascending coefficients of monic finite factors.
    let indices = [0, 1, 7, 6, 5, 4, 3, 2];
    [0, 8].iter().fold(0, |mut output, offset| {
        for (old, new) in indices.iter().enumerate() {
            output |= ((vector >> (old + offset)) & 1) << (new + offset);
        }
        output
    })
}

#[test]
fn actual_quadrics_and_all_three_hasse_failures_match_independent_mod_two_span() {
    let output = stdout(example("hyperelliptic-k3-pencil-two-covering"));
    let mw_span = span(&[0x1e0d, 0x6d1f, 0x0112, 0x490a].map(permute_places));
    let mw_report = output.lines().next().unwrap();
    assert!(mw_report.contains("Mordell-Weil rank=2, image dimension=4, dim Sha[2]=2"));
    let reported_rows: Vec<u16> = mw_report
        .split_once("echelon masks=[")
        .unwrap()
        .1
        .trim_end_matches(']')
        .split(", ")
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(mw_span.len(), 16);
    assert_eq!(span(&reported_rows), mw_span);

    let masks = [0x0a0a, 0x0615, 0x0a0a ^ 0x0615].map(permute_places);
    let expected = masks.map(|mask| {
        mw_span
            .iter()
            .map(|row| row ^ mask)
            .min_by_key(|v| v.reverse_bits())
            .unwrap()
    });
    assert_eq!(expected.into_iter().collect::<BTreeSet<_>>().len(), 3);
    assert!(expected.iter().all(|&class| class != 0));
    for (name, coset) in ["A", "Bcover", "Sum"].into_iter().zip(expected) {
        assert!(output.contains(&format!("rvK3PencilCoverStatus[{name}, M] = everywhere locally soluble, outside certified global image, rational point=no, Hasse-principle failure, Sha coset mask={coset}")));
    }
    assert!(output.contains("in certified global image, rational point=yes, Sha coset mask=0"));
    assert!(output.contains("rvK3PencilCoverQuadrics[A] = diagonal (U0²,U1²,U2²,W²)=[[[[0], [13], [1]], [[0], [10], [22]], [], [[0], [17], [4], [22], [6]]], [[[0], [13], [1]], [], [[22]], [[0], [14], [10], [11], [22]]]], product square root=[[0], [13], [1]]"));
}

#[test]
fn literal_points_projectivize_and_cover_map_and_lifts_certify_actual_sections() {
    let source = format!(
        r"{}
        let p={POINT};
        let origin=zero[RVK3PencilSection(P)];
        let trivial=rvK3PencilKummer2[P]{{[4],[8],[4]}}{{[9]}}{{[1]}};
        assert trivial==zero[RVK3PencilKummer2(P)];
        rvk3pencilcover T[P](trivial);
        let infinity=rvK3PencilCoverPoint[T]{{[1]}}{{[1]}}{{[1]}}{{}};
        let infinityScaled=rvK3PencilCoverPoint[T]{{[0],[3]}}{{[0],[3]}}{{[0],[3]}}{{}};
        assert infinity==infinityScaled;
        assert rvK3PencilCoverMap(infinity)==origin;
        assert rvK3PencilCoverLift[T](origin)==infinity;
        let twice=add(p,p);
        let denominatorLift=rvK3PencilCoverLift[T](twice);
        assert rvK3PencilCoverMap(denominatorLift)==twice;
        let b=rvK3PencilDescent2(p);
        rvk3pencilcover H[P](b);
        let lifted=rvK3PencilCoverLift[H](p);
        let literal=rvK3PencilCoverPoint[H]{{[2],[1]}}{{[0],[6],[18]}}{{[15],[19]}}{{[6]}};
        let scaled=rvK3PencilCoverPoint[H]{{[10],[11],[3]}}{{[0],[7],[16],[8]}}{{[6],[2],[11]}}{{[7],[18]}};
        assert literal==scaled;
        assert literal==lifted;
        assert rvK3PencilCoverMap(scaled)==p;
        assert rvK3PencilCoverMap(rvK3PencilCoverLift[H](neg(p)))==neg(p);
        let t0=rvK3PencilTorsion2At[P,0];
        rvk3pencilcover R[P](rvK3PencilDescent2(t0));
        assert rvK3PencilCoverMap(rvK3PencilCoverLift[R](t0))==t0;
        show lifted;
        ",
        setup(),
    );
    let output = stdout(interpret(&source));
    assert!(output.starts_with("RVK3PencilCoverPoint(H)["));
}

#[test]
fn local_obstructions_are_distinct_from_uncertified_global_status() {
    let source = format!(
        r"{}
        let p={POINT};
        let q={SECOND};
        rvk3pencilmw2 M[P]{{p,q}};
        rvk3pencilmw2 Reversed[P]{{q,p}};
        let a=rvK3PencilKummer2[P]{{[0],[-10],[1]}}{{[0],[-10],[1]}}{{[1]}};
        rvk3pencilcover A[P](a);
        rvK3PencilCoverStatus[A];
        assert rvK3PencilShaClass[A,M]==rvK3PencilShaClass[A,Reversed];
        rvk3pencilcover Outside[P](rvK3PencilKummer2[P]{{[-1],[1]}}{{[-1],[1]}}{{[1]}});
        rvk3pencilcover Inside[P](rvK3PencilKummer2[P]{{[5]}}{{[5]}}{{[1]}});
        rvK3PencilCoverStatus[Outside,M];
        rvK3PencilCoverStatus[Inside,M];
        rvK3PencilCoverStatus[Inside];
        ",
        setup(),
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvK3PencilCoverStatus[A] = everywhere locally soluble, global image uncertified\n\
         rvK3PencilCoverStatus[Outside, M] = local obstruction, rational point=no\n\
         rvK3PencilCoverStatus[Inside, M] = local obstruction, rational point=no\n\
         rvK3PencilCoverStatus[Inside] = local obstruction, rational point=no\n"
    );
}

#[test]
fn torsion_lifts_and_local_solubility_work_over_a_characteristic_three_extension() {
    let output = stdout(interpret(
        r"
        rvunram F[3,1]{1,0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let i=rvUnramElem[F]{0,1};
        rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
        let r0=rvJacPoint[C](z,z);
        let r1=rvJacPoint[C](o,z);
        let r2=rvJacPoint[C](neg(o),z);
        let r3=rvJacPoint[C](i,z);
        rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
        rvk3pencil P[B,0,1];
        let t=rvK3PencilTorsion2At[P,1];
        rvk3pencilcover T[P](rvK3PencilDescent2(t));
        assert rvK3PencilCoverMap(rvK3PencilCoverLift[T](t))==t;
        rvK3PencilCoverStatus[T];
        ",
    ));
    assert_eq!(
        output,
        "rvK3PencilCoverStatus[T] = everywhere locally soluble, global image uncertified\n"
    );
}

#[test]
fn invalid_points_classes_certificates_and_nominal_mixes_report_errors() {
    let source = format!(
        "{} let p={POINT}; let q={SECOND}; rvk3pencilmw2 M[P]{{p,q}};\n\
         let trivial=zero[RVK3PencilKummer2(P)];\n\
         rvk3pencilcover T[P](trivial); rvk3pencilcover U[P](trivial);\n\
         let A=rvK3PencilCoverPoint[T]{{[1]}}{{[1]}}{{[1]}}{{}};\n\
         let Bpoint=rvK3PencilCoverPoint[U]{{[1]}}{{[1]}}{{[1]}}{{}};\n\
         rvk3pencilcover N[P](rvK3PencilKummer2[P]{{[5]}}{{[5]}}{{[1]}});\n\
         rvk3pencil P2[B,14,0];",
        setup(),
    );
    for (operation, expected) in [
        (
            "show rvK3PencilKummer2[P]{}{[1]}{[1]};",
            "representative must be nonzero",
        ),
        (
            "show rvK3PencilKummer2[P]{[5]}{[1]}{[1]};",
            "must have square product",
        ),
        (
            "show rvK3PencilCoverPoint[T]{}{}{}{};",
            "projective point must be nonzero",
        ),
        (
            "show rvK3PencilCoverPoint[T]{[1]}{[1]}{}{};",
            "does not lie on both",
        ),
        (
            "show rvK3PencilCoverPoint[N]{[1]}{[1]}{[1]}{};",
            "does not lie on both",
        ),
        (
            "show rvK3PencilCoverLift[T](p);",
            "section is outside the covering image",
        ),
        ("show rvK3PencilShaClass[N,M];", "has a local obstruction"),
        (
            "rvk3pencilmw2 M2[P]{p,p};",
            "image has dimension 3, expected 4",
        ),
        (
            "rvk3pencilcover T2[P2](trivial);",
            "expected type RVK3PencilKummer2(P2)",
        ),
        ("assert A==Bpoint;", "expected type RVK3PencilCoverPoint(T)"),
        (
            "show rvK3PencilCoverLift[T](zero[RVK3PencilSection(P2)]);",
            "expected type RVK3PencilSection(P)",
        ),
        (
            "rvk3pencilmw2 M2[P2]{p,q};",
            "expected type RVK3PencilSection(P2)",
        ),
        (
            "rvk3pencilcover Z[P2](zero[RVK3PencilKummer2(P2)]); rvK3PencilCoverStatus[Z,M];",
            "require the same named pencil",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{source}\n{operation}")));
        assert!(error.contains(expected), "{operation}: {error}");
    }
}
