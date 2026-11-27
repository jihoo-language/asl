mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const RP3_FIXTURE: &str = include_str!("../examples/finite-homotopy-path-integral.asl");
const LENS_FIXTURE: &str = include_str!("../examples/dijkgraaf-witten-lens.asl");

fn rp3_site() -> &'static str {
    RP3_FIXTURE.split_once("rvgroup H[4]").unwrap().0
}

fn rp3_data() -> (Vec<&'static str>, Vec<[&'static str; 2]>) {
    let mut points = Vec::new();
    let mut covers = Vec::new();
    let data = rp3_site()
        .split_once('{')
        .unwrap()
        .1
        .split_once('}')
        .unwrap()
        .0;
    for entry in data
        .split(',')
        .map(str::trim)
        .filter(|entry| !entry.is_empty())
    {
        if let Some((a, b)) = entry.split_once('<') {
            covers.push([a, b]);
        } else {
            points.push(entry);
        }
    }
    (points, covers)
}

fn signed_face_label(marked: usize, axes: &[usize]) -> String {
    let mut axes = axes.to_vec();
    axes.sort_unstable();
    let first = (marked >> axes[0]) & 1;
    let mut label = String::new();
    for i in axes {
        label.push(char::from(b'a' + i as u8));
        label.push(if (marked >> i) & 1 == first { 'p' } else { 'm' });
    }
    label
}

fn oriented_rp3_cycle(points: &[&str]) -> Vec<i32> {
    let mut flags = BTreeMap::new();
    for marked in (0usize..16).step_by(2) {
        for a in 0..4 {
            for b in (0..4).filter(|&b| b != a) {
                for c in (0..4).filter(|&c| c != a && c != b) {
                    let d = 6 - a - b - c;
                    let axes = [a, b, c, d];
                    let inversions = (0..4)
                        .map(|i| axes[i + 1..].iter().filter(|&&j| axes[i] > j).count())
                        .sum::<usize>();
                    let flag = [1, 2, 3, 4].map(|size| {
                        let name = signed_face_label(marked, &axes[..size]);
                        points.iter().position(|&p| p == name).unwrap()
                    });
                    // The outward facet orientation is the product of vertex signs.
                    let sign = if (inversions + marked.count_ones() as usize).is_multiple_of(2) {
                        1
                    } else {
                        -1
                    };
                    assert!(flags.insert(flag, sign).is_none());
                }
            }
        }
    }
    let mut boundary = BTreeMap::new();
    for (flag, &sign) in &flags {
        for i in 0..4 {
            let triangle = [&flag[..i], &flag[i + 1..]].concat();
            *boundary.entry(triangle).or_insert(0) +=
                sign * if i.is_multiple_of(2) { 1 } else { -1 };
        }
    }
    assert!(boundary.values().all(|&value| value == 0));
    assert_eq!(flags.len(), 192);
    flags.into_values().collect()
}

fn list<T: ToString>(values: impl IntoIterator<Item = T>) -> String {
    values
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn cycle() -> String {
    let (points, _) = rp3_data();
    list(oriented_rp3_cycle(&points))
}

#[test]
fn oriented_lens_fixture_has_nonreal_exact_value_and_conjugate() {
    let actual = stdout(example("dijkgraaf-witten-lens"));
    let lines = actual.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 7);
    assert_eq!(
        &lines[..3],
        [
            "RVCycloQ(3){[1, 2] / 3}",
            "RVCycloQ(3){[-1, -2] / 3}",
            "RVCycloQ(3){[1] / 1}",
        ]
    );
    assert!(lines[3].starts_with("rvDWComponent[Forward,0] = root=x000, vertices=120"));
    assert!(lines[3].ends_with(
        "generators=1, relators=1, normalized_solutions=3, conjugacy_classes=3, value=[1, 2]/3"
    ));
    for (line, phase) in lines[4..].iter().zip([0, 1, 1]) {
        assert!(line.ends_with(&format!("centralizer=3, orbit_size=1, phase={phase}")));
    }
}

#[test]
fn projective_space_cancellation_and_gauge_independence_use_an_integral_cycle() {
    let (points, covers) = rp3_data();
    let cycle = list(oriented_rp3_cycle(&points));
    let mut source = format!(
        "{}rvgroup H[4];rvgroup E[4];
        rvcrossed C[H,E]{{0,2,0,2}}{{{{0,1,2,3}},{{0,3,2,1}},{{0,1,2,3}},{{0,3,2,1}}}};
        rvpostnikov P[C,2];let k=rvPostnikovK[P];
        rvdw Twisted[RP3,P](k){{{cycle}}};show rvDWValue[Twisted];
        rvdw Untwisted[RP3,P](zero[RVPostnikovCochain(P,3)]){{{cycle}}};show rvDWValue[Untwisted];
        let lifted=rvPostnikovLifted[P]{{0,3}}{{{{0,0}},{{0,3}}}};
        rvdw Lifted[RP3,P](lifted){{{cycle}}};assert rvDWValue[Lifted]==rvDWValue[Twisted];",
        rp3_site()
    );
    for (name, change_gauge) in [("Cover", false), ("Changed", true)] {
        write!(source, "rvflat {name}[RP3,P]{{").unwrap();
        for [a, b] in &covers {
            let first = a.as_bytes()[0];
            let a_sign = a.as_bytes()[1];
            let b_sign = b.as_bytes()[b
                .as_bytes()
                .iter()
                .position(|&letter| letter == first)
                .unwrap()
                + 1];
            let label = usize::from(a_sign != b_sign);
            let label = if change_gauge {
                let i = points.iter().position(|&point| point == *a).unwrap();
                let j = points.iter().position(|&point| point == *b).unwrap();
                (label + i + j) % 2
            } else {
                label
            };
            write!(source, "{a}<{b}:{label},").unwrap();
        }
        writeln!(source, "}}; rvDWPhase[{name}](k){{{cycle}}};").unwrap();
    }
    let actual = stdout(interpret(&source));
    assert_eq!(
        actual,
        concat!(
            "RVCycloQ(2){[] / 1}\n",
            "RVCycloQ(2){[1] / 1}\n",
            "rvDWPhase[Cover] = 1\n",
            "rvDWPhase[Changed] = 1\n",
        )
    );
}

fn s3_mul(x: usize, y: usize) -> usize {
    let (a, b, c, d) = (x % 3, x / 3, y % 3, y / 3);
    (a + if b == 0 { c } else { 3 - c }) % 3 + 3 * ((b + d) % 2)
}

#[test]
fn nonabelian_twist_sums_conjugacy_classes_with_their_actual_centralizers() {
    let cycle = cycle();
    let mut source = format!("{}rvgroup G{{", rp3_site());
    for a in 0..6 {
        write!(source, "{{{}}},", list((0..6).map(|b| s3_mul(a, b)))).unwrap();
    }
    let pullback =
        list((1..6).flat_map(|a| {
            (1..6).flat_map(move |b| (1..6).map(move |c| (a / 3) * (b / 3) * (c / 3)))
        }));
    write!(
        source,
        "}};rvdwcoeff P[G,2];let k=rvPostnikovCochain[P,3]{{{pullback}}};
        rvdw A[RP3,P](k){{{cycle}}};show rvDWValue[A];rvDWComponent[A,0];
        rvDWClass[A,0,0];rvDWClass[A,0,1];
        rvdw U[RP3,P](zero[RVPostnikovCochain(P,3)]){{{cycle}}};show rvDWValue[U];",
    )
    .unwrap();
    let actual = stdout(interpret(&source));
    let lines = actual.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "RVCycloQ(2){[-1] / 3}");
    assert!(lines[1].ends_with("normalized_solutions=4, conjugacy_classes=2, value=[-1]/3"));
    assert!(lines[2].ends_with("centralizer=6, orbit_size=1, phase=0"));
    assert!(lines[3].ends_with("centralizer=2, orbit_size=3, phase=1"));
    assert_eq!(lines[4], "RVCycloQ(2){[2] / 3}");
}

#[test]
fn nonreal_disconnected_product_and_empty_space_use_exact_cyclotomic_arithmetic() {
    let end = LENS_FIXTURE.find("rvdw Reverse").unwrap();
    let first_integral = LENS_FIXTURE[..end].replacen("};", " extra\n};", 1);
    let source = format!(
        "{first_integral}show rvDWValue[Forward];rvDWComponent[Forward,1];
         space Empty{{}};rvdw U[Empty,P](k){{}};show rvDWValue[U];"
    );
    let actual = stdout(interpret(&source));
    let lines = actual.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "RVCycloQ(3){[1, 2] / 9}");
    assert!(lines[1].starts_with("rvDWComponent[Forward,1] = root=extra, vertices=1"));
    assert!(lines[1].ends_with("value=[1]/3"));
    assert_eq!(lines[2], "RVCycloQ(3){[1] / 1}");
}

#[test]
fn rejects_noncycles_noncocycles_twisted_actions_and_mistyped_inputs() {
    let base = "space X{x};rvgroup G[3];rvdwcoeff P[G,3];";
    let bad = [
        (
            "rvdw A[X,P](rvPostnikovCochain[P,3]{1,0,0,0,0,0,0,0}){};",
            "fails the twisted cocycle",
        ),
        (
            "rvdw A[X,P](zero[RVPostnikovCochain(P,3)]){1};",
            "requires 0 integral tetrahedron coefficients",
        ),
        (
            "rvdw A[X,P](zero[RVPostnikovCochain(P,2)]){};",
            "expected type RVPostnikovCochain(P,3)",
        ),
        (
            "rvgroup H[3];rvgroup E[2];rvcrossed C[H,E]{0,0,0}{{0,1,2},{0,2,1}};
          rvpostnikov Q[C,1];rvdw A[X,Q](rvPostnikovK[Q]){};",
            "require trivial coefficient action",
        ),
    ];
    for (suffix, expected) in bad {
        let error = diagnostic(interpret(&format!("{base}{suffix}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
    let nonintegral = list([1; 192]);
    let source = format!(
        "{}rvgroup G[2];rvdwcoeff P[G,2];rvdw A[RP3,P](rvPostnikovK[P]){{{nonintegral}}};",
        rp3_site()
    );
    // The all-positive chain is closed modulo two and fails the integral test.
    assert!(diagnostic(interpret(&source)).contains("integral chain is not a 3-cycle"));
}
