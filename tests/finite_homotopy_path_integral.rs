mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const FIXTURE: &str = include_str!("../examples/finite-homotopy-path-integral.asl");
const TWISTED: &str = "
    rvgroup H[3];rvgroup E[2];
    rvcrossed C[H,E]{0,0,0}{{0,1,2},{0,2,1}};
    rvpostnikov P[C,1];let k=rvPostnikovK[P];
";

fn rp3() -> &'static str {
    FIXTURE.split_once("rvgroup H[4]").unwrap().0
}

fn product(x: usize, y: usize) -> usize {
    let (a, b, c, d) = (x % 3, x / 3, y % 3, y / 3);
    (a + if b == 0 { c } else { 3 - c }) % 3 + 3 * ((b + d) % 2)
}

fn table<T>(rows: impl IntoIterator<Item = T>) -> String
where
    T: IntoIterator<Item = usize>,
{
    let mut out = "{".to_owned();
    for row in rows {
        out.push('{');
        for x in row {
            write!(&mut out, "{x},").unwrap();
        }
        out.push_str("},");
    }
    out.push('}');
    out
}

#[test]
fn rp3_partition_function_sums_both_normalized_backgrounds() {
    let actual = stdout(example("finite-homotopy-path-integral"));
    let lines: Vec<_> = actual.lines().collect();
    assert_eq!(lines.len(), 7);
    assert_eq!(
        lines[0],
        "rvHomotopyCardinality[Nonsplit] = 1, components=1"
    );
    assert!(lines[1].ends_with(
        "generators=1, relators=1, normalized_solutions=2, conjugacy_classes=2, cardinality=1"
    ));
    assert!(lines[2].starts_with(
        "rvHomotopyClass[Nonsplit,0,0] = labels=[], centralizer=2, orbit_size=1, betti=[1, 1, 1]"
    ));
    assert!(lines[2].ends_with("obstructed=false, fixed_cardinality=2, weight=1"));
    assert!(lines[3].contains("dp<cpdm:1"));
    assert!(lines[3].ends_with("obstructed=true, fixed_cardinality=0, weight=0"));
    assert_eq!(lines[4], "rvHomotopyCardinality[Split] = 2, components=1");
    assert_eq!(lines[5], "rvHomotopyCardinality[Loops] = 1, components=1");
    assert!(lines[6].ends_with(
        "generators=1, relators=0, normalized_solutions=2, conjugacy_classes=2, cardinality=1"
    ));
}

#[test]
fn nonabelian_pullback_crossed_module_counts_transpositions_in_one_orbit() {
    // Number the fiber product S3 ×_{C2} C4 by the S3 index and then n.
    let elements: Vec<_> = (0..6)
        .flat_map(|s| (0..4).filter(move |&n| n % 2 == s / 3).map(move |n| (s, n)))
        .collect();
    let e_table = table(elements.iter().map(|&(s, m)| {
        elements
            .iter()
            .map(|&(t, n)| {
                elements
                    .iter()
                    .position(|&x| x == (product(s, t), (m + n) % 4))
                    .unwrap()
            })
            .collect::<Vec<_>>()
    }));
    let action = table(
        elements
            .iter()
            .map(|&(_, n)| (0..4).map(move |h| if n % 2 == 0 { h } else { (4 - h) % 4 })),
    );
    let trivial = table(elements.iter().map(|_| 0..4));
    let source = format!(
        "{}rvgroup H[4];rvgroup Ext{e_table};
        rvcrossed C[H,Ext]{{0,1,0,1}}{action};
        rvcrossed C0[H,Ext]{{0,1,0,1}}{trivial};
        rvpostnikov P[C,2];rvpostnikov P0[C0,2];
        let k=rvPostnikovK[P];let k0=rvPostnikovK[P0];
        rvPostnikovCoeff[3,4,5](k);
        rvhomotopy I[RP3,P](k);rvHomotopyCardinality[I];rvHomotopyComponent[I,0];
        rvHomotopyClass[I,0,0];rvHomotopyClass[I,0,1];
        rvhomotopy J[RP3,P0](k0);rvHomotopyCardinality[J];rvHomotopyClass[J,0,1];",
        rp3()
    );
    let actual = stdout(interpret(&source));
    let lines: Vec<_> = actual.lines().collect();
    assert_eq!(lines[0], "rvPostnikovCoeff[3, 4, 5] = 1");
    assert_eq!(lines[1], "rvHomotopyCardinality[I] = 1/3, components=1");
    assert!(lines[2].ends_with(
        "generators=1, relators=1, normalized_solutions=4, conjugacy_classes=2, cardinality=1/3"
    ));
    assert!(lines[3].contains("labels=[], centralizer=6, orbit_size=1, betti=[1, 1, 1]"));
    assert!(lines[3].ends_with("obstructed=false, fixed_cardinality=2, weight=1/3"));
    assert!(lines[4].contains("centralizer=2, orbit_size=3, betti=[1, 1, 1]"));
    assert!(lines[4].ends_with("obstructed=true, fixed_cardinality=0, weight=0"));
    assert_eq!(lines[5], "rvHomotopyCardinality[J] = 4/3, components=1");
    assert!(lines[6].ends_with("obstructed=false, fixed_cardinality=2, weight=1"));
}

#[test]
fn torus_agrees_with_an_independent_count_of_commuting_pairs() {
    let e = table((0..6).map(|x| (0..6).map(move |y| product(x, y))));
    let trivial = table((0..6).map(|_| 0..2));
    let mut source = format!(
        "rvgroup H[2];rvgroup E{e};rvcrossed C[H,E]{{0,0}}{trivial};
        rvpostnikov P[C,1];let k=rvPostnikovK[P];space Torus{{"
    );
    for i in 0..4 {
        for j in 0..4 {
            write!(&mut source, "v{i}{j},").unwrap();
        }
    }
    for i in 0..4 {
        for x in 0..2 {
            for y in 2..4 {
                write!(&mut source, "v{i}{x}<v{i}{y},v{x}{i}<v{y}{i},").unwrap();
            }
        }
    }
    source
        .push_str("};rvhomotopy Q[Torus,P](k);rvHomotopyCardinality[Q];rvHomotopyComponent[Q,0];");
    for i in 0..8 {
        write!(&mut source, "rvHomotopyClass[Q,0,{i}];").unwrap();
    }
    let commuting_pairs = (0..6)
        .flat_map(|a| (0..6).map(move |b| (a, b)))
        .filter(|&(a, b)| product(a, b) == product(b, a))
        .count();
    let actual = stdout(interpret(&source));
    let lines: Vec<_> = actual.lines().collect();
    assert_eq!(
        lines[0],
        format!(
            "rvHomotopyCardinality[Q] = {}, components=1",
            commuting_pairs / 6
        )
    );
    assert!(lines[1].ends_with(&format!("generators=2, relators=1, normalized_solutions={commuting_pairs}, conjugacy_classes=8, cardinality=3")));
    assert_eq!(
        lines[2..]
            .iter()
            .filter(|line| line.contains("centralizer=2"))
            .count(),
        3
    );
    assert_eq!(
        lines[2..]
            .iter()
            .filter(|line| line.contains("centralizer=3"))
            .count(),
        4
    );
    assert_eq!(
        lines[2..]
            .iter()
            .filter(|line| line.contains("centralizer=6"))
            .count(),
        1
    );
    assert!(lines[2..]
        .iter()
        .all(|line| line.contains("betti=[1, 2, 1], obstructed=false, fixed_cardinality=1")));
}

#[test]
fn twisted_coefficients_and_gauge_changes_and_components_are_counted() {
    let source = format!(
        "{}{TWISTED}
        let b=rvPostnikovCochain[P,2]{{1}};
        let changed=rvPostnikovGauge(k,b);
        rvhomotopy Q[RP3,P](changed);
        rvHomotopyCardinality[Q];rvHomotopyClass[Q,0,0];rvHomotopyClass[Q,0,1];
        rvhomotopy Plain[RP3,P](k);rvHomotopyCardinality[Plain];
        space Empty{{}};rvhomotopy Void[Empty,P](k);rvHomotopyCardinality[Void];
        space Sum{{x,y,a,b,c,d,a<c,a<d,b<c,b<d}};
        rvhomotopy Disjoint[Sum,P](k);rvHomotopyCardinality[Disjoint];rvHomotopyComponent[Disjoint,2];",
        rp3()
    );
    let actual = stdout(interpret(&source));
    let lines: Vec<_> = actual.lines().collect();
    assert_eq!(lines[0], "rvHomotopyCardinality[Q] = 2, components=1");
    assert!(
        lines[1].ends_with("betti=[1, 0, 0], obstructed=false, fixed_cardinality=3, weight=3/2")
    );
    assert!(
        lines[2].ends_with("betti=[0, 0, 0], obstructed=false, fixed_cardinality=1, weight=1/2")
    );
    assert_eq!(lines[3], "rvHomotopyCardinality[Plain] = 2, components=1");
    assert_eq!(lines[4], "rvHomotopyCardinality[Void] = 1, components=0");
    assert_eq!(
        lines[5],
        "rvHomotopyCardinality[Disjoint] = 9/4, components=3"
    );
    assert!(lines[6].ends_with("normalized_solutions=2, conjugacy_classes=2, cardinality=1"));
}

#[test]
fn invalid_queries_and_assignment_limits_are_diagnosed() {
    let book = "space X{x};rvhomotopy Q[X,P](k);";
    for (code, expected) in [
        (
            "rvhomotopy Q[Missing,P](k);",
            "requires a named finite space",
        ),
        (
            "rvhomotopy Q[X,P](rvPostnikovCochain[P,2]{0});",
            "expected type RVPostnikovCochain(P,3)",
        ),
        (
            "rvHomotopyComponent[Q,1];",
            "component index is out of range",
        ),
        (
            "rvHomotopyClass[Q,0,1];",
            "conjugacy-class index is out of range",
        ),
        (
            "rvhomotopy Q[X,P](k);",
            "homotopy integral is already defined",
        ),
    ] {
        assert!(diagnostic(interpret(&format!("{TWISTED}{book}{code}"))).contains(expected));
    }
    let mut source = format!("{TWISTED}space X{{a,b,");
    for i in 0..20 {
        write!(&mut source, "a<x{i},b<x{i},").unwrap();
    }
    source.push_str("};rvhomotopy Q[X,P](k);");
    let error = diagnostic(interpret(&source));
    assert!(
        error.contains("normalized assignment limit after elimination"),
        "{error}"
    );
    let vertices = (0..90)
        .map(|i| format!("p{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let error = diagnostic(interpret(&format!(
        "{TWISTED}space X{{{vertices}}};rvhomotopy Q[X,P](k);"
    )));
    assert!(
        error.contains("exact signed 128-bit rational range"),
        "{error}"
    );
}
