mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const TYPE_III: &str = include_str!("../examples/twisted-double-lagrangian-boundaries.asl");
const KLEIN: &str = "rvgroup G{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
    rvdwcoeff P[G,2];";
const SEMION: &str = "rvgroup G[2];rvdwcoeff P[G,2];
    let k=rvPostnikovCochain[P,3]{1};rvdouble D[P](k);";
const TORSION: &str = "0,(1,2),(1,2),0,0,0,0,(1,2),(1,2)";

fn type_iii() -> &'static str {
    TYPE_III.split_once("// The span").unwrap().0
}

#[test]
fn type_iii_subgroup_and_alternating_cycle_match_the_oracle() {
    let output = stdout(example("twisted-double-lagrangian-boundaries"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines[0], "rvDoubleBoundaryData[Plane] = subgroup=[0, 1, 2, 3], obstructed=false, quantum_dimension=8, sectors=5, vacuum=1, modular_fixed=true");
    assert_eq!(
        lines[1],
        "RVDouble(D){[(0, 1), (1, 1), (8, 1), (10, 1), (12, 1)]}"
    );
    assert!(lines[2].ends_with("obstructed=true, witness_terms=6, pairing=1/2"));
    assert!(lines[3].starts_with("rvDoubleBoundaryWitness[Whole] = cycle=["));
    assert!(lines[3].ends_with("pairing=1/2"));
    assert_eq!(
        lines[4],
        "rvDoubleBoundaryCycle[D] = terms=6, pairing=1/2, obstructs=true"
    );
    assert_eq!(
        lines[5..],
        [
            "RVDouble(U){[(0, 1), (4, 1), (8, 1), (12, 1)]}",
            "RVDouble(U){[(0, 1), (5, 1), (10, 1), (15, 1)]}",
            "RVQZ(1/2)",
            "rvDoubleAnnulus[Base,Base] = 4",
            "rvDoubleAnnulus[Torsion,Torsion] = 4",
            "rvDoubleAnnulus[Base,Torsion] = 1",
        ]
    );
    let line = format!(
        "{}rvdoubleboundary Line[D]{{0,1}};
         show rvDoubleCondensate[Line];rvDoubleAnnulus[Line,Line];",
        type_iii()
    );
    assert_eq!(
        stdout(interpret(&line)),
        "RVDouble(D){[(0, 1), (1, 1), (2, 1), (3, 1), (8, 2)]}\nrvDoubleAnnulus[Line,Line] = 8\n"
    );
}

#[test]
fn the_same_bulk_vector_can_use_different_boundary_torsion() {
    let source = format!(
        "{}rvdoubleboundary A[D]{{3,2,0,1}};
         rvdoubleboundaryshift B[A]{{{TORSION}}};
         assert rvDoubleCondensate[A]==rvDoubleCondensate[B];
         assert rvDoubleBoundaryMu[A,1,2]==zero[RVQZ];
         assert rvDoubleBoundaryMu[B,1,2]==rvQZ[1,2];
         assert rvDoubleBoundaryMu[B,2,1]==zero[RVQZ];
         assert rvDoubleBoundaryMu[B,0,1]==zero[RVQZ];
         rvDoubleAnnulus[A,B];",
        type_iii()
    );
    assert_eq!(stdout(interpret(&source)), "rvDoubleAnnulus[A,B] = 5\n");
}

fn bockstein() -> String {
    let lift = |a: i32, b: i32| (a & 1) * ((b >> 1) & 1);
    let mut coeff = Vec::new();
    let mut known = Vec::new();
    for a in 1..4 {
        for b in 1..4 {
            known.push(format!("({},4)", lift(a, b)));
            for c in 1..4 {
                let delta = lift(b, c) - lift(a ^ b, c) + lift(a, b ^ c) - lift(a, b);
                assert_eq!(delta % 2, 0);
                coeff.push((delta / 2).rem_euclid(2).to_string());
            }
        }
    }
    format!(
        "{KLEIN}let k=rvPostnikovCochain[P,3]{{{}}};rvPostnikovClass(k);
         rvdouble D[P](k);rvdoubleboundary Auto[D]{{0,1,2,3}};
         rvDoubleBoundaryData[Auto];
         rvdoubleboundary Supplied[D]{{0,1,2,3}}{{{}}};
         rvDoubleBoundaryData[Supplied];",
        coeff.join(","),
        known.join(",")
    )
}

#[test]
fn integral_smith_lifts_a_nonzero_field_class_to_denominator_four() {
    let mut source = bockstein();
    for a in 1..4 {
        for b in 1..4 {
            writeln!(source, "show rvDoubleBoundaryMu[Auto,{a},{b}];").unwrap();
        }
    }
    source.push_str(
        "assert rvDoubleBoundaryMu[Supplied,1,2]==rvQZ[1,4];
         show rvDoubleCondensate[Auto];show rvDoubleCondensate[Supplied];",
    );
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].starts_with("rvPostnikovClass = nonzero=true"));
    for (line, name) in lines[1..3].iter().zip(["Auto", "Supplied"]) {
        assert_eq!(*line,format!(
            "rvDoubleBoundaryData[{name}] = subgroup=[0, 1, 2, 3], obstructed=false, quantum_dimension=4, sectors=4, vacuum=1, modular_fixed=true"
        ));
    }
    assert!(lines[3..12]
        .iter()
        .any(|x| matches!(*x, "RVQZ(1/4)" | "RVQZ(3/4)")));
    assert!(lines[12..]
        .iter()
        .all(|x| x.starts_with("RVDouble(D){[(0, 1),")));
}

#[test]
fn double_semion_and_odd_carry_have_genuine_u1_obstructions() {
    let source = format!(
        "{SEMION}
         rvdoubleboundary Full[D]{{0,1}};
         rvdoubleboundary Unit[D]{{0}}{{}};
         rvDoubleBoundaryWitness[Full];show rvDoubleCondensate[Unit];
         rvgroup C[3];rvdwcoeff P3[C,3];let k3=rvPostnikovCochain[P3,3]{{0,1,1,1,0,2,2,2}};
         rvdouble B[P3](k3);rvdoubleboundary All[B]{{0,1,2}};
         rvDoubleBoundaryData[All];
         rvDoubleBoundaryCycle[B]{{0,1,2}}{{(1,1,1,1),(1,2,1,1)}};"
    );
    let actual = stdout(interpret(&source));
    let lines = actual.lines().collect::<Vec<_>>();
    assert_eq!(
        lines[0],
        "rvDoubleBoundaryWitness[Full] = cycle=[([1, 1, 1], 1)], pairing=1/2"
    );
    assert_eq!(lines[1], "RVDouble(D){[(0, 1), (1, 1)]}");
    assert!(lines[2].contains("obstructed=true"));
    assert_eq!(
        lines[3],
        "rvDoubleBoundaryCycle[B] = terms=2, pairing=1/3, obstructs=true"
    );
}

#[test]
fn ordinary_qz_coefficients_can_have_order_outside_the_bulk_conductor() {
    let source = "
        rvgroup C[3];rvdwcoeff P[C,3];rvdouble D[P](rvPostnikovK[P]);
        rvdoubleboundary A[D]{0,1,2};
        rvdoubleboundary B[D]{0,1,2}{0,(3,17),(3,17),(3,17)};
        assert rvDoubleBoundaryMu[B,1,2]==rvQZ[3,17];
        assert rvDoubleCondensate[A]==rvDoubleCondensate[B];
        rvDoubleAnnulus[A,B];";
    assert_eq!(stdout(interpret(source)), "rvDoubleAnnulus[A,B] = 3\n");
}

#[test]
fn wrong_subgroups_cochains_cycles_and_nominal_annuli_are_checked() {
    for (suffix, expected) in [
        ("rvdoubleboundary A[D]{1};","must contain identity"),
        ("rvdoubleboundary A[D]{0,0};","distinct elements"),
        ("rvdoubleboundary A[D]{0,2};","element is out of range"),
        ("rvdoubleboundary A[D]{0,1}{(1,4)};","cochain fails dµ"),
        ("rvdoubleboundary A[D]{0,1}{};","requires 1 normalized cochain"),
        ("rvdoubleboundary A[D]{0,1};show rvDoubleCondensate[A];","has an obstructed subgroup"),
        ("rvdoubleboundary A[D]{0};rvDoubleBoundaryWitness[A];","has a trivializing cochain"),
        ("rvDoubleBoundaryCycle[D]{0,1}{(0,1,1,1)};","use nonidentity triples"),
        ("rvdoubleboundary A[D]{0};show rvDoubleBoundaryMu[A,0,1];","outside the subgroup"),
        ("rvdouble E[P](k);rvdoubleboundary A[D]{0};rvdoubleboundary B[E]{0};rvDoubleAnnulus[A,B];","same nominal bulk"),
    ] {
        let error = diagnostic(interpret(&format!("{SEMION}{suffix}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
    let source = format!(
        "{}rvDoubleBoundaryCycle[D]{{0,1,2,3,4,5,6,7}}{{(1,2,4,1)}};",
        type_iii()
    );
    assert!(diagnostic(interpret(&source)).contains("integral subgroup bar chain is not a cycle"));
    let source = "
        rvgroup G[3];rvdwcoeff P[G,3];rvdouble D[P](rvPostnikovK[P]);rvdoubleboundary A[D]{0,1};";
    assert!(diagnostic(interpret(source)).contains("not closed under multiplication"));
    let source = format!(
        "{KLEIN}rvdouble D[P](rvPostnikovK[P]);rvdoubleboundary A[D]{{0,1,2,3}};
         rvdoubleboundaryshift B[A]{{0,(1,2),0,0,0,0,0,0,0}};"
    );
    assert!(diagnostic(interpret(&source)).contains("cochain fails dµ"));
}
