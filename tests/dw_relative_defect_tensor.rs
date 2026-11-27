mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-relative-tensor.asl");
const FOUR: &str = "sectors=[(0, 1), (1, 1), (2, 1), (3, 1)]";

fn type_iii_inputs() -> &'static str {
    EXAMPLE.split_once("let product=").unwrap().0
}

#[test]
fn nontrivial_mixed_associator_and_cocycle_phase_are_on_actual_balanced_bases() {
    let result = stdout(example("twisted-double-relative-tensor"));
    assert!(result.contains(
        "dimension=16, conductor=4, grades=[(4, 4), (5, 4), (6, 4), (7, 4)], sectors=[(4, 2)]"
    ));
    assert!(result.contains("rvDoubleRelativeHom[4] = multiplicity=2, rows=4, columns=2"));
    assert!(result.contains("rvDoubleRelativeBasis[4,0] = grades=[0, 4], columns=[0, 0]"));
    assert!(result
        .contains("rvDoubleRelativeAction[0,1,4,0] = grade=5, target=1, phase=0, conductor=4"));
    assert!(result.contains(
        "rvDoubleRelativeInducedAction[0,1,0,0] = grade=1, target=0, phase=0, conductor=4"
    ));
    assert!(result
        .contains("dimension=16, conductor=4, nonzero_phases=4, sectors=[(0, 4)], verified=true"));
    assert!(result.contains("rvDoubleRelativeAssocImage[0,3] = target=3, phase=2, conductor=4"));
    assert_eq!(result.matches("dimension=16, conductor=4, grades=[(0, 4), (1, 4), (2, 4), (3, 4)], sectors=[(0, 4)]").count(), 2);
    assert!(result.contains("rvDoubleRelativeAssocImage[7,0] = target=0, phase=2, conductor=4"));
}

#[test]
fn multiplicity_two_has_independent_inclusions_and_exact_hom_projector() {
    let mut source = format!(
        "{}let product=rvDoubleRelative[PP](y0,x4);",
        type_iii_inputs()
    );
    for b in 0..2 {
        for r in 0..4 {
            for s in 0..2 {
                writeln!(source, "rvDoubleRelativeInclusion[4,{b},{r},{s}](product);").unwrap();
            }
        }
    }
    for r in 0..8 {
        for s in 0..8 {
            writeln!(source, "rvDoubleRelativeProjector[4,{r},{s}](product);").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let scalars: Vec<_> = output
        .lines()
        .map(|line| match line.split_once(" = ").unwrap().1 {
            "RVCycloQ(4){[] / 1}" => (0, 1),
            "RVCycloQ(4){[1] / 1}" => (1, 1),
            "RVCycloQ(4){[-1] / 1}" => (-1, 1),
            "RVCycloQ(4){[1] / 4}" => (1, 4),
            "RVCycloQ(4){[-1] / 4}" => (-1, 4),
            other => panic!("unexpected coefficient {other}"),
        })
        .collect();
    let inclusions = [&scalars[..8], &scalars[8..16]];
    for (b, u) in inclusions.iter().enumerate() {
        for (c, v) in inclusions.iter().enumerate() {
            for i in 0..2 {
                for j in 0..2 {
                    let inner = (0..4).map(|r| u[2 * r + i].0 * v[2 * r + j].0).sum::<i32>();
                    assert_eq!(inner, 2 * i32::from(b == c && i == j));
                }
            }
        }
    }
    for r in 0..8 {
        for s in 0..8 {
            let (n, d) = scalars[16 + 8 * r + s];
            let expected = inclusions.iter().map(|u| u[r].0 * u[s].0).sum::<i32>();
            assert_eq!(4 * n, expected * d);
        }
    }
}

#[test]
fn mixed_sector_products_match_independent_counts_and_same_boundary_fusion() {
    let source = format!(
        "{}let y4=rvDoubleDefectBasis[Y,1];
        let x0=rvDoubleDefectBasis[X,0];
        let m=rvDoubleDefectBasis[PP,4];
        assert rvDoubleDefectOpposite[Y](x4)==y4;
        assert rvDoubleDefectOpposite[X](y0)==x0;
        rvDoubleRelativeData(rvDoubleRelative[PP](y0,x0));
        rvDoubleRelativeData(rvDoubleRelative[PP](y4,x4));
        rvDoubleRelativeData(rvDoubleRelative[PP](y4,x0));
        rvDoubleRelativeData(rvDoubleRelative[PP](m,m));
        rvDoubleDefectFuse(m,m);
        rvDoubleRelativeData(rvDoubleRelative[Y](m,y4));
        rvDoubleRelativeData(rvDoubleRelative[QQ](x4,y4));
        let product=rvDoubleRelative[PP](y0,x4);
        assert product==rvDoubleRelative[PP](y0,x4);",
        type_iii_inputs()
    );
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    for i in [0, 1, 3, 4, 6] {
        assert!(lines[i].ends_with(FOUR), "{}", lines[i]);
    }
    assert!(lines[2].ends_with("sectors=[(4, 2)]"));
    assert!(lines[5].ends_with("sectors=[(0, 2)]"));
    assert!(lines[6].contains("grades=[(0, 1), (1, 1), (2, 1), (3, 1)]"));
}

#[test]
fn cyclic_point_and_full_boundaries_have_all_grades_and_all_characters() {
    for p in [2, 3] {
        let members = (0..p).map(|i| i.to_string()).collect::<Vec<_>>().join(",");
        let source = format!(
            "rvgroup G[{p}];rvdwcoeff P[G,{p}];rvdouble D[P](rvPostnikovK[P]);
            rvdoubleboundary Point[D]{{0}};rvdoubleboundary Full[D]{{{members}}};
            rvdoubledefects A[Point,Full];rvdoubledefects B[Full,Point];
            rvdoubledefects C[Point,Point];rvdoubledefects E[Full,Full];
            let x=rvDoubleDefectBasis[A,0];let y=rvDoubleDefectBasis[B,0];
            assert rvDoubleDefectOpposite[B](x)==y;
            rvDoubleRelativeData(rvDoubleRelative[C](x,y));
            rvDoubleRelativeData(rvDoubleRelative[E](y,x));"
        );
        let output = stdout(interpret(&source));
        let lines: Vec<_> = output.lines().collect();
        let sectors: Vec<_> = (0..p).map(|i| (i, 1)).collect();
        let full: Vec<_> = (0..p).map(|i| (i, p)).collect();
        assert_eq!(lines[0], format!("rvDoubleRelativeData = dimension={p}, conductor={}, grades={sectors:?}, sectors={sectors:?}", p*p));
        assert_eq!(lines[1], format!("rvDoubleRelativeData = dimension={}, conductor={}, grades={full:?}, sectors={sectors:?}", p*p, p*p));
    }
}

#[test]
fn odd_carry_cocycle_uses_the_complex_inverse_associator_phase() {
    let source = "rvgroup G[3];rvdwcoeff P[G,3];
        let k=rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2};rvdouble D[P](k);
        rvdoubleboundary Point[D]{0};rvdoubledefects C[Point,Point];
        let x=rvDoubleDefectBasis[C,1];let y=rvDoubleDefectBasis[C,2];
        rvDoubleRelativeAssocImage[C,1,0](x,x,y);";
    assert_eq!(
        stdout(interpret(source)),
        "rvDoubleRelativeAssocImage[1,0] = target=0, phase=6, conductor=9\n"
    );
}

#[test]
fn products_distinguish_nominal_middle_and_outer_boundaries() {
    let setup = format!(
        "{}rvdoubleboundary Copy[D]{{0}};
        rvdoubledefects XC[Copy,Plane];let wrong=rvDoubleDefectBasis[XC,1];",
        type_iii_inputs()
    );
    for (operation, expected) in [
        (
            "show rvDoubleRelative[PP](y0,wrong);",
            "identical middle boundary",
        ),
        (
            "show rvDoubleRelative[QQ](y0,x4);",
            "ordered outer boundary",
        ),
        (
            "rvDoubleRelativeData(y0);",
            "requires a graded defect bimodule",
        ),
        (
            "rvDoubleRelativeBasis[0,0](rvDoubleRelative[PP](y0,x4));",
            "grade or basis column is out of range",
        ),
        (
            "rvDoubleRelativeInclusion[4,2,0,0](rvDoubleRelative[PP](y0,x4));",
            "inclusion index is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{operation}")));
        assert!(error.contains(expected), "{expected}: {error}");
    }
}

#[test]
fn quarter_trivializations_survive_balancing_and_reversed_boundaries() {
    let lift = |a: i32, b: i32| (a & 1) * ((b >> 1) & 1);
    let k = (1..4)
        .flat_map(|a| {
            (1..4).flat_map(move |b| {
                (1..4).map(move |c| {
                    ((lift(b, c) - lift(a ^ b, c) + lift(a, b ^ c) - lift(a, b)) / 2)
                        .rem_euclid(2)
                        .to_string()
                })
            })
        })
        .collect::<Vec<_>>()
        .join(",");
    let mu = (1..4)
        .flat_map(|a| (1..4).map(move |b| format!("({},4)", lift(a, b))))
        .collect::<Vec<_>>()
        .join(",");
    let source = format!(
        "rvgroup G{{{{0,1,2,3}},{{1,0,3,2}},{{2,3,0,1}},{{3,2,1,0}}}};
        rvdwcoeff P[G,2];rvdouble D[P](rvPostnikovCochain[P,3]{{{k}}});
        rvdoubleboundary A[D]{{0,1,2,3}};
        rvdoubleboundary B[D]{{0,1,2,3}}{{{mu}}};
        rvdoubledefects AB[A,B];rvdoubledefects BA[B,A];
        rvdoubledefects AA[A,A];rvdoubledefects BB[B,B];
        let x=rvDoubleDefectBasis[AB,0];let y=rvDoubleDefectBasis[BA,0];
        assert rvDoubleDefectOpposite[BA](x)==y;
        let product=rvDoubleRelative[AA](x,y);
        rvDoubleRelativeData(product);
        rvDoubleRelativeData(rvDoubleRelative[BB](y,x));
        let unit=rvDoubleDefectBasis[AA,0];
        rvDoubleRelativeData(rvDoubleRelative[AA](product,unit));
        rvDoubleRelativeAssociator[AB](x,y,x);
        rvDoubleRelativeAssociator[AA](product,unit,unit);"
    );
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    for line in &lines[..3] {
        assert_eq!(*line, format!("rvDoubleRelativeData = dimension=16, conductor=8, grades=[(0, 4), (1, 4), (2, 4), (3, 4)], {FOUR}"));
    }
    assert!(lines[3].contains("dimension=32, conductor=8"));
    assert!(lines[3].ends_with("sectors=[(0, 4)], verified=true"));
    assert!(lines[4].contains("dimension=16, conductor=8"));
    assert!(lines[4].ends_with(&format!("{FOUR}, verified=true")));
}
