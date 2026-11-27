mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-bulk-central-functor.asl");

fn type_iii() -> &'static str {
    EXAMPLE.split_once("// Induce").unwrap().0
}

#[test]
fn example_checks_the_bulk_fusion_table_and_non_scalar_central_move() {
    let output = stdout(example("twisted-double-bulk-central-functor"));
    assert!(output.contains("grades=[(4, 2), (5, 2), (6, 2), (7, 2)], sectors=[(4, 1)]"));
    assert!(output.contains(
        "bulk_simples=22, boundary_simples=5, vacuum=0, fusion_pairs=484, verified=true"
    ));
    assert!(output.contains(
        "bulk_simples=22, boundary_simples=8, vacuum=0, fusion_pairs=484, verified=true"
    ));
    assert!(output.contains("rvDoubleBulkBoundaryFusion[14,14] = sectors=[(0, 1), (1, 1), (2, 1), (3, 1)], verified=true"));
    assert!(output.contains("rvDoubleCentralImage[6,0] = target=1, phase=0, conductor=4"));
    assert!(output.contains("rvDoubleCentralBasis[0,6,0] = grades=[4, 2], columns=[0, 0]"));
    assert!(output.contains("rvDoubleCentralUnit[14] = columns=8, conductor=4, verified=true"));
    assert!(output.contains("rvDoubleCentralWhisker = dimension=32, conductor=4, nonzero_phases=0, moved_columns=24, verified=true"));
    assert!(output.contains("rvDoubleCentralImage[2,0] = target=0, phase=8, conductor=9"));
}

#[test]
fn induced_bimodules_branch_at_the_checked_representative_and_compose() {
    let source = format!(
        "{}rvDoubleRelativeData(rvDoubleBulkBoundary[PP,4]);
         rvDoubleRelativeData(rvDoubleBulkBoundary[PP,8]);
         rvDoubleRelativeData(rvDoubleBulkBoundary[PP,9]);
         rvDoubleRelativeData(rvDoubleBulkBoundary[PP,12]);
         rvDoubleRelativeData(rvDoubleBulkBoundary[QQ,8]);
         let b=rvDoubleBulkBoundary[PP,14];
         rvDoubleRelativeBasis[5,0](b);
         rvDoubleRelativeAction[0,2,4,0](b);
         rvDoubleRelativeAction[1,2,4,0](b);
         rvDoubleRelativeData(rvDoubleRelative[PP](b,rvDoubleBulkBoundary[PP,0]));",
        type_iii()
    );
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    for (line, sectors) in lines[..5].iter().zip([
        "[(2, 1)]",
        "[(0, 1), (1, 1)]",
        "[(2, 1), (3, 1)]",
        "[(0, 1), (3, 1)]",
        "[(1, 2)]",
    ]) {
        assert!(line.ends_with(sectors), "{line}");
    }
    assert_eq!(
        lines[5],
        "rvDoubleRelativeBasis[5,0] = grades=[4, 1], columns=[0, 0]"
    );
    assert_eq!(
        lines[6],
        "rvDoubleRelativeAction[0,2,4,0] = grade=6, target=1, phase=0, conductor=4"
    );
    assert_eq!(
        lines[7],
        "rvDoubleRelativeAction[1,2,4,0] = grade=6, target=0, phase=0, conductor=4"
    );
    assert!(lines[8].ends_with("sectors=[(4, 1)]"));
}

#[test]
fn exact_point_moves_use_the_inverse_of_each_independent_bulk_operator() {
    let mut program = type_iii().to_owned();
    for x in 0..8 {
        let g = x ^ 4;
        writeln!(
            program,
            "let m{x}=rvDoubleCentralMove[QQ,QQ,14](rvDoubleDefectBasis[QQ,{x}]);"
        )
        .unwrap();
        for i in 0..2 {
            for j in 0..2 {
                writeln!(program,
                    "assert rvDoubleCentralEntry[{g},{i},{j}](m{x})==rvCycloQGalois[-1](rvDoubleRepEntry[D,14,{x},{j},{i}]);
                     assert rvDoubleCentralInvEntry[{g},{i},{j}](m{x})==rvDoubleRepEntry[D,14,{x},{i},{j}];"
                ).unwrap();
            }
        }
    }
    assert_eq!(stdout(interpret(&program)), "");
    let odd = "rvgroup G[3];rvdwcoeff P[G,3];
        rvdouble D[P](rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2});
        rvdoubleboundary Z[D]{0};rvdoubledefects E[Z,Z];
        let a=rvDoubleDefectBasis[E,1];
        let f=rvDoubleCentralMove[E,E,3](a);
        assert rvDoubleCentralEntry[2,0,0](f)==rvCycloRat[1](rvCyclo[9,8]);
        assert rvDoubleCentralInvEntry[2,0,0](f)==rvCycloRat[1](rvCyclo[9,1]);
        rvDoubleCentralWhisker[E,E,E,E,3](a,a);";
    assert_eq!(stdout(interpret(odd)),
        "rvDoubleCentralWhisker = dimension=1, conductor=9, nonzero_phases=1, moved_columns=0, verified=true\n");
}

#[test]
fn whiskering_is_coherent_across_three_distinct_boundaries_and_graded_products() {
    let source = format!(
        "{}rvdoubleboundary Line[D]{{0,1}};
         rvdoubledefects LL[Line,Line];rvdoubledefects PL[Plane,Line];rvdoubledefects LQ[Line,Point];
         let p=rvDoubleDefectBasis[PL,0];let q=rvDoubleDefectBasis[LQ,2];
         rvDoubleCentralWhisker[PP,LL,QQ,Y,14](p,q);
         let product=rvDoubleRelative[Y](p,q);
         let c=rvDoubleCentralMove[PP,QQ,14,4](product);
         rvDoubleCentralData(c);
         rvDoubleCentralWhisker[PP,QQ,PP,PP,14](product,x4);",
        type_iii()
    );
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines[0], "rvDoubleCentralWhisker = dimension=8, conductor=4, nonzero_phases=2, moved_columns=4, verified=true");
    assert!(lines[1].contains("defect=graded, dimension=8, conductor=4"));
    assert!(lines[2].contains("dimension=32, conductor=4"));
    assert!(lines[2].ends_with("verified=true"));
}

fn quarter_boundaries() -> String {
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
    format!("rvgroup G{{{{0,1,2,3}},{{1,0,3,2}},{{2,3,0,1}},{{3,2,1,0}}}};
         rvdwcoeff P[G,2];rvdouble D[P](rvPostnikovCochain[P,3]{{{k}}});
         rvdoubleboundary A[D]{{0,1,2,3}};rvdoubleboundary B[D]{{0,1,2,3}}{{{mu}}};
         rvdoubleboundary Z[D]{{0}};
         rvdoubledefects AA[A,A];rvdoubledefects BB[B,B];rvdoubledefects ZZ[Z,Z];
         rvdoubledefects AB[A,B];rvdoubledefects BA[B,A];rvdoubledefects BZ[B,Z];rvdoubledefects AZ[A,Z];
         let a=rvDoubleDefectBasis[AB,0];let b=rvDoubleDefectBasis[BZ,0];
         let product=rvDoubleRelative[AA](a,rvDoubleDefectBasis[BA,0]);")
}

#[test]
fn quarter_boundary_phases_survive_balancing_and_dynamic_conductor_checks() {
    let setup = quarter_boundaries();
    let source = format!(
        "{setup}
        rvDoubleBulkBoundaryValidate[AA]; rvDoubleBulkBoundaryValidate[BB];
        rvDoubleCentralData(rvDoubleCentralMove[AA,BB,4](a));
        rvDoubleCentralWhisker[AA,BB,ZZ,AZ,4](a,b);
        rvDoubleCentralData(rvDoubleCentralMove[AA,AA,4,8](product));
        rvDoubleCentralUnit[BB,5];"
    );
    let text = stdout(interpret(&source));
    assert_eq!(
        text.matches("bulk_simples=16, boundary_simples=4, vacuum=0, fusion_pairs=256")
            .count(),
        2
    );
    assert!(text.contains("bulk=4, defect=0, dimension=8, conductor=8"));
    assert!(text.contains("rvDoubleCentralWhisker = dimension=8, conductor=8, nonzero_phases=8, moved_columns=0, verified=true"));
    assert!(text.contains("bulk=4, defect=graded, dimension=16, conductor=8"));
    assert!(text.contains("rvDoubleCentralUnit[5] = columns=4, conductor=4, verified=true"));
    for (command, error) in [
        (
            "show rvDoubleCentralMove[AA,AA,4](product);",
            "requires an explicit conductor",
        ),
        (
            "show rvDoubleCentralMove[AA,AA,4,4](product);",
            "does not contain the graded operand's phases",
        ),
    ] {
        assert!(diagnostic(interpret(&format!("{setup}{command}"))).contains(error));
    }
}

#[test]
fn central_maps_reject_wrong_nominal_boundaries_indices_and_field_choices() {
    let setup = format!(
        "{}rvdoubleboundary Copy[D]{{0}};rvdoubledefects Other[Copy,Copy];
         let f=rvDoubleCentralMove[PP,QQ,14](y0);",
        type_iii()
    );
    assert!(stdout(interpret(&format!("{setup}show f;"))).starts_with("RVDoubleCentralMap(Y,4)"));
    for (command, expected) in [
        ("show rvDoubleBulkBoundary[Y,0];", "nominal self-boundary"),
        (
            "show rvDoubleBulkBoundary[PP,100];",
            "simple index is out of range",
        ),
        (
            "show rvDoubleCentralMove[PP,Other,14](y0);",
            "exact ordered defect boundaries",
        ),
        (
            "show rvDoubleCentralMove[PP,QQ,14,3](y0);",
            "positive multiple of the boundary conductor",
        ),
        (
            "rvDoubleCentralImage[0,0](f);",
            "grade or column is out of range",
        ),
        (
            "show rvDoubleCentralEntry[6,0,0](y0);",
            "requires a bulk central map",
        ),
        (
            "rvDoubleCentralWhisker[PP,QQ,PP,QQ,14](y0,x4);",
            "ordered outer boundary",
        ),
    ] {
        let err = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(err.contains(expected), "{expected} in {err}");
    }
}
