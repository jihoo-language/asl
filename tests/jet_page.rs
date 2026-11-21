mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn tangent_and_transverse_arcs_have_distinct_page_differentials() {
    assert_eq!(
        stdout(example("jet-pages")),
        concat!(
            "JetHom(One, Cusp, Axis, 1){u^0*e(0,0, 3,0) + u^1*e(0,0, 7,0)}\n",
            "JetPage(One, Cusp, Axis, 2, 0, 2){g0}\n",
            "jetdim[One, Cusp, Transverse, 2, 1, 0] = 0\n",
            "jetdim[One, Cusp, Transverse, 3, 1, 0] = 0\n",
            "JetPage(One, Cusp, Curve, 7, 1, 0){g0}\n",
        )
    );
}

#[test]
fn stabilized_bivariant_page_recovers_jet_homology_dimension() {
    let declarations = include_str!("../examples/jet-pages.asl")
        .split("\nlet a1")
        .next()
        .unwrap();
    for (from, arc, expected) in [
        ("One", "Axis", 6),
        ("One", "Transverse", 4),
        ("One", "Curve", 18),
        ("Cusp", "Axis", 20),
        ("Cusp", "Curve", 52),
    ] {
        let mut program = format!("{declarations}\nrvjet[{from},Cusp,{arc}];\n");
        for parity in 0..2 {
            for filtration in 0..=7 {
                writeln!(
                    program,
                    "jetdim[{from},Cusp,{arc},7,{parity},{filtration}];"
                )
                .unwrap();
            }
        }
        let result = stdout(interpret(&program));
        let mut lines = result.lines();
        assert!(lines
            .next()
            .unwrap()
            .contains(&format!("dimension {expected},")));
        let sum: usize = lines
            .map(|line| line.rsplit_once(" = ").unwrap().1.parse::<usize>().unwrap())
            .sum();
        assert_eq!(sum, expected, "{from} to Cusp along {arc}");
    }
}

#[test]
fn field_coefficients_and_filtered_composition_satisfy_leibniz() {
    let source = r"
        qmodule T[0] {a:1,b:0} {0(a)=>[b]};
        rvpath Line[[0],7,2,[[0,2]]];
        let f = jetClass[1,0](jetterm[One,0,0,T,1,0,Line,0]);
        let g = jetClass[1,0](jetterm[T,0,0,One,0,0,Line,0]);
        let h = jetClass[1,0](jetterm[T,1,0,One,0,0,Line,0]);
        let e = jetgen[One,One,Line,1,0,0,0];
        let u = jetgen[One,One,Line,1,0,1,0];
        assert ecomp(h,f) == e;
        assert ecomp(u,u) == zero[JetPage(One,One,Line,1,0,2)];
        assert ecomp(rvScale[2](h),rvScale[3](f)) == e;
        assert jetD(jetD(f)) == zero[JetPage(One,T,Line,1,1,2)];
        assert ecomp(g,jetD(f)) == rvScale[2](u);
        assert add(ecomp(jetD(g),f),ecomp(g,jetD(f))) == jetD(ecomp(g,f));
        assert jetD(rvScale[3](f)) == rvScale[3](jetD(f));
        assert jetClass[1,0](jetRep(rvScale[2](f))) == rvScale[2](f);
        assert jetNext(jetD(f)) == zero[JetPage(One,T,Line,2,0,1)];
        assert jetPrev(jetNext(e)) == e;
        let id = \x : JetPage(One,One,Line,1,0,0) => x;
        assert id(rvScale[2](e)) == rvScale[2](e);
        card(\x : JetPage(One,One,Line,1,0,0) => S.hi);
    ";
    assert_eq!(stdout(interpret(source)), "card = 4 / 4\n");
}

#[test]
fn page_condition_and_transition_errors_are_local() {
    let declarations = r"
        qmodule T[0] {a:1,b:0} {0(a)=>[b]};
        rvpath Line[[0],7,2,[[0,2]]];
        let a = jetterm[One,0,0,T,1,0,Line,0];
        let e1 = jetClass[1,0](a);
    ";
    for (term, reason) in [
        ("jetClass[2,0](a)", "jet page cycle condition"),
        ("jetClass[1,1](a)", "below the requested jet filtration"),
        ("jetNext(e1)", "requires a cycle on the current page"),
        (
            "jetPrev(jetClass[0,0](a))",
            "requires a positive page number",
        ),
        (
            "jetgen[One,T,Line,3,1,0,0]",
            "jetgen index 0 is out of range for rank 0",
        ),
        (
            "ecomp(e1,jetgen[One,One,Line,0,0,0,0])",
            "JetPage classes on the same page",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{declarations}\nshow {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
}

#[test]
fn page_zero_is_the_actual_residue_differential() {
    let source = r"
        qmodule T[0] {a:1,b:0} {0(a)=>[b]};
        rvpath Unit[[0],7,2,[[3,2]]];
        let x = jetgen[One,T,Unit,0,1,1,0];
        let y = jetgen[One,T,Unit,0,0,1,0];
        assert jetD(x) == rvScale[3](y);
        assert jetNext(y) == zero[JetPage(One,T,Unit,1,0,1)];
        jetdim[One,T,Unit,1,1,1];
        jetdim[One,T,Unit,1000000,0,1];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "jetdim[One, T, Unit, 1, 1, 1] = 0\njetdim[One, T, Unit, 1000000, 0, 1] = 0\n"
    );
}

#[test]
fn zero_differential_pages_can_read_a_slice_at_a_large_order() {
    let value = usize::MAX;
    let source = format!(
        "rvpath Z[[],3,{value},[]]; jetdim[One,One,Z,{value},0,0]; jetdim[One,One,Z,{value},0,{}];",
        value - 1
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "jetdim[One, One, Z, {value}, 0, 0] = 1\njetdim[One, One, Z, {value}, 0, {}] = 1\n",
            value - 1
        )
    );
}
