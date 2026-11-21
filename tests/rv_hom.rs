mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cusp_specialization_preserves_extension_coefficients_and_composition() {
    assert_eq!(
        stdout(example("rv-specialization")),
        concat!(
            "RVHom(One, Cusp, [0, 1, 2], 7, [1, 3, 1]){g0}\n",
            "RVHom(One, Cusp, [0, 1, 2], 7, [1, 3, 1]){3*g0}\n",
            "RVHom(One, Cusp, [0, 1, 2], 7, [1, 3, 1]){2*g0}\n",
            "RVHom(Cusp, Cusp, [0, 1, 2], 7, [2, 3, 1]){0}\n",
        ),
    );
}

#[test]
fn raw_specialization_checks_the_evaluated_chain_map_at_each_point() {
    let module = "qmodule T[0] { a:1,b:0 } { 0(a)=>[b] };";
    for raw in ["qterm[One,0,0,T,1,0,[0]]{0}", "qterm[T,0,0,One,0,0,[0]]{0}"] {
        let failed = format!("{module}\nshow rvSpec[7,[3]]({raw});");
        let error = diagnostic(interpret(&failed));
        assert!(
            error.contains("requires a chain map at the given Margolis point"),
            "{error}"
        );
    }
    let valid = format!(
        "{module}\n{}",
        r"
        let a = qterm[One,0,0,T,1,0,[0]]{0};
        assert rvSpec[7,[0]](a) == rvgen[One,T,[0],7,[0],1];
        assert rvSpec[7,[3]](qD(a)) == zero[RVHom(One,T,[0],7,[3])];
        let weighted = qterm[One,0,0,T,1,0,[0]]{4};
        assert rvSpec[7,[0]](weighted) == zero[RVHom(One,T,[0],7,[0])];
    "
    );
    assert_eq!(stdout(interpret(&valid)), "");
}

#[test]
fn margolis_maps_are_first_class_with_full_field_digits() {
    let program = r"
        let unit = rvId[One,[0],7,[2]];
        let twice = rvScale[2](unit);
        let thrice = rvScale[3](unit);
        assert twice == add(unit, thrice);
        assert ecomp(twice, thrice) == unit;
        assert ecomp(twice, twice) == thrice;
        let translate = \x : RVHom(One,One,[0],7,[2]) => add(x, twice);
        assert translate(translate(unit)) == unit;
        assert translate(unit) == thrice;
        assert rvgen[One,One,[0],7,[2],0] == unit;
        betti RVHom(One,One,[0],7,[2]);
        let unitEmpty = rvSpec[11,[]](qid[One,[]]);
        assert unitEmpty == rvId[One,[],11,[]];
        betti RVHom(One,One,[],11,[]);
        qmodule Two[] { a:1,b:0 } {};
        betti RVHom(Two,One,[],7,[]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "betti(RVHom(One, One, [0], 7, [2])) = [4]\n",
            "betti(RVHom(One, One, [], 11, [])) = [8]\n",
            "betti(RVHom(Two, One, [], 7, [])) = [16]\n",
        )
    );
}

#[test]
fn affine_points_and_laurent_units_are_checked_before_evaluation() {
    for (program, expected) in [
        (
            "show rvSpec[7,[0]](vLocal(qid[One,[0]]));",
            "nonzero top coefficient",
        ),
        (
            "show rvSpec[7,[0]](vLocal(qRep(qid[One,[0]])));",
            "nonzero top coefficient",
        ),
        (
            "show ecomp(rvId[One,[0],7,[1]],rvId[One,[0],7,[2]]);",
            "same affine Margolis field point",
        ),
        (
            "show ecomp(rvId[One,[0],7,[1]],rvId[S,[0],7,[1]]);",
            "composable Margolis endpoints",
        ),
        (
            "show rvScale[4](rvId[One,[0],7,[1]]);",
            "field coefficient 4 is out of range",
        ),
        ("show rvgen[One,One,[0],7,[1],1];", "generator index"),
        (
            "show rvSpec[7,[1]](extgen[1,1,0]);",
            "requires a polynomial or localized Koszul",
        ),
        ("show rvId[One,[0],5,[1]];", "reducible over F2"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}

#[test]
fn topological_specialization_uses_the_steenrod_action() {
    let plane = include_str!("../examples/projective-plane.asl")
        .split("\nbetti RP2;")
        .next()
        .unwrap();
    let program = format!(
        "{plane}\n{}",
        r"
        let point = qterm[One,0,0,RP2,0,0,[0,1]]{0,0};
        let project = qterm[RP2,0,0,One,0,0,[0,1]]{0,0};
        let middle = qterm[One,0,0,RP2,1,0,[0,1]]{0,0};
        let top = qterm[One,0,0,RP2,2,0,[0,1]]{0,0};
        let atPoint = rvSpec[7,[3,0]](qClass(point));
        let atProject = rvSpec[7,[3,0]](qClass(project));
        assert ecomp(atProject, atPoint) == rvId[One,[0,1],7,[3,0]];
        assert ecomp(atPoint, atProject) == rvId[RP2,[0,1],7,[3,0]];
        assert rvSpec[7,[3,0]](qClass(middle)) == zero[RVHom(One,RP2,[0,1],7,[3,0])];
        assert rvSpec[7,[0,1]](top) == rvgen[One,RP2,[0,1],7,[0,1],2];
    "
    );
    assert_eq!(stdout(interpret(&program)), "");
}
