mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cyclic_cover_changes_components_and_global_sections() {
    assert_eq!(
        stdout(example("set-families")),
        concat!(
            "betti(Cover) = [1, 1]\n",
            "betti(Three) = [3, 3]\n",
            "setsections[Twisted] = 0\n",
            "setsections[Trivial] = 3\n",
        )
    );
}

#[test]
fn derived_posets_support_tables_and_natural_maps_are_finite_typed_values() {
    let output = stdout(interpret(
        r"
        space Chain { x < y < z };
        setfamily D[Chain] {x:2,y:1,z:2}
            {x<y:[0,0],y<z:[1],x<z:[1,1]};
        setfamily T[Chain] {x:1,y:1,z:1} {x<y:[0],y<z:[0]};
        elements E[D];
        elements A[T];
        setnat collapse[D,T] {y:[0],z:[0,0],x:[0,0]};
        setnat chosen[T,D] {x:[1],y:[0],z:[1]};
        assert setComp(collapse,chosen) == setId[T];
        setnat constant[D,D] {x:[1,1],y:[0],z:[1,1]};
        assert setComp(chosen,collapse) == constant;
        assert setProjection[E] == tab E -> Chain
            {x_0 => Chain.x, x_1 => Chain.x, y_0 => Chain.y, z_0 => Chain.z, z_1 => Chain.z};
        let spaceMap = setTotal[E,A](collapse);
        assert spaceMap(E.x_1) == A.x_0;
        assert setTotal[A,E](chosen)(A.z_0) == E.z_1;
        let typed = \f : SetNat(D,T) => f;
        assert typed(collapse) == collapse;
        let all = \f : SetNat(T,D) => S.hi;
        card all;
        setsections[D];
    ",
    ));
    assert_eq!(output, "card = 2 / 2\nsetsections[D] = 2\n");
}

#[test]
fn empty_fibers_and_sites_have_correct_section_and_map_types() {
    let output = stdout(interpret(
        r"
        space Empty {};
        setfamily Vacuous[Empty] {} {};
        elements V[Vacuous];
        setsections[Vacuous];
        assert setProjection[V] == tab V -> Empty {};
        assert setsect[V,0] == tab Empty -> V {};
        setfamily F[S] {lo:0,hi:1} {lo<hi:[]};
        setfamily G[S] {lo:0,hi:0} {lo<hi:[]};
        setsections[F];
        elements E[F];
        assert setProjection[E] == tab E -> S {hi_0 => S.hi};
        let empty = \f : SetNat(F,G) => S.hi;
        card empty;
    ",
    ));
    assert_eq!(
        output,
        "setsections[Vacuous] = 1\nsetsections[F] = 0\ncard = 0 / 0\n",
    );
}

#[test]
fn declarations_reject_inconsistent_diagrams_and_maps() {
    let diamond = "space X {x<y<w,x<z<w};";
    let cases = [
        (
            format!(
                "{diamond}setfamily F[X] {{x:1,y:1,z:1,w:2}} {{x<y:[0],x<z:[0],y<w:[0],z<w:[1]}};"
            ),
            "transitions disagree",
        ),
        (
            "setfamily F[S] {lo:1,hi:2} {};".into(),
            "missing set-family transition",
        ),
        (
            "setfamily F[S] {lo:1,hi:0} {lo<hi:[0]};".into(),
            "outside its target fiber",
        ),
        (
            "setfamily F[S] {lo:2,hi:2} {lo<hi:[0,1]}; setnat u[F,F] {lo:[0,1],hi:[1,0]};".into(),
            "does not commute",
        ),
        (
            "setfamily F[S] {lo:2,hi:2} {lo<hi:[0,1]}; elements E[F]; show setsect[E,2];".into(),
            "outside its range",
        ),
        (
            "setfamily F[One] {pt:513} {}; elements E[F];".into(),
            "512 points",
        ),
        ("setfamily F[One] {pt:4097} {};".into(), "4096 occurrences"),
        (
            "setfamily F[One] {pt:1} {}; setfamily G[S] {lo:0,hi:0} {lo<hi:[]}; setnat u[F,G] {};"
                .into(),
            "same named site",
        ),
    ];
    for (script, expected) in cases {
        let error = diagnostic(interpret(&script));
        assert!(error.contains(expected), "{error}");
    }
}
