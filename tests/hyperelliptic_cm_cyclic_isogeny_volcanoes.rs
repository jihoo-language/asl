mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn volcano_walk_and_exceptional_unit_multiplicities() {
    let output = stdout(example("hyperelliptic-cm-cyclic-isogeny-volcanoes"));
    assert_eq!(
        output,
        concat!(
            "rvCMNeighbors[2] = ascending=0, horizontal=2, descending=1, targets=[",
            "{discriminant=-23, form=[2, -1, 3], direction=horizontal, multiplicity=1}, ",
            "{discriminant=-23, form=[2, 1, 3], direction=horizontal, multiplicity=1}, ",
            "{discriminant=-92, form=[1, 0, 23], direction=descending, multiplicity=1}]\n",
            "rvCMNeighborData[2,1] = address=1, discriminant=-23, form=[1, 1, 6], direction=ascending, content=4, multiplicity=1\n",
            "rvCMNeighbors[2] = ascending=1, horizontal=0, descending=2, targets=[",
            "{discriminant=-23, form=[1, 1, 6], direction=ascending, multiplicity=1}, ",
            "{discriminant=-368, form=[1, 0, 92], direction=descending, multiplicity=1}, ",
            "{discriminant=-368, form=[4, 0, 23], direction=descending, multiplicity=1}]\n",
            "rvCMNeighborCounts[2] = ascending=1, horizontal=0, descending=2\n",
            "rvCMNeighbors[2] = ascending=0, horizontal=1, descending=2, targets=[",
            "{discriminant=-4, form=[1, 0, 1], direction=horizontal, multiplicity=1}, ",
            "{discriminant=-16, form=[1, 0, 4], direction=descending, multiplicity=2}]\n",
            "rvCMNeighbors[2] = ascending=0, horizontal=0, descending=3, targets=[",
            "{discriminant=-12, form=[1, 0, 3], direction=descending, multiplicity=3}]\n",
        ),
    );
}

#[test]
fn horizontal_edges_match_independent_ideal_products_and_keep_target_names() {
    let source = "rvcmorder Top[-23]; rvcmorder Twin[-23]; rvcmorder Down[-92];
        let g=rvCMClass[Top,2,1,3];
        let h=rvCMInverse(g);
        assert rvCMNeighbor[Top,2,1](g)==rvCMCompose(g,g);
        assert rvCMNeighbor[Top,2,2](g)==rvCMCompose(g,h);
        assert rvCMNeighbor[Twin,2,2](g)==rvCMPrincipal[Twin];
        let child=rvCMNeighbor[Down,2,0](g);
        assert child==rvCMClass[Down,3,-2,8];
        rvCMNeighborData[2,2](g);";
    assert_eq!(
        stdout(interpret(source)),
        "rvCMNeighborData[2,2] = address=infinity, discriminant=-23, form=[1, 1, 6], direction=horizontal, content=2, multiplicity=1\n"
    );
}

#[test]
fn ramified_conductor_and_inert_odd_prime_neighbors() {
    let output = stdout(interpret(
        "rvcmorder E[-3]; rvcmorder ER[-27]; rvcmorder G[-4]; rvcmorder GI[-36];
         let e=rvCMPrincipal[E];
         let g=rvCMPrincipal[G];
         rvCMNeighbors[3](e); rvCMNeighbors[3](g);
         let er=rvCMNeighbor[ER,3,0](e);
         let gi=rvCMNeighbor[GI,3,0](g);
         rvCMNeighborCounts[3](er); rvCMNeighborCounts[3](gi);",
    ));
    assert_eq!(
        output,
        concat!(
            "rvCMNeighbors[3] = ascending=0, horizontal=1, descending=3, targets=[",
            "{discriminant=-3, form=[1, 1, 1], direction=horizontal, multiplicity=1}, ",
            "{discriminant=-27, form=[1, 1, 7], direction=descending, multiplicity=3}]\n",
            "rvCMNeighbors[3] = ascending=0, horizontal=0, descending=4, targets=[",
            "{discriminant=-36, form=[1, 0, 9], direction=descending, multiplicity=2}, ",
            "{discriminant=-36, form=[2, 2, 5], direction=descending, multiplicity=2}]\n",
            "rvCMNeighborCounts[3] = ascending=1, horizontal=0, descending=3\n",
            "rvCMNeighborCounts[3] = ascending=1, horizontal=0, descending=3\n",
        ),
    );
}

#[test]
fn wide_reports_and_precise_input_diagnostics() {
    let setup = format!(
        "rvcmorder O[-23]; rvcmorder Wrong[-20]; rvcmorder Wide[{}];",
        i128::MIN,
    );
    assert_eq!(
        stdout(interpret(&format!("{setup} rvCMNeighborData[3,3](rvCMPrincipal[Wide]);"))),
        concat!(
            "rvCMNeighborData[3,3] = address=infinity, discriminant=-1531270651144223085585185733442956951552, ",
            "form=[1, 0, 382817662786055771396296433360739237888], direction=descending, content=1, multiplicity=1\n"
        ),
    );
    for (command, expected) in [
        (
            "rvCMNeighborCounts[1](rvCMPrincipal[O]);",
            "isogeny degree must be prime",
        ),
        ("rvCMNeighborCounts[9](rvCMPrincipal[O]);", "9 is not prime"),
        (
            "rvCMNeighborData[2,3](rvCMPrincipal[O]);",
            "address must be between 0 and ell inclusive",
        ),
        (
            "show rvCMNeighbor[Wrong,2,0](rvCMPrincipal[O]);",
            "discriminant -23 differs from the named target order",
        ),
        (
            "show rvCMNeighbor[Missing,2,0](rvCMPrincipal[O]);",
            "unknown CM ring order",
        ),
        (
            "rvCMNeighborCounts[4294967291](rvCMPrincipal[O]);",
            "neighbor enumeration exceeds the computation limit",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}
