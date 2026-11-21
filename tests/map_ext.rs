mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn composition_tracks_signed_degrees_and_multiplies_homology_matrices() {
    assert_eq!(
        stdout(example("map-ext")),
        "mapextdim[Two, Two, 0, 0] = 4\nmapextdim[Circle, One, 0, -1] = 1\nMapExt(Circle, Circle, 0, 0){g1}\n"
    );
}

#[test]
fn bivariant_cobar_detects_the_steenrod_linkage_of_the_projective_plane() {
    let program = format!(
        "{}\n{}",
        include_str!("../examples/projective-plane.asl"),
        r"
        mapextdim[RP2, RP2, 0, 0];
        mapextdim[Wedge, Wedge, 0, 0];
        mapextdim[RP2, RP2, 1, 1];
        assert eid[RP2] == emap(\x : RP2 => x);
        assert eid[Wedge] == emap(\x : Wedge => x);
        let id = \x : MapExt(RP2, RP2, 1, 1) => x;
        let right = \x : MapExt(RP2, RP2, 1, 1) => ecomp(x, eid[RP2]);
        let left = \x : MapExt(RP2, RP2, 1, 1) => ecomp(eid[RP2], x);
        assert left == id;
        assert right == id;
        betti MapExt(RP2, RP2, 1, 1);
        "
    );
    assert!(stdout(interpret(&program)).ends_with(concat!(
        "mapextdim[RP2, RP2, 0, 0] = 2\n",
        "mapextdim[Wedge, Wedge, 0, 0] = 3\n",
        "mapextdim[RP2, RP2, 1, 1] = 3\n",
        "betti(MapExt(RP2, RP2, 1, 1)) = [8]\n",
    )));
}

#[test]
fn sphere_maps_recover_ext_and_point_sources_recover_coefficient_ranks() {
    let program = r"
        let h0 = mapextgen[One, One, 1, 1, 0];
        let h1 = mapextgen[One, One, 1, 2, 0];
        let h2 = mapextgen[One, One, 1, 4, 0];
        assert ecomp(h0, h1) == zero[MapExt(One, One, 2, 3)];
        assert ecomp(h1, h0) == zero[MapExt(One, One, 2, 3)];
        let cube = ecomp(h1, ecomp(h1, h1));
        assert cube == ecomp(h2, ecomp(h0, h0));
        assert cube == mapextgen[One, One, 3, 6, 0];
        mapextdim[One, One, 3, 6];
        extdim[3, 6];
        space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };
        mapextdim[One, Circle, 1, 2];
        xextdim[Circle, 1, 2];
        let inc = emap(\x : One => Circle.a);
        let collapse = emap(\x : Circle => One.pt);
        assert ecomp(collapse, ecomp(inc, h0)) == h0;
        assert ecomp(ecomp(inc, h0), eid[One]) == ecomp(inc, h0);
        let select = mapextgen[Circle, One, 0, -1, 0];
        assert ecomp(h0, select) == mapextgen[Circle, One, 1, 0, 0];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "mapextdim[One, One, 3, 6] = 1\nextdim[3, 6] = 1\n",
            "mapextdim[One, Circle, 1, 2] = 2\nxextdim[Circle, 1, 2] = 2\n",
        )
    );
}

#[test]
fn map_ext_checks_endpoints_and_signed_degree_arithmetic() {
    let program = format!(
        "space Void {{}};\n\
         assert eid[Void] == zero[MapExt(Void, Void, 0, 0)];\n\
         assert emap(\\x : Void => One.pt) == zero[MapExt(Void, One, 0, 0)];\n\
         mapextdim[One, One, 0, {}];\n\
         mapextdim[One, One, 0, {}];\n\
         mapextdim[One, One, 2, {}];\n\
         mapextdim[A(100) -> S, One, 0, {}];",
        i64::MAX,
        i64::MIN,
        i64::MIN,
        i64::MAX,
    );
    assert_eq!(stdout(interpret(&program)), format!(
        "mapextdim[One, One, 0, {}] = 0\nmapextdim[One, One, 0, {}] = 0\nmapextdim[One, One, 2, {}] = 0\nmapextdim[A(100) -> S, One, 0, {}] = 0\n",
        i64::MAX, i64::MIN, i64::MIN, i64::MAX,
    ));

    for (command, message) in [
        ("show ecomp(eid[S], eid[One]);".into(), "inner target One"),
        ("show emap(S.hi);".into(), "emap requires a function"),
        (
            "show mapextgen[One, One, 0, -1, 0];".into(),
            "generator index",
        ),
        (
            format!(
                "show ecomp(zero[MapExt(One, One, 0, {})], zero[MapExt(One, One, 0, 1)]);",
                i64::MAX,
            ),
            "signed internal degree is out of range",
        ),
        (
            format!(
                "show ecomp(zero[MapExt(One, One, 0, {})], zero[MapExt(One, One, 0, -1)]);",
                i64::MIN,
            ),
            "signed internal degree is out of range",
        ),
        (
            format!(
                "show ecomp(zero[MapExt(One, One, {}, 0)], zero[MapExt(One, One, 1, 0)]);",
                usize::MAX,
            ),
            "degree is out of range",
        ),
        (
            format!("mapextdim[One, One, 0, {}];", i128::from(i64::MAX) + 1),
            "signed internal degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&command));
        assert!(error.contains(message), "{command}\n{error}");
    }
}
