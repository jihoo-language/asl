mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn coefficient_coactions_and_leibniz_are_visible_on_the_associated_graded() {
    assert_eq!(
        stdout(example("map-may")),
        concat!(
            "mapmaydim[One, RP2, 0, 0, 2, 0] = 1\n",
            "mapmaydim[One, RP2, 1, 0, 2, 0] = 0\n",
            "mapmaydim[One, Wedge, 1, 0, 2, 0] = 1\n",
            "mapmaydim[RP2, One, 1, 0, -1, 0] = 0\n",
            "mapmaydim[Wedge, One, 1, 0, -1, 0] = 1\n",
            "mapmaydim[Wedge, One, 9, 0, -1, 0] = 1\n",
            "mapextdim[Wedge, One, 0, -1] = 1\n",
            "mapmaydim[Susp(Wedge), One, 9, 1, -1, 0] = 1\n",
            "mapmaydim[Susp(Wedge), One, 9, 1, -1, 1] = 1\n",
            "mapextdim[Susp(Wedge), One, 1, -1] = 2\n",
        )
    );
}

#[test]
fn sphere_morphisms_agree_with_the_scalar_may_differentials_and_lifts() {
    let program = r"
        let xi2 = cword{[0, 1]};
        let a = mayClass[1, 0](mword[One, 0, 0, One, 0, 0](xi2));
        let h0 = mapmaygen[One, One, 1, 1, 1, 0, 0];
        let h1 = mapmaygen[One, One, 1, 1, 2, 1, 0];
        assert mayD(a) == ecomp(h0, h1);
        assert mayRep(mayD(a)) == mword[One, 0, 0, One, 0, 0](mayRep(mayD(mayClass[1, 0](xi2))));

        let sphere = maygen[7, 4, 12, 0, 0];
        let map = mapmaygen[One, One, 7, 4, 12, 0, 0];
        assert mayRep(map) == mword[One, 0, 0, One, 0, 0](mayRep(sphere));
        let image = mayD(map);
        assert image == mapmaygen[One, One, 7, 5, 12, 7, 0];
        assert mayRep(image) == mword[One, 0, 0, One, 0, 0](mayRep(mayD(sphere)));
        assert mayD(image) == zero[MapMay(One, One, 7, 6, 12, 14)];
        mapmaydim[One, One, 8, 4, 12, 0];
        mapmaydim[One, One, 8, 5, 12, 7];
        betti MapMay(One, One, 1, 1, 1, 0);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "mapmaydim[One, One, 8, 4, 12, 0] = 0\n",
            "mapmaydim[One, One, 8, 5, 12, 7] = 0\n",
            "betti(MapMay(One, One, 1, 1, 1, 0)) = [2]\n",
        )
    );
}

#[test]
fn map_may_checks_filtered_cycles_endpoints_pages_and_signed_grades() {
    let max = usize::MAX;
    let min_internal = i64::MIN;
    let max_internal = i64::MAX;
    assert_eq!(
        stdout(interpret(&format!(
            "mapmaydim[One, One, 0, 0, {max_internal}, 0];\n\
             mapmaydim[One, One, 1, 1, {min_internal}, 0];\n\
             mapmaydim[One, One, 0, 1, {max_internal}, {max}];"
        ))),
        format!(
            "mapmaydim[One, One, 0, 0, {max_internal}, 0] = 0\n\
             mapmaydim[One, One, 1, 1, {min_internal}, 0] = 0\n\
             mapmaydim[One, One, 0, 1, {max_internal}, {max}] = 0\n"
        )
    );
    for (command, expected) in [
        (
            "show mayClass[0, 1](mword[One, 0, 0, One, 0, 0](cword{[1]}));".into(),
            "below the requested May filtration",
        ),
        (
            "show mayClass[2, 0](mword[One, 0, 0, One, 0, 0](cword{[0, 1]}));".into(),
            "May page cycle condition",
        ),
        (
            "show ecomp(zero[MapMay(One, One, 0, 0, 0, 0)], zero[MapMay(One, One, 1, 0, 0, 0)]);".into(),
            "same page",
        ),
        (
            "show ecomp(zero[MapMay(One, One, 0, 0, 0, 0)], zero[MapMay(One, S, 0, 0, 0, 0)]);".into(),
            "endpoint",
        ),
        (
            format!("show ecomp(zero[MapMay(One, One, 0, 0, {max_internal}, 0)], zero[MapMay(One, One, 0, 0, 1, 0)]);"),
            "internal degree is out of range",
        ),
        (
            format!("show mayD(zero[MapMay(One, One, {max}, 1, 2, 1)]);"),
            "mayD degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&command));
        assert!(error.contains(expected), "{error}");
    }
}
