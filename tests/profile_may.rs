mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_full_seventh_differential_restricts_to_the_surviving_a1_periodicity_class() {
    assert_eq!(
        stdout(example("profile-may")),
        concat!(
            "mayandim[0, 8, 4, 4, 0] = 1\n",
            "MayAn(1, 1, 2, 3, 1){g0}\n",
            "May(7, 5, 12, 7){g0}\n",
            "MayAn(1, 7, 5, 12, 7){0}\n",
            "MayAn(1, 8, 4, 12, 0){g0}\n",
            "ExtAn(1, 4, 12){g0}\n",
            "mayandim[1, 20, 4, 12, 0] = 1\n",
        )
    );
}

#[test]
fn restriction_commutes_with_products_differentials_and_iteration() {
    let program = r"
        let h1 = mayClass[1, 1](cword{[2]});
        let h1square = mayMul(h1, h1);
        let a1h1 = resA[1](h1);
        let a1square = resA[1](h1square);
        assert a1square == mayMul(a1h1, a1h1);
        assert a1square == mayangen[1, 1, 2, 4, 2, 0];
        assert a1square == resA[1](resA[2](h1square));
        assert a1square == resA[1](a1square);
        assert mayD(a1square) == resA[1](mayD(h1square));
        let unit = mayClass[1, 0](resA[1](cword{}));
        assert mayMul(unit, a1square) == a1square;
        assert mayMul(a1square, unit) == a1square;
        assert mayRep(a1square) == resA[1](cword{[2], [2]});
        let translate = \x : MayAn(1, 1, 2, 4, 2) => add(x, a1square);
        assert translate(a1square) == zero[MayAn(1, 1, 2, 4, 2)];
        assert translate(zero[MayAn(1, 1, 2, 4, 2)]) == a1square;
        betti MayAn(1, 1, 2, 4, 2);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(MayAn(1, 1, 2, 4, 2)) = [2]\n"
    );
}

#[test]
fn finite_page_queries_use_profile_bounds_and_preserve_types() {
    let max = usize::MAX;
    assert_eq!(
        stdout(interpret(&format!(
            "mayandim[0, {max}, 128, 128, 0];\n\
             mayandim[0, 2, 1, {max}, 0];\n\
             mayandim[{max}, 0, 1, 2, 1];"
        ))),
        format!(
            "mayandim[0, {max}, 128, 128, 0] = 1\n\
             mayandim[0, 2, 1, {max}, 0] = 0\n\
             mayandim[{max}, 0, 1, 2, 1] = 1\n"
        )
    );

    for (command, expected) in [
        (
            "show resA[2](zero[MayAn(1, 0, 0, 0, 0)]);".into(),
            "resA cannot restrict from a smaller Steenrod algebra",
        ),
        (
            "show mayMul(zero[May(0, 1, 1, 0)], zero[MayAn(1, 0, 1, 1, 0)]);".into(),
            "same Steenrod algebra profile",
        ),
        (
            "show mayangen[0, 2, 1, 2, 1, 0];".into(),
            "generator index 0 is out of range",
        ),
        (
            "show mayClass[2, 0](resA[1](cword{[0, 1]}));".into(),
            "May page cycle condition",
        ),
        (
            format!("show mayD(zero[MayAn(1, {max}, 1, 2, 1)]);"),
            "mayD degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&command));
        assert!(error.contains(expected), "{error}");
    }
}
