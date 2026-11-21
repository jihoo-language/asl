mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn restricted_may_lift_detects_the_nonzero_v1_fourth_class() {
    assert_eq!(
        stdout(example("chromatic-restriction")),
        concat!(
            "extdim[4, 12] = 0\n",
            "extandim[1, 4, 12] = 1\n",
            "extandim[0, 4, 12] = 0\n",
            "ExtAn(1, 4, 12){g0}\n",
        )
    );
}

#[test]
fn finite_profiles_form_cobar_complexes_and_restrict_transitively() {
    let program = r"
        let x = cword{[3, 1]};
        let xi2 = cword{[0, 1]};
        let a = cword{[1]};
        assert resA[0](xi2) == zero[CobarAn(0, 1, 3)];
        assert cobarD(resA[1](x)) == resA[1](cobarD(x));
        assert cobarD(cobarD(resA[1](x))) == zero[CobarAn(1, 3, 6)];
        assert cobarD(resA[1](xi2)) == resA[1](cword{[2], [1]});
        let exact = cobarD(resA[1](xi2));
        assert extClass(exact) == zero[ExtAn(1, 2, 3)];
        assert cobarD(extPrimitive(exact)) == exact;
        assert resA[1](concat(a, xi2)) == concat(resA[1](a), resA[1](xi2));
        assert resA[0](resA[1](concat(a, xi2))) == resA[0](concat(a, xi2));
        assert resA[0](resA[1](a)) == resA[0](a);
        assert resA[1](resA[2](x)) == resA[1](x);
        assert resA[1](resA[1](a)) == resA[1](a);
        assert add(resA[1](x), resA[1](x)) == zero[CobarAn(1, 1, 6)];
        betti CobarAn(1, 1, 3);
        betti CobarAn(1, 100, 600);
        betti CobarAn(1, 100, 601);
        betti CobarAn(0, 4, 4);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "betti(CobarAn(1, 1, 3)) = [4]\n",
            "betti(CobarAn(1, 100, 600)) = [2]\n",
            "betti(CobarAn(1, 100, 601)) = [1]\n",
            "betti(CobarAn(0, 4, 4)) = [2]\n",
        )
    );
}

#[test]
fn finite_adams_classes_support_towers_products_and_representatives() {
    let program = r"
        let h0 = extgen[1, 1, 0];
        let h1 = extgen[1, 2, 0];
        let h2 = extgen[1, 4, 0];
        let a = resA[1](h0);
        let b = resA[1](h1);
        assert a == extangen[1, 1, 1, 0];
        assert b == extangen[1, 1, 2, 0];
        assert resA[1](h2) == zero[ExtAn(1, 1, 4)];
        assert resA[0](resA[1](h0)) == resA[0](h0);
        assert resA[0](b) == zero[ExtAn(0, 1, 2)];
        assert yoneda(a, b) == zero[ExtAn(1, 2, 3)];
        assert yoneda(b, a) == zero[ExtAn(1, 2, 3)];
        let square = yoneda(b, b);
        assert square == extangen[1, 2, 4, 0];
        assert square == resA[1](yoneda(h1, h1));
        assert yoneda(square, b) == zero[ExtAn(1, 3, 6)];
        assert yoneda(extangen[0, 5, 5, 0], extangen[0, 2, 2, 0]) == extangen[0, 7, 7, 0];
        assert extRep(resA[0](h0)) == resA[0](cword{[1]});
        let roundtrip = \x : ExtAn(1, 2, 4) => extClass(extRep(x));
        assert roundtrip == (\x : ExtAn(1, 2, 4) => x);
        extandim[0, 100, 100];
        extandim[0, 100, 101];
        extandim[1, 3, 7];
        extandim[2, 3, 6];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "extandim[0, 100, 100] = 1\nextandim[0, 100, 101] = 0\n",
            "extandim[1, 3, 7] = 1\nextandim[2, 3, 6] = 1\n",
        )
    );
}

#[test]
fn profile_mismatches_and_impossible_grades_are_bounded() {
    for (term, message) in [
        ("resA[1](resA[0](cword{[1]}))", "smaller Steenrod algebra"),
        ("resA[1](extangen[0, 1, 1, 0])", "smaller Steenrod algebra"),
        ("resA[1](One.pt)", "cochain or Ext class"),
        (
            "concat(resA[0](cword{[1]}), resA[1](cword{[1]}))",
            "same Steenrod algebra profile",
        ),
        (
            "yoneda(extgen[1, 1, 0], extangen[1, 1, 1, 0])",
            "same Steenrod algebra profile",
        ),
        (
            "add(extangen[0, 1, 1, 0], extangen[1, 1, 1, 0])",
            "expected type",
        ),
        ("extangen[1, 1, 4, 0]", "generator index"),
        ("extClass(resA[1](cword{[0, 1]}))", "closed cobar cochain"),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(message), "{term}\n{error}");
    }

    let n = usize::MAX;
    let program = format!(
        "extandim[0, 1, {n}];\n\
         extandim[{n}, 1, 4];\n\
         assert resA[1](cword{{[{n}]}}) == zero[CobarAn(1, 1, {n})];"
    );
    assert_eq!(
        stdout(interpret(&program)),
        format!("extandim[0, 1, {n}] = 0\nextandim[{n}, 1, 4] = 1\n")
    );
}
