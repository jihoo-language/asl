mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn a_genuine_intermediate_square_is_nonzero() {
    assert_eq!(stdout(example("cobar-squares")), "Ext(3, 10){g0}\n");
}

#[test]
fn cup_i_uses_normalization_cofaces_and_the_chain_homotopy_identities() {
    let program = "
        let u = cword{};
        let a = cword{[0, 1]};
        let b = cword{[3], [1]};
        assert cCup[0](a, b) == concat(a, b);
        assert cCup[0](u, b) == b;
        assert cCup[1](a, a) == cword{[0, 2]};
        assert cCup[1](cword{[1]}, cword{[2]}) == cword{[3]};
        assert cCup[2](a, b) == zero[Cobar(1, 7)];
        assert cobarD(cCup[1](a, b)) ==
            add(add(cCup[1](cobarD(a), b), cCup[1](a, cobarD(b))),
                add(concat(a, b), concat(b, a)));
        assert cobarD(cCup[2](a, b)) ==
            add(add(cCup[2](cobarD(a), b), cCup[2](a, cobarD(b))),
                add(cCup[1](a, b), cCup[1](b, a)));
        assert cCup[1](add(a, cword{[3]}), b) ==
            add(cCup[1](a, b), cCup[1](cword{[3]}, b));
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn algebraic_squares_are_independent_of_cobar_representatives() {
    let program = "
        let r = cword{[1], [4]};
        let delta = cobarD(cword{[5]});
        assert delta == add(cword{[1], [4]}, cword{[4], [1]});
        let swapped = add(r, delta);
        let x = extClass(r);
        assert extClass(swapped) == x;
        assert extClass(cCup[1](r, r)) == extSq[1](x);
        assert extClass(cCup[1](swapped, swapped)) == extSq[1](x);
        assert extClass(cCup[1](delta, delta)) == zero[Ext(3, 10)];
        assert extSq[4](x) == zero[Ext(6, 10)];
        assert extSq[0](extgen[0, 0, 0]) == extgen[0, 0, 0];
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn squares_commute_with_restriction_and_satisfy_cartan_and_adem() {
    let program = "
        let h0 = extgen[1, 1, 0];
        let h1 = extgen[1, 2, 0];
        let h2 = extgen[1, 4, 0];
        let x = yoneda(h0, h2);
        let y = yoneda(h0, yoneda(h0, h0));
        let middle = extSq[1](x);
        assert middle ==
            add(yoneda(extSq[1](h0), extSq[0](h2)),
                yoneda(extSq[0](h0), extSq[1](h2)));
        assert extSq[1](extSq[2](y)) == extSq[3](extSq[0](y));
        assert extSq[0](resA[1](h0)) == resA[1](h1);
        assert extSq[0](resA[1](h1)) == zero[ExtAn(1, 1, 4)];
        assert extSq[1](resA[2](x)) == resA[2](middle);
        assert extSq[1](resA[3](x)) == resA[3](middle);

        let a = cword{[0, 1]};
        let b = cword{[3], [1]};
        assert resA[1](cCup[1](a, b)) == cCup[1](resA[1](a), resA[1](b));
        let power = resA[1](cword{[2]});
        assert cCup[1](power, power) == zero[CobarAn(1, 1, 4)];
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn cup_i_and_squares_check_profiles_bidegrees_and_combinatorial_work() {
    for (expression, expected) in [
        ("cCup[3](cword{[1]}, cword{[1]})", "negative"),
        ("cCup[1](cword{[1]}, resA[0](cword{[1]}))", "profile"),
        ("cCup[1](cword{[1]}, S.hi)", "cobar"),
        ("extSq[1](cword{[1]})", "Ext"),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(expected), "{error}");
    }
    let max = usize::MAX;
    let empty = format!("show cCup[{max}](zero[Cobar({max}, 0)], zero[Cobar(1, 0)]);");
    assert_eq!(stdout(interpret(&empty)), "Cobar(1, 0){0}\n");
    let error = diagnostic(interpret(&format!("show extSq[1](zero[Ext({max}, 0)]);")));
    assert!(error.contains("out of range"), "{error}");
    let factors = vec!["[1]"; 24].join(", ");
    let error = diagnostic(interpret(&format!(
        "let a = cword{{{factors}}}; show cCup[12](a, a);"
    )));
    assert!(error.contains("limit"), "{error}");
}
