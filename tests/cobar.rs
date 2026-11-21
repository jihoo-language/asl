mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn milnor_coproduct_and_derivation_use_the_normalized_cobar_convention() {
    assert_eq!(
        stdout(example("cobar")),
        concat!(
            "Cobar(2, 3){cword{[2], [1]}}\n",
            "Cobar(2, 7){cword{[4], [0, 1]} + cword{[0, 2], [1]}}\n",
            "betti(Cobar(1, 3)) = [4]\n",
        )
    );

    let program = r"
        let a = cword{[4]};
        let x = cword{[1, 1]};
        assert cobarD(a) == zero[Cobar(2, 4)];
        assert cobarD(cword{[0, 2]}) == cword{[4], [2]};
        assert cobarD(cword{[3]}) == add(cword{[1], [2]}, cword{[2], [1]});
        assert cobarD(x) ==
            add(add(cword{[3], [1]}, cword{[1], [0, 1]}),
                add(cword{[0, 1], [1]}, cword{[2], [2]}));
        assert cobarD(add(x, a)) == add(cobarD(x), cobarD(a));
        assert cobarD(cobarD(cword{[3, 5, 1], [1, 1]})) == zero[Cobar(4, 29)];
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn cobar_cochains_are_first_class_and_discrete() {
    let program = r"
        let x = cword{[3, 0, 0]};
        let y = cword{[0, 1]};
        let sum = add(x, y);
        let id = \u : Cobar(1, 3) => u;
        let translate = \u : Cobar(1, 3) => add(u, sum);
        assert id(sum) == sum;
        assert translate(translate(y)) == y;
        assert x == cword{[3]};
        assert add(sum, sum) == zero[Cobar(1, 3)];
        assert concat(x, cword{}) == x;
        assert concat(x, zero[Cobar(1, 3)]) == zero[Cobar(2, 6)];
        assert concat(zero[Cobar(5, 1)], x) == zero[Cobar(6, 4)];
        assert concat(sum, sum) ==
            add(add(cword{[3], [3]}, cword{[3], [0, 1]}),
                add(cword{[0, 1], [3]}, cword{[0, 1], [0, 1]}));
        assert concat(concat(x, y), x) == concat(x, concat(y, x));
        assert cobarD(cword{}) == zero[Cobar(1, 0)];
        assert cobarD(zero[Cobar(0, 9)]) == zero[Cobar(1, 9)];
        show cword{};
        show cword{[1, 0, 0], [0, 1, 0]};
        show zero[Cobar(1, 3)];
        betti Cobar(0, 0);
        betti Cobar(0, 30);
        betti Cobar(5, 4);
        betti Cobar(2, 4);
        card (\u : Cobar(3, 4) => S.hi);
        homotopy (\p : One => x) ~ (\p : One => y);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "Cobar(0, 0){cword{}}\nCobar(2, 4){cword{[1], [0, 1]}}\nCobar(1, 3){0}\n",
            "betti(Cobar(0, 0)) = [2]\nbetti(Cobar(0, 30)) = [1]\n",
            "betti(Cobar(5, 4)) = [1]\nbetti(Cobar(2, 4)) = [32]\n",
            "card = 8 / 8\nhomotopy = no\n",
        )
    );
}

#[test]
fn cobar_names_are_contextual() {
    let program = r"
        space Cobar { a };
        let cword = \p : Cobar => p;
        let cobarD = cword;
        let concat = cword;
        assert concat(cobarD(Cobar.a)) == Cobar.a;
        assert (\p : Cobar(0, 0) => p)(cword{}) == cword{};
        assert concat(cword{}, cword{}) == cword{};
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn cobar_rejects_degenerate_factors_and_checks_both_grades() {
    for literal in ["cword{[]}", "cword{[0, 0]}", "cword{[1], [0]}"] {
        assert!(diagnostic(interpret(&format!("show {literal};"))).contains("positive degree"));
    }
    for input in [
        "add(cword{[2]}, cword{[1], [1]})",
        "add(cword{[2]}, cword{[1]})",
        "concat(sq[0], cword{})",
        "cobarD(S.hi)",
    ] {
        let error = diagnostic(interpret(&format!("\nshow {input};")));
        assert!(error.contains(":2:"), "{input}\n{error}");
    }
    for input in [
        format!("cword{{[{}, 1]}}", usize::MAX),
        format!("cword{{[{}], [1]}}", usize::MAX),
        format!("concat(cword{{[{}]}}, cword{{[1]}})", usize::MAX),
        format!("concat(zero[Cobar({}, 0)], cword{{[1]}})", usize::MAX),
        format!("cobarD(zero[Cobar({}, 0)])", usize::MAX),
    ] {
        let error = diagnostic(interpret(&format!("show {input};")));
        assert!(error.contains("out of range"), "{input}\n{error}");
    }
}

#[test]
fn large_sparse_products_work_without_enumeration_and_dense_work_is_bounded() {
    let n = 1usize << (usize::BITS - 1);
    assert_eq!(
        stdout(interpret(&format!(
            "let x = cword{{[{n}]}};
             assert cobarD(x) == zero[Cobar(2, {n})];
             assert concat(cword{{}}, x) == x;
             show x;"
        ))),
        format!("Cobar(1, {n}){{cword{{[{n}]}}}}\n"),
    );
    assert!(diagnostic(interpret(&format!("betti Cobar(1, {n});"))).contains("enumeration limit"));
    assert!(
        diagnostic(interpret(&format!("betti Cobar({n}, {n});"))).contains("word storage limit")
    );
    assert!(diagnostic(interpret(&format!(
        "show cobarD(cword{{[{}]}});",
        usize::MAX
    )))
    .contains("limit"));
}
