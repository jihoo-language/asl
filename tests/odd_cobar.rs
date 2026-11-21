mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn odd_primary_cobar_detects_both_koszul_signs() {
    let output = stdout(example("odd-cobar"));
    assert!(output.contains("PCobar(3, 2, 8){2 pword[3]"), "{output}");
    assert!(
        output.ends_with("betti(PCobar(3, 1, 16)) = [9]\n"),
        "{output}"
    );
}

#[test]
fn five_primary_differential_respects_leibniz_and_squares_to_zero() {
    let program = r"
        let q = pword[5]{([1], [])};
        let a = pword[5]{([], [0])};
        let b = pword[5]{([], [1])};
        let s = pword[5]{([2], [])};
        let c = pword[5]{([1], [0, 1, 2])};
        let u = pword[5]{};
        assert cobarD(b) == concat(q, a);
        assert cobarD(s) == pscale[2](concat(q, q));
        assert cobarD(pword[5]{([5], [])}) == zero[PCobar(5, 2, 40)];
        assert cobarD(concat(q, b)) == pscale[4](concat(q, cobarD(b)));
        assert cobarD(concat(a, b)) == concat(a, cobarD(b));
        assert cobarD(cobarD(c)) == zero[PCobar(5, 3, 67)];
        assert cobarD(cobarD(concat(c, b))) == zero[PCobar(5, 4, 76)];
        assert concat(u, c) == c;
        assert concat(c, u) == c;
        assert add(c, pscale[4](c)) == zero[PCobar(5, 1, 67)];
        assert pscale[11](c) == c;
        assert pscale[10](c) == zero[PCobar(5, 1, 67)];
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn the_finite_space_uses_all_prime_field_digits_in_coordinate_order() {
    let program = r"
        let xi2 = pword[3]{([0, 1], [])};
        let power = pword[3]{([4], [])};
        let triple = \x : PCobar(3, 1, 16) => pscale[3](x);
        let addTwice = \x : PCobar(3, 1, 16) => add(x, x);
        let change = \x : PCobar(3, 1, 16) => add(x, add(xi2, pscale[2](power)));
        assert addTwice(add(xi2, power)) == add(pscale[2](xi2), pscale[2](power));
        assert change(change(change(power))) == power;
        assert triple(change(xi2)) == zero[PCobar(3, 1, 16)];
        betti PCobar(3, 1, 16);
        betti PCobar(5, 0, 0);
        betti PCobar(3, 2, 3);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(PCobar(3, 1, 16)) = [9]\nbetti(PCobar(5, 0, 0)) = [5]\nbetti(PCobar(3, 2, 3)) = [1]\n"
    );
}

#[test]
fn sparse_large_words_and_large_prime_coefficients_avoid_dense_enumeration() {
    let calculations = "
        let primitive = pword[3]{([59049], [])};
        let long = concat(primitive, primitive);
        assert cobarD(primitive) == zero[PCobar(3, 2, 236196)];
        assert cobarD(long) == zero[PCobar(3, 3, 472392)];
        assert long == pword[3]{([59049], []), ([59049], [])};
        let unit = pword[4294967291]{};
        let minus = pscale[4294967290](unit);
        assert concat(minus, minus) == unit;
        assert add(minus, unit) == zero[PCobar(4294967291, 0, 0)];
    ";
    assert!(stdout(interpret(calculations)).is_empty());
    let error = diagnostic(interpret("betti PCobar(4294967291, 0, 0);"));
    assert!(error.contains("limit"), "{error}");
}

#[test]
fn malformed_primes_and_word_factors_and_grades_are_diagnosed() {
    let max = usize::MAX;
    let expressions = [
        ("pword[2]{}".into(), "odd prime"),
        ("zero[PCobar(9, 0, 0)]".into(), "not prime"),
        ("pword[4294967295]{}".into(), "not prime"),
        ("pword[3]{([], [])}".into(), "positive degree"),
        ("pword[3]{([0, 0], [])}".into(), "positive degree"),
        ("pword[3]{([], [1, 0])}".into(), "increasing"),
        ("pword[3]{([], [0, 0])}".into(), "increasing"),
        (format!("pword[3]{{([{max}], [])}}"), "degree"),
        (format!("pword[3]{{([], [{max}])}}"), "degree"),
        (format!("cobarD(zero[PCobar(3, {max}, 0)])"), "degree"),
        ("concat(pword[3]{}, pword[5]{})".into(), "prime"),
        ("add(pword[3]{}, pword[5]{})".into(), "type"),
        ("pscale[2](cword{})".into(), "PCobar"),
    ];
    for (expression, expected) in expressions {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(expected), "{expression}\n{error}");
    }
    if usize::BITS > 32 {
        let error = diagnostic(interpret("show pword[4294967296]{};"));
        assert!(error.contains("u32"), "{error}");
    }
}
