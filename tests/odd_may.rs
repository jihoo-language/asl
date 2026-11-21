mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn odd_primary_associated_graded_and_first_differential_keep_coefficients() {
    let source = "
        let xi = pword[3]{([1], [])};
        let square = pword[3]{([2], [])};
        let leading = mayClass[0, 6](square);
        assert mayD(leading) == pscale[2](mayClass[0, 6](concat(xi, xi)));
        assert mayClass[0, 6](mayRep(pscale[2](leading))) == pscale[2](leading);
        assert mayNext(mayD(leading)) == zero[PMay(3, 1, 2, 8, 6)];
        let t = pword[3]{([], [1])};
        let e0 = mayClass[0, 2](t);
        let e1 = pmaygen[3, 1, 1, 5, 2, 0];
        assert mayD(e0) == zero[PMay(3, 0, 2, 5, 2)];
        assert mayNext(e0) == e1;
        assert mayPrev(e1) == e0;
        assert mayNext(mayPrev(pscale[2](e1))) == pscale[2](e1);
        let image = mayD(e1);
        assert image == mayClass[1, 3](pword[3]{([1], []), ([], [0])});
        assert mayD(image) == zero[PMay(3, 1, 3, 5, 4)];
        assert mayNext(image) == zero[PMay(3, 2, 2, 5, 3)];
    ";
    assert!(stdout(interpret(source)).is_empty());
}

#[test]
fn the_second_and_fifth_differentials_remove_classes_on_both_ends() {
    let output = stdout(example("odd-may"));
    assert!(output.starts_with("PMay(3, 1, 2, 5, 3){g0}\n"));
    assert!(output.ends_with("pmaydim[3, 6, 3, 15, 6] = 0\npmaydim[3, 6, 4, 15, 11] = 0\n"));
    let source = "
        let a = pmaygen[3, 2, 3, 14, 7, 0];
        let z = pword[3]{([], [0]), ([], [0])};
        let b = add(pword[3]{([1], []), ([2], [])}, pword[3]{([2], []), ([1], [])});
        let target = mayClass[2, 9](concat(z, b));
        assert mayD(a) == pscale[2](target);
        assert mayNext(mayD(a)) == zero[PMay(3, 3, 4, 14, 9)];
        pmaydim[3, 3, 3, 14, 7];
        pmaydim[3, 3, 4, 14, 9];
        let higher = pmaygen[3, 5, 3, 15, 6, 0];
        assert mayNext(mayPrev(higher)) == higher;
        assert mayClass[5, 6](mayRep(pscale[2](higher))) == pscale[2](higher);
        let permanent = mayClass[100, 9](b);
        assert permanent == pmaygen[3, 100, 2, 12, 9, 0];
        assert extClass(mayRep(permanent)) == pextgen[3, 2, 12, 0];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "pmaydim[3, 3, 3, 14, 7] = 0\npmaydim[3, 3, 4, 14, 9] = 0\n"
    );
}

#[test]
fn five_primary_page_products_see_the_odd_stem_sign_and_respect_transitions() {
    let source = r"
        let x = mayClass[1, 7](pword[5]{([1], [])});
        let y = mayClass[1, 39](pword[5]{([5], [])});
        let xy = mayMul(x, y);
        assert xy == pscale[4](mayMul(y, x));
        assert xy == pmaygen[5, 1, 2, 48, 46, 0];
        assert mayMul(pscale[3](x), pscale[2](y)) == xy;
        assert mayNext(xy) == mayMul(mayNext(x), mayNext(y));
        assert mayNext(xy) == zero[PMay(5, 2, 2, 48, 46)];
        let passage = \a : PMay(5, 1, 1, 8, 7) => mayPrev(mayNext(a));
        assert passage == (\a : PMay(5, 1, 1, 8, 7) => a);
        betti PMay(5, 1, 2, 48, 46);
    ";
    assert_eq!(
        stdout(interpret(source)),
        "betti(PMay(5, 1, 2, 48, 46)) = [5]\n"
    );
}

#[test]
fn individual_classes_work_over_large_prime_fields_without_enumeration() {
    let source = "
        let unit = pmaygen[4294967291, 1, 0, 0, 0, 0];
        let minus = pscale[4294967290](unit);
        assert add(unit, minus) == zero[PMay(4294967291, 1, 0, 0, 0)];
        assert mayMul(minus, minus) == unit;
        assert mayPrev(mayNext(minus)) == minus;
    ";
    assert!(stdout(interpret(source)).is_empty());
    assert!(diagnostic(interpret("betti PMay(4294967291, 1, 0, 0, 0);")).contains("limit"));
}

#[test]
fn prime_page_filtration_and_cycle_checks_report_errors() {
    for (expression, expected) in [
        ("zero[PMay(2, 1, 0, 0, 0)]".into(), "odd prime"),
        ("zero[PMay(9, 1, 0, 0, 0)]".into(), "not prime"),
        ("pmaygen[3, 1, 1, 5, 2, 1]".into(), "generator index"),
        ("mayClass[1, 1](pword[3]{([], [0])})".into(), "below"),
        (
            "mayClass[1, 6](pword[3]{([2], [])})".into(),
            "cycle condition",
        ),
        (
            "mayClass[2, 2](pword[3]{([], [1])})".into(),
            "cycle condition",
        ),
        ("mayNext(pmaygen[3, 1, 1, 5, 2, 0])".into(), "cycle for d_1"),
        (
            "mayPrev(pmaygen[3, 0, 1, 5, 2, 0])".into(),
            "greater than zero",
        ),
        (
            "mayMul(zero[PMay(3, 1, 0, 0, 0)], zero[PMay(5, 1, 0, 0, 0)])".into(),
            "prime",
        ),
        (
            "mayMul(zero[PMay(3, 1, 0, 0, 0)], zero[PMay(3, 2, 0, 0, 0)])".into(),
            "page",
        ),
        (
            format!("mayNext(zero[PMay(3, {}, 0, 0, 0)])", usize::MAX),
            "out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(expected), "{expression}\n{error}");
    }
}
