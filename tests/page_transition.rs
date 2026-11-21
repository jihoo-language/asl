mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn next_page_corrects_lifts_and_quotients_by_incoming_differentials() {
    assert_eq!(
        stdout(example("may-survival")),
        "May(3, 2, 6, 0){g0}\nMay(3, 3, 6, 3){g0 + g1}\nMay(4, 3, 6, 3){0}\n"
    );
    let program = include_str!("../examples/may-survival.asl");
    let error = diagnostic(interpret(&format!("{program}\nshow mayNext(b3);")));
    assert!(error.contains("cycle for d_3"), "{error}");
    let error = diagnostic(interpret(&format!(
        "{program}\nshow mayClass[3, 0](mayRep(b1));"
    )));
    assert!(error.contains("May page cycle condition"), "{error}");
}

#[test]
fn genuine_seventh_differential_source_can_be_followed_across_pages() {
    let program = r"
        let e7 = maygen[7, 4, 12, 0, 0];
        let e6 = mayPrev(e7);
        let e5 = mayPrev(e6);
        let e4 = mayPrev(e5);
        assert mayNext(e4) == e5;
        assert mayNext(e5) == e6;
        assert mayNext(e6) == e7;
        let image = mayD(e7);
        assert image == maygen[7, 5, 12, 7, 0];
        assert mayNext(image) == zero[May(8, 5, 12, 7)];
    ";
    assert!(stdout(interpret(program)).is_empty());
    let error = diagnostic(interpret(&format!("{program}\nshow mayNext(e7);")));
    assert!(error.contains("cycle for d_7"), "{error}");
}

#[test]
fn transitions_preserve_linear_functions_and_products_in_both_families() {
    let program = r"
        let forward = \x : May(1, 1, 2, 1) => mayNext(x);
        let backward = \y : May(2, 1, 2, 1) => mayPrev(y);
        let h = maygen[1, 1, 2, 1, 0];
        let h2 = forward(h);
        assert h2 == maygen[2, 1, 2, 1, 0];
        assert forward(backward(h2)) == h2;
        assert forward(backward(zero[May(2, 1, 2, 1)])) == zero[May(2, 1, 2, 1)];
        assert mayNext(mayMul(h, h)) == mayMul(h2, h2);

        let mapH = mayClass[1, 1](mword[One, 0, 0, One, 0, 0](mayRep(h)));
        let mapH2 = mayNext(mapH);
        assert mayRep(mapH2) == mword[One, 0, 0, One, 0, 0](mayRep(h2));
        assert mayNext(mayPrev(mapH2)) == mapH2;
        assert mayNext(ecomp(mapH, mapH)) == ecomp(mapH2, mapH2);
        assert mayRep(mayNext(mayD(mayClass[1, 0](mword[One, 0, 0, One, 0, 0](cword{[0, 1]})))))
            == zero[MapCobar(One, One, 2, 3)];

        let word = cword{[0, 1], [0, 1]};
        let b = mayNext(mayNext(mayClass[1, 0](word)));
        let mapB = mayNext(mayNext(mayClass[1, 0](mword[One, 0, 0, One, 0, 0](word))));
        assert mayRep(mapB) == mword[One, 0, 0, One, 0, 0](mayRep(b));
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn bivariant_page_zero_detects_the_coaction_and_preserves_signed_grades() {
    let program = format!(
        "{}\n{}",
        include_str!("../examples/projective-plane.asl"),
        r"
        let rp = mapmaygen[One, RP2, 0, 0, 2, 0, 0];
        let wedge = mapmaygen[One, Wedge, 0, 0, 2, 0, 0];
        assert mayNext(wedge) == mapmaygen[One, Wedge, 1, 0, 2, 0, 0];
        assert mayNext(mayD(rp)) == zero[MapMay(One, RP2, 1, 1, 2, 0)];
        let other = mapmaygen[Wedge, One, 0, 0, -1, 0, 0];
        let earlier = mayPrev(mayNext(other));
        assert earlier == other;
        "
    );
    assert!(stdout(interpret(&program)).ends_with("Stable(RP2, 1000000, 1000002){g0}\n"));
    let error = diagnostic(interpret(&format!("{program}\nshow mayNext(rp);")));
    assert!(error.contains("cycle for d_0"), "{error}");
}

#[test]
fn page_transitions_preserve_finite_profiles_and_commute_with_restriction() {
    let program = r"
        let full = mayClass[1, 1](cword{[2]});
        let finite = resA[1](full);
        let next = mayNext(finite);
        assert next == mayangen[1, 2, 1, 2, 1, 0];
        assert next == resA[1](mayNext(full));
        assert mayNext(mayPrev(next)) == next;
        assert mayNext(resA[0](finite)) == zero[MayAn(0, 2, 1, 2, 1)];
        let killed = mayD(mayClass[1, 0](resA[1](cword{[0, 1]})));
        assert mayNext(killed) == zero[MayAn(1, 2, 2, 3, 1)];

        let periodic = resA[1](maygen[7, 4, 12, 0, 0]);
        let e8 = mayNext(periodic);
        assert e8 == mayangen[1, 8, 4, 12, 0, 0];
        assert mayNext(mayPrev(e8)) == e8;
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn page_families_and_page_number_limits_are_checked() {
    for (expression, reason) in [
        ("mayPrev(zero[May(0, 1, 1, 0)])".into(), "greater than zero"),
        (
            "mayPrev(zero[MapMay(One, One, 0, 1, 1, 0)])".into(),
            "greater than zero",
        ),
        ("mayNext(sq[0])".into(), "a May or MapMay class"),
        (
            format!("mayNext(zero[May({}, 1, 1, 0)])", usize::MAX),
            "page is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }

    let last = usize::MAX;
    let earlier = last - 1;
    let min = i64::MIN;
    let program = format!(
        "assert mayNext(maygen[{earlier}, 1, 1, 0, 0]) == maygen[{last}, 1, 1, 0, 0];
         assert mayPrev(maygen[{last}, 1, 1, 0, 0]) == maygen[{earlier}, 1, 1, 0, 0];
         assert mayNext(zero[MapMay(One, One, {earlier}, 1, {min}, 0)])
             == zero[MapMay(One, One, {last}, 1, {min}, 0)];"
    );
    assert!(stdout(interpret(&program)).is_empty());
}
