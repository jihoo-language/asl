mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn composition_and_the_differential_act_on_certificates() {
    let output = stdout(example("map-cobar"));
    assert_eq!(
        output,
        concat!(
            "MapCobar(One, One, 2, 3){mword[One, 0, 0, One, 0, 0](cword{[2], [1]})}\n",
            "MapCobar(Circle, One, 0, -1){mword[Circle, 1, 0, One, 0, 0](cword{})}\n"
        )
    );
}

#[test]
fn raw_morphisms_are_discrete_first_class_values() {
    let program = r"
        space Two { a, b };
        let u = mword[Two, 0, 0, Two, 0, 1](cword{});
        let v = mword[Two, 0, 1, Two, 0, 0](cword{});
        let swap = add(u, v);
        let id = \x : MapCobar(Two, Two, 0, 0) => x;
        let translate = \x : MapCobar(Two, Two, 0, 0) => add(x, swap);
        assert id(swap) == swap;
        assert translate(translate(u)) == u;
        assert ecomp(swap, swap) == eRep(eid[Two]);
        assert eClass(eRep(eClass(swap))) == eClass(swap);
        assert eClass(eRep(mapextgen[One, One, 1, 2, 0])) == mapextgen[One, One, 1, 2, 0];
        assert mword[One, 0, 0, One, 0, 0](add(cword{[0, 1]}, cword{[3]})) ==
            add(mword[One, 0, 0, One, 0, 0](cword{[0, 1]}), mword[One, 0, 0, One, 0, 0](cword{[3]}));
        card (\x : MapCobar(Two, Two, 0, 0) => S.hi);
        homotopy (\x : One => u) ~ (\x : One => v);
        betti MapCobar(One, One, 1, 3);
        betti MapCobar(One, One, 0, -1);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "card = 16 / 16\nhomotopy = no\n",
            "betti(MapCobar(One, One, 1, 3)) = [4]\n",
            "betti(MapCobar(One, One, 0, -1)) = [1]\n",
        )
    );
}

#[test]
fn raw_types_check_indices_endpoints_cycles_and_overflow() {
    for (expression, reason) in [
        ("mword[One, 0, 1, One, 0, 0](cword{})", "rank 1"),
        ("mword[One, 0, 0, One, 1, 0](cword{})", "rank 0"),
        ("mword[One, 0, 0, One, 0, 0](sq[0])", "cobar cochain"),
        (
            "eClass(mword[One, 0, 0, One, 0, 0](cword{[0, 1]}))",
            "closed",
        ),
        ("ecomp(eRep(eid[S]), eRep(eid[One]))", "inner target"),
        ("ecomp(eid[One], eRep(eid[One]))", "MapExt"),
        ("eD(eid[One])", "MapCobar"),
        ("eRep(eRep(eid[One]))", "MapExt"),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }
    for expression in [
        format!("eD(zero[MapCobar(One, One, {}, 0)])", usize::MAX),
        format!(
            "ecomp(zero[MapCobar(One, One, 0, {})], zero[MapCobar(One, One, 0, 1)])",
            i64::MAX
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains("out of range"), "{expression}\n{error}");
    }
}

#[test]
fn elementary_morphisms_and_zero_do_not_enumerate_the_ambient_grade() {
    let n = 1usize << (usize::BITS - 3);
    let program = format!(
        "let u = mword[One, 0, 0, One, 0, 0](cword{{[{n}]}});
         let z = zero[MapCobar(One, One, 1, {n})];
         assert eD(u) == zero[MapCobar(One, One, 2, {n})];
         assert ecomp(u, eRep(eid[One])) == u;
         assert eClass(z) == zero[MapExt(One, One, 1, {n})];
         assert eRep(eClass(z)) == z;
         show u;"
    );
    assert_eq!(
        stdout(interpret(&program)),
        format!("MapCobar(One, One, 1, {n}){{mword[One, 0, 0, One, 0, 0](cword{{[{n}]}})}}\n")
    );
    let error = diagnostic(interpret(&format!("betti MapCobar(One, One, 1, {n});")));
    assert!(error.contains("limit"), "{error}");
}
