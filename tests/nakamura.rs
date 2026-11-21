mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn corrected_squares_construct_the_third_and_seventh_may_transgressions() {
    assert_eq!(
        stdout(example("nakamura-may")),
        concat!(
            "May(3, 3, 6, 3){g0 + g1}\n",
            "May(7, 5, 12, 7){g0}\n",
            "ExtAn(1, 4, 12){g0}\n",
        )
    );
    let source = include_str!("../examples/nakamura-may.asl");
    let continuation = format!("{source}\nshow b; show cobarD(b); show c; show cobarD(c);");
    let output = stdout(interpret(&continuation));
    let sizes: Vec<_> = output
        .lines()
        .skip(3)
        .map(|line| line.matches("cword{").count())
        .collect();
    assert_eq!(sizes, [4, 2, 98, 22]);

    let error = diagnostic(interpret(&format!("{source}\nshow mayClass[7, 0](c);")));
    assert!(error.contains("May page cycle condition"), "{error}");
    let error = diagnostic(interpret(&format!("{source}\nshow mayLiftTo[8, 0](c);")));
    assert!(
        error.contains("nonzero obstruction below filtration 8"),
        "{error}"
    );
}

#[test]
fn raw_square_commutes_with_the_differential_and_agrees_with_square_on_classes() {
    let program = r"
        let u = cword{[0, 1]};
        let v = cword{[3]};
        let a = add(u, v);
        let d = cobarD(a);
        assert cNakamura[0](a) == add(cword{[0, 2]}, cword{[6]});
        assert cNakamura[1](a) == add(concat(a, a), cCup[1](a, d));
        assert cNakamura[2](a) == concat(a, d);
        assert cNakamura[3](a) == zero[Cobar(4, 6)];
        assert cobarD(cNakamura[0](a)) == cNakamura[0](d);
        assert cobarD(cNakamura[1](a)) == cNakamura[1](d);
        assert cobarD(cNakamura[2](a)) == cNakamura[2](d);
        assert cobarD(cNakamura[3](a)) == cNakamura[3](d);

        let x = cword{[1]};
        let h = extClass(x);
        assert extClass(cNakamura[0](x)) == extSq[0](h);
        assert extClass(cNakamura[1](x)) == extSq[1](h);
        assert extClass(cNakamura[2](x)) == extSq[2](h);
        assert cNakamura[0](cword{}) == cword{};
        assert cNakamura[1](cword{}) == zero[Cobar(1, 0)];
        assert cNakamura[0](resA[1](u)) == zero[CobarAn(1, 1, 6)];
        assert cNakamura[1](resA[1](u)) == resA[1](cNakamura[1](u));
        let rawSquare = \x : Cobar(1, 3) => cNakamura[1](x);
        assert rawSquare(a) == cNakamura[1](a);
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn lifts_work_for_scalar_and_bivariant_cochains_and_preserve_valid_input() {
    let program = r"
        let u = cword{[0, 1]};
        let b = cNakamura[1](u);
        let c = cNakamura[2](b);
        let m = mword[One, 0, 0, One, 0, 0](c);
        assert mayLiftTo[7, 0](m) == mword[One, 0, 0, One, 0, 0](mayLiftTo[7, 0](c));
        assert mayLiftTo[6, 0](m) == m;
        assert mayLiftTo[0, 0](u) == u;
        assert mayLiftTo[1, 0](u) == u;
        assert mayLiftTo[2, 0](cword{[3]}) == zero[Cobar(1, 3)];
        space Circle { a < ab, b < ab, b < bc, c < bc, c < ca, a < ca };
        let projection = mword[Circle, 1, 0, One, 0, 0](cword{});
        assert mayLiftTo[12, 0](projection) == projection;
    ";
    assert!(stdout(interpret(program)).is_empty());
    let non_cycle = format!(
        "{}\nshow mayLiftTo[1, 0](mword[One, 0, 0, RP2, 2, 0](cword{{}}));",
        include_str!("../examples/projective-plane.asl")
    );
    let error = diagnostic(interpret(&non_cycle));
    assert!(error.contains("filtration 0 terms fixed"), "{error}");
}

#[test]
fn corrections_and_raw_squares_check_degrees_and_keep_sparse_zero_and_cycles() {
    for (expression, reason) in [
        ("cNakamura[1](sq[1])", "cobar cochain"),
        ("mayLiftTo[1, 0](sq[1])", "cobar cochain"),
        (
            "mayLiftTo[0, 1](cword{[1]})",
            "below requested filtration 1",
        ),
        ("mayLiftTo[2, 0](cword{[0, 1]})", "nonzero obstruction"),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }
    let max = usize::MAX;
    for (expression, reason) in [
        (
            format!("cNakamura[{max}](cword{{[1]}})"),
            "length is out of range",
        ),
        (
            format!("cNakamura[0](zero[Cobar(0, {max})])"),
            "internal degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }

    let high = 1usize << (usize::BITS - 3);
    let source = format!(
        "let u = cword{{[{high}]}};
         assert mayLiftTo[{max}, 0](u) == u;
         assert mayLiftTo[{max}, {max}](zero[Cobar({max}, {max})]) == zero[Cobar({max}, {max})];
         assert cNakamura[0](u) == cword{{[{}]}};",
        2 * high,
    );
    assert!(stdout(interpret(&source)).is_empty());
}
