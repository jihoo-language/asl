mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn may_differentials_detect_the_first_and_third_cobar_filtration_jumps() {
    assert_eq!(
        stdout(example("may")),
        concat!(
            "May(1, 2, 3, 1){g0}\n",
            "May(3, 3, 6, 3){g0 + g1}\n",
            "maydim[3, 2, 6, 0] = 1\n",
            "maydim[4, 2, 6, 0] = 0\n",
            "maydim[3, 3, 6, 3] = 2\n",
            "maydim[4, 3, 6, 3] = 1\n",
        )
    );
}

#[test]
fn associated_graded_pages_and_later_cycles_enforce_their_conditions() {
    let program = r"
        let xi1cube = cword{[3]};
        let column = add(cword{[1], [2]}, cword{[2], [1]});
        let x = mayClass[0, 1](xi1cube);
        let dx = mayD(x);
        assert dx == mayClass[0, 1](column);
        assert dx == add(maygen[0, 2, 3, 1, 0], maygen[0, 2, 3, 1, 1]);
        assert mayD(dx) == zero[May(0, 3, 3, 1)];
        assert mayClass[1, 1](column) == zero[May(1, 2, 3, 1)];
        assert mayClass[0, 0](xi1cube) == zero[May(0, 1, 3, 0)];
    ";
    assert_eq!(stdout(interpret(program)), "");

    let outside = diagnostic(interpret("show mayClass[0, 1](cword{[0, 1]});"));
    assert!(outside.contains("below the requested May filtration"));
    let uncorrected = diagnostic(interpret("show mayClass[3, 0](cword{[0, 1], [0, 1]});"));
    assert!(uncorrected.contains("May page cycle condition"));
    let deeper = diagnostic(interpret("show mayClass[3, 0](cword{[3]});"));
    assert!(deeper.contains("May page cycle condition"));
}

#[test]
fn may_spaces_support_finite_functions_and_large_pages() {
    let program = r"
        let unit = mayClass[30, 0](cword{});
        let h0 = maygen[30, 1, 1, 0, 0];
        assert mayMul(unit, h0) == h0;
        assert mayMul(h0, unit) == h0;
        assert mayD(unit) == zero[May(30, 1, 0, 30)];
        assert mayD(h0) == zero[May(30, 2, 1, 30)];
        let translate = \x : May(30, 1, 1, 0) => add(x, h0);
        assert translate(h0) == zero[May(30, 1, 1, 0)];
        assert translate(zero[May(30, 1, 1, 0)]) == h0;
        betti May(30, 1, 1, 0);
        maydim[30, 3, 6, 0];
        maydim[30, 3, 6, 1];
        maydim[30, 3, 6, 2];
        maydim[30, 3, 6, 3];
        extdim[3, 6];
        maydim[30, 0, 1, 0];
        maydim[30, 2, 1, 0];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "betti(May(30, 1, 1, 0)) = [2]\n",
            "maydim[30, 3, 6, 0] = 0\n",
            "maydim[30, 3, 6, 1] = 0\n",
            "maydim[30, 3, 6, 2] = 0\n",
            "maydim[30, 3, 6, 3] = 1\n",
            "extdim[3, 6] = 1\n",
            "maydim[30, 0, 1, 0] = 0\n",
            "maydim[30, 2, 1, 0] = 0\n",
        )
    );
}

#[test]
fn may_computes_a_nonzero_seventh_differential_from_a_cobar_lift() {
    assert_eq!(
        stdout(example("may-seventh")),
        "May(7, 5, 12, 7){g0}\nmaydim[8, 4, 12, 0] = 0\nmaydim[8, 5, 12, 7] = 0\n"
    );
}

#[test]
fn may_checks_pages_and_bounds_filtration_arithmetic_and_basis_work() {
    let mismatch = diagnostic(interpret(
        "show mayMul(zero[May(0, 1, 1, 0)], zero[May(1, 1, 1, 0)]);",
    ));
    assert!(mismatch.contains("mayMul requires classes on the same page"));
    let index = diagnostic(interpret("show maygen[0, 1, 1, 0, 1];"));
    assert!(index.contains("generator index 1 is out of range"));
    let budget = diagnostic(interpret("maydim[0, 1, 1000000, 0];"));
    assert!(budget.contains("cobar basis exceeds the enumeration limit"));

    let max = usize::MAX;
    assert_eq!(
        stdout(interpret(&format!("maydim[{max}, 3, 6, 3];"))),
        format!("maydim[{max}, 3, 6, 3] = 1\n")
    );
    let overflow = diagnostic(interpret(&format!("show mayD(zero[May({max}, 1, 2, 1)]);")));
    assert!(overflow.contains("mayD degree is out of range"));
}
