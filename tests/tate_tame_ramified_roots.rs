mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn tame_quadratic_and_cubic_coordinates_realize_the_galois_relation() {
    let capture = stdout(example("tate-tame-ramified-roots"));
    let lines: Vec<_> = capture.lines().collect();
    assert_eq!(lines.len(), 4);
    assert!(lines[0].contains("slope=(1, 2, 2), ring=Q, normalization=2, y_precision=6, y_moduli=[27, 27], root_precision=7, root_moduli=[81, 27], verified_mod_pi=8"), "{capture}");
    assert!(
        lines[0].contains("roots=[[[0], [26]], [[0], [1]]]"),
        "{capture}"
    );
    assert!(
        lines[0].ends_with("frobenius=[0, 1], zeta=[80], inertia=[1, 0], relation=verified"),
        "{capture}"
    );
    assert!(lines[1].contains("slope=(1, 3, 3), ring=U, normalization=3, y_precision=9, y_moduli=[125, 125, 125], root_precision=10, root_moduli=[625, 125, 125], verified_mod_pi=12"), "{capture}");
    assert!(lines[1].contains("roots=[[[0, 0], [0, 111], [0, 0]], [[0, 0], [14, 14], [0, 0]], [[0, 0], [111, 0], [0, 0]]]"), "{capture}");
    assert!(
        lines[1]
            .ends_with("frobenius=[1, 0, 2], zeta=[0, 1], inertia=[1, 2, 0], relation=verified"),
        "{capture}"
    );
    assert_eq!(&lines[2..], ["rvTameVal = 1", "rvTameVal = at_least=12"]);

    for (modulus, expected) in [(5, 1), (25, 11), (125, 111), (625, 486)] {
        let roots: Vec<_> = (0_i64..modulus)
            .filter(|&a| (a.pow(3) - 6).rem_euclid(modulus) == 0)
            .collect();
        assert_eq!(roots, [expected]);
    }
}

#[test]
fn mixed_newton_slopes_and_the_original_amice_chart() {
    let source = "
        rvunram U[5,4]{1,1,1};
        let mixed = rvTatePoly[5,4](150,0,-30,-5,0,1);
        rvTateTame[U,0](mixed); rvTateTame[U,1](mixed);
        rvTateTame[U,0](rvTatePoly[5,4](-150,0,0,1));

        rvheckepacket Delta(rvHeckeT[12,2]);
        rvheckeeigenperiod P[Delta,plus,2];
        rvmazurordinary O[P,11,4];
        rvunram Q[11,2]{0,1};
        let d = rvMazurDisc[0,1](rvMazurAmice[O,3,0,0]);
        rvTateTame[Q,0](d);
        rvTateTame[Q,0](rvMazurDiscRestrict[5,1](d));";
    let capture = stdout(interpret(source));
    let lines: Vec<_> = capture.lines().collect();
    assert!(
        lines[0].contains("slope=(1, 3, 3), ring=U, normalization=5, y_precision=7"),
        "{capture}"
    );
    assert!(lines[0].contains("removed_zero_power=2"), "{capture}");
    assert!(
        lines[0]
            .ends_with("frobenius=[1, 0, 2], zeta=[0, 1], inertia=[1, 2, 0], relation=verified"),
        "{capture}"
    );
    assert!(
        lines[1].contains("slope=(1, 2, 2), ring=U, normalization=4, y_precision=4"),
        "{capture}"
    );
    assert!(
        lines[1].contains("roots=[[[0, 0], [24, 0]], [[0, 0], [1, 0]]]"),
        "{capture}"
    );
    assert!(
        lines[2].contains("slope=(2, 3, 3), ring=U, normalization=6, y_precision=6"),
        "{capture}"
    );
    assert!(
        lines[2]
            .ends_with("frobenius=[1, 0, 2], zeta=[0, 1], inertia=[2, 0, 1], relation=verified"),
        "{capture}"
    );
    assert!(
        lines[3].contains("coordinate=RVMazurDisc(O,3,0,0,0,1)"),
        "{capture}"
    );
    assert!(lines[3].contains("roots=[[[115]]]"), "{capture}");
    assert!(
        lines[4].contains("coordinate=RVMazurDisc(O,3,0,0,55,2)"),
        "{capture}"
    );
    assert!(lines[4].contains("roots=[[[10]]]"), "{capture}");
}

#[test]
fn repeated_wild_outside_and_precision_deficits_remain_distinct() {
    let source = "
        rvunram Q[3,4]{0,1};
        rvunram A[5,4]{0,1};
        rvunram B[5,1]{1,1,1};
        rvTateTame[Q,0](rvTatePoly[3,4](9,0,-6,0,1));
        rvTateTame[Q,0](rvTatePoly[3,4](-3,0,0,1));
        rvTateTame[Q,0](rvTatePoly[3,4](0,0,1));
        let cubic = rvTatePoly[5,4](-30,0,0,1);
        rvTateTame[A,0](cubic); rvTateTame[B,0](cubic);
        rvTateTame[Q,0](rvTatePoly[3,4](0));";
    let capture = stdout(interpret(source));
    let lines: Vec<_> = capture.lines().collect();
    assert!(
        lines[0]
            .contains("roots=[], repeated=[([[1], [1]], 2), ([[2], [1]], 2)], outside_field=[]"),
        "{capture}"
    );
    assert!(lines[1].ends_with("wild, unresolved"), "{capture}");
    assert!(
        lines[2].ends_with("valuation_floor=(2, 1, 2), unresolved"),
        "{capture}"
    );
    assert!(
        lines[3].contains(
            "roots=[[[0], [111], [0]]], repeated=[], outside_field=[([[1], [1], [1]], 1)]"
        ),
        "{capture}"
    );
    assert!(lines[3].ends_with("inertia=unavailable"), "{capture}");
    assert!(
        lines[4].ends_with(
            "normalization=3, available_uniformizer_precision=3, insufficient_precision"
        ),
        "{capture}"
    );
    assert!(lines[5].ends_with("precision=4, unresolved"), "{capture}");
}

#[test]
fn typed_tame_arithmetic_respects_each_row_precision_and_base_identity() {
    let output = stdout(interpret(
        "
        rvunram U[5,4]{1,1,1};
        let p = rvTameElem[U,3,5]([],[1]);
        let a = rvTameElem[U,3,5]([26,27],[26,27],[6,7]);
        assert a == rvTameElem[U,3,5]([1,2],[1,2],[1,2]);
        assert add(a,neg(a)) == zero[RVTame(U,3,5)];
        let carry = rvTameMul(rvTameMul(p,p),p);
        show carry;
        rvTameVal(carry);
        assert rvTameFrob[2](a) == a;
        assert rvTameFrob[1](rvTameInertia(a)) ==
            rvTameInertia(rvTameInertia(rvTameFrob[1](a)));",
    ));
    assert_eq!(
        output,
        "RVTame(U,3,5){[[5, 0], [0, 0], [0, 0]]}\nrvTameVal = 3\n"
    );

    let header = "rvunram U[5,4]{1,1,1}; rvunram V[5,4]{1,1,1}; rvunram A[5,4]{0,1};";
    for (source, expected) in [
        ("show rvTameElem[U,5,8]([1]);", "is wild"),
        ("show rvTameElem[U,3,13]([1]);", "uniformizer precision"),
        ("show rvTameElem[U,2,2]([0],[0],[0]);", "more rows"),
        (
            "show rvTameEval(rvTameElem[U,3,12]([1]),rvTatePoly[5,3](1));",
            "insufficient",
        ),
        (
            "show rvTameEval(rvTameElem[U,3,4]([1]),rvTatePoly[3,2](1));",
            "primes do not match",
        ),
        (
            "show rvTameTeich[U,3,4](rvUnramElem[V]{0,1});",
            "RVUnram(U)",
        ),
        (
            "show rvTameInertia(rvTameElem[A,3,4]([0],[1]));",
            "does not contain the required tame roots",
        ),
        (
            "rvTateTame[U,1](rvTatePoly[5,4](-30,0,0,1));",
            "Newton segment index",
        ),
        (
            "rvTateTame[U,0](rvTatePoly[3,4](-3,0,1));",
            "primes do not match",
        ),
        (
            "rvTateTame[Missing,0](rvTatePoly[5,4](0));",
            "unknown unramified ring",
        ),
    ] {
        let result = diagnostic(interpret(&format!("{header} {source}")));
        assert!(result.contains(expected), "{result}");
    }
}

#[test]
fn hensel_digits_use_mixed_rows_and_inertia_can_find_a_new_generator() {
    let capture = stdout(interpret(
        "
        rvunram Q[3,4]{0,1};
        let p = rvTatePoly[3,4](-3,3,1);
        rvTateTame[Q,0](p);
        assert rvTameEval(rvTameElem[Q,2,8]([39],[7]),p) == zero[RVTame(Q,2,8)];
        assert rvTameEval(rvTameElem[Q,2,8]([39],[20]),p) == zero[RVTame(Q,2,8)];
        rvunram W[3,4]{2,1,1};
        rvTateTame[W,0](rvTatePoly[3,4](-3,0,0,0,1));",
    ));
    let lines: Vec<_> = capture.lines().collect();
    assert!(
        lines[0].contains("y=[[[20], [13]], [[7], [13]]], roots=[[[39], [20]], [[39], [7]]]"),
        "{capture}"
    );
    assert!(
        lines[1].contains("slope=(1, 4, 4), ring=W, normalization=4"),
        "{capture}"
    );
    assert!(
        lines[1].ends_with(
            "frobenius=[0, 3, 2, 1], zeta=[56, 31], inertia=[3, 0, 1, 2], relation=verified"
        ),
        "{capture}"
    );

    // In Q the normalized equation is y² + pi*y - 1. Its two rows
    // can be evaluated independently as integer congruences modulo 27.
    let mut solutions = Vec::new();
    for a in 0_i32..27 {
        for b in 0_i32..27 {
            if (a * a + 3 * b * b + 3 * b - 1).rem_euclid(27) == 0
                && (2 * a * b + a).rem_euclid(27) == 0
            {
                solutions.push((a, b));
            }
        }
    }
    assert_eq!(solutions, [(7, 13), (20, 13)]);
}
