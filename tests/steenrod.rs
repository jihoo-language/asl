mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn cup_square_and_its_suspended_operation_distinguish_the_projective_plane() {
    assert_eq!(
        stdout(example("projective-plane")),
        concat!(
            "betti(RP2) = [1, 1, 1]\nbetti(Wedge) = [1, 1, 1]\n",
            "H(RP2, 2){g0}\nH(Wedge, 2){0}\n",
            "betti(Susp(RP2)) = [1, 0, 1, 1]\nbetti(Susp(Wedge)) = [1, 0, 1, 1]\n",
            "H(Susp(RP2), 3){g0}\nH(Susp(Wedge), 3){0}\n",
            "Stable(RP2, 1000000, 1000002){g0}\n"
        )
    );
}

#[test]
fn adem_reduces_low_and_high_degree_operations_to_admissible_words() {
    assert_eq!(
        stdout(example("adem")),
        concat!(
            "A(5){Sq^5 + Sq^4 Sq^1}\n",
            "A(256){Sq^255 Sq^1 + Sq^254 Sq^2 + Sq^252 Sq^4 + Sq^248 Sq^8 + ",
            "Sq^240 Sq^16 + Sq^224 Sq^32 + Sq^192 Sq^64}\n"
        )
    );
}

#[test]
fn operation_values_can_be_passed_and_composed_noncommutatively() {
    let program = r"
        let u = mul(sq[2], sq[1]);
        let v = mul(sq[1], sq[2]);
        assert v == sq[3];
        assert mul(sq[0], u) == mul(u, sq[0]);
        assert mul(sq[0], u) == u;
        let swap = \p : A(1) * A(2) => mul(snd(p), fst(p));
        assert swap((sq[1], sq[2])) == u;
        let translate = \x : A(3) => add(x, u);
        assert translate(translate(v)) == v;
        assert mul(sq[2], add(u, v)) == add(mul(sq[2], u), mul(sq[2], v));
        assert mul(mul(sq[2], sq[2]), sq[2]) == mul(sq[2], mul(sq[2], sq[2]));
        show add(v, u);
        show sq[0];
        show add(v, v);
        betti A(0);
        betti A(3);
        betti A(7);
        homotopy (\p : One => u) ~ (\p : One => v);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "A(3){Sq^3 + Sq^2 Sq^1}\nA(0){1}\nA(3){0}\n",
            "betti(A(0)) = [2]\nbetti(A(3)) = [4]\nbetti(A(7)) = [16]\nhomotopy = no\n"
        )
    );
}

#[test]
fn geometric_suspension_has_two_poles_and_its_equator_is_continuous() {
    let program = r"
        space Void {};
        betti Susp(Void);
        betti Susp(Susp(Void));
        betti Susp(Susp(Susp(Void)));
        betti Susp(A(3));
        show north[S];
        show south[S];
        show equator[S](S.hi);
        homotopy (\u : One => north[Void]) ~ (\u : One => south[Void]);
        homotopy (\u : One => north[S]) ~ (\u : One => south[S]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "betti(Susp(Void)) = [2]\nbetti(Susp(Susp(Void))) = [1, 1]\n",
            "betti(Susp(Susp(Susp(Void)))) = [1, 0, 1]\nbetti(Susp(A(3))) = [1, 3]\n",
            "north[S]\nsouth[S]\nequator[S](S.hi)\nhomotopy = no\nhomotopy = yes (2 steps)\n"
        )
    );
}

#[test]
fn degree_zero_suspension_and_shift_with_a_disjoint_basepoint_are_distinct() {
    let program = r"
        space D { a, b };
        let a = cocycle[D, 0] { a };
        let b = cocycle[D, 0] { b };
        let one = add(a, b);
        assert susp(a) == susp(b);
        assert susp(a) == gen[Susp(D), 1, 0];
        assert susp(one) == zero[Susp(D), 1];
        assert pull(equator[D], susp(a)) == zero[D, 1];
        assert susp(susp(a)) == gen[Susp(Susp(D)), 2, 0];
        assert add(shift[1](a), shift[1](b)) == shift[1](one);
        betti Stable(D, 1, 1);
        show shift[1](one);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(Stable(D, 1, 1)) = [4]\nStable(D, 1, 1){g0 + g1}\n"
    );
}

#[test]
fn stable_action_uses_the_base_degree_and_respects_iterated_shifts() {
    let program = r"
        let x = gen[One, 0, 0];
        let w = shift[1000000](x);
        assert shift[999999](shift[1](x)) == w;
        assert shift[0](w) == w;
        assert act(sq[0], w) == w;
        assert act(sq[1], w) == zero[Stable(One, 1000000, 1000001)];
        let z = zero[Stable(One, 5, 3)];
        assert act(sq[2], z) == zero[Stable(One, 5, 5)];
        assert shift[2](z) == zero[Stable(One, 7, 5)];
        assert act(zero[A(0)], w) == zero[Stable(One, 1000000, 1000000)];
        betti Stable(One, 5, 3);
        betti H(Stable(One, 5, 5), 0);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(Stable(One, 5, 3)) = [1]\nbetti(H(Stable(One, 5, 5), 0)) = [4]\n"
    );
}

#[test]
fn nontop_square_satisfies_cartan_on_a_product_of_finite_spaces() {
    let calculations = r"
        space C { p < a, q < a, p < b, q < b };
        let left = \z : RP2 * C => fst(z);
        let right = \z : RP2 * C => snd(z);
        let a = pull(left, alpha);
        let b = pull(right, gen[C, 1, 0]);
        let ab = cup(a, b);
        let resultCartan = act(sq[1], ab);
        assert resultCartan == add(cup(act(sq[1], a), b), cup(a, act(sq[1], b)));
        assert resultCartan == gen[RP2 * C, 3, 0];
        assert act(sq[1], a) == pull(left, act(sq[1], alpha));
        assert act(mul(sq[1], sq[1]), ab) == act(sq[1], act(sq[1], ab));
        show resultCartan;
    ";
    let program = format!(
        "{}\n{calculations}",
        include_str!("../examples/projective-plane.asl")
    );
    assert!(stdout(interpret(&program)).ends_with("H(RP2 * C, 3){g0}\n"));
}

#[test]
fn general_cup_i_detects_sq_two_on_a_named_model_of_suspended_rp4() {
    // A face of the five-dimensional cross-polytope selects signed axes.
    // Identify antipodal faces by making the first selected sign positive.
    fn face(support: u32, negatives: u32) -> String {
        let first_axis = 1 << support.trailing_zeros();
        let mut signs = negatives & support;
        if signs & first_axis != 0 {
            signs ^= support;
        }
        format!("f{support}_{signs}")
    }
    let mut relations = Vec::new();
    for support in 1u32..32 {
        let first_axis = 1 << support.trailing_zeros();
        for negatives in 0u32..32 {
            if negatives & !support != 0 || negatives & first_axis != 0 {
                continue;
            }
            let upper = face(support, negatives);
            if support.count_ones() > 1 {
                for axis in (0..5).map(|i| 1 << i).filter(|bit| support & bit != 0) {
                    let lower = face(support ^ axis, negatives);
                    relations.push(format!("{lower} < {upper}"));
                }
            }
            if support == 31 {
                relations.push(format!("{upper} < north, {upper} < south"));
            }
        }
    }
    // The orbit face poset has order complex a subdivision of RP⁴.
    // The two extra maxima suspend it. If t generates H¹(RP⁴),
    // Sq²(σt²) = σt⁴, and Sq²Sq¹ and Sq¹Sq² act differently on σt.
    let calculation = "
        let t = gen[P, 2, 0];
        let u = gen[P, 3, 0];
        let v = gen[P, 5, 0];
        assert act(sq[1], t) == u;
        assert act(sq[2], u) == v;
        assert act(mul(sq[2], sq[1]), t) == v;
        assert act(mul(sq[1], sq[2]), t) == zero[P, 5];";
    let program = format!("space P {{ {} }};\n{calculation}", relations.join(", "));
    assert!(stdout(interpret(&program)).is_empty());
}

#[test]
fn formal_operations_and_stable_classes_keep_their_grades_and_spaces() {
    for expression in [
        "add(sq[1], sq[2])",
        "mul(sq[1], zero[S, 0])",
        "act(zero[S, 0], sq[0])",
        "act(sq[1], One.pt)",
        "susp(sq[1])",
        "shift[2](S.lo)",
        "equator[One](S.lo)",
        "add(shift[2](gen[One, 0, 0]), zero[Stable(S, 2, 2)])",
        "add(shift[2](gen[One, 0, 0]), zero[Stable(One, 1, 2)])",
        "add(shift[0](gen[One, 0, 0]), gen[One, 0, 0])",
    ] {
        let error = diagnostic(interpret(&format!("\nshow {expression};")));
        assert!(error.contains(":2:"), "{expression}\n{error}");
    }
    for expression in [
        format!("mul(sq[{}], sq[1])", usize::MAX),
        format!("act(sq[{}], zero[One, 1])", usize::MAX),
        format!("shift[1](shift[{}](gen[One, 0, 0]))", usize::MAX),
    ] {
        let error = diagnostic(interpret(&format!("\nshow {expression};")));
        assert!(error.contains(":2:"), "{expression}\n{error}");
    }
    assert_eq!(
        stdout(interpret("show sq[1000000];")),
        "A(1000000){Sq^1000000}\n"
    );
    assert!(diagnostic(interpret("betti A(1000000);")).contains("limit"));
}

#[test]
fn earlier_bindings_can_use_the_new_contextual_names() {
    let program = r"
        space A { p };
        space Susp { p };
        space Stable { p };
        let susp = \x : A => x;
        let shift = susp;
        let act = shift;
        let mul = act;
        let sq = mul;
        assert sq(mul(act(shift(susp(A.p))))) == A.p;
        assert (\x : Susp => x)(Susp.p) == Susp.p;
        assert (\x : Stable => x)(Stable.p) == Stable.p;
        show sq[1];
    ";
    assert_eq!(stdout(interpret(program)), "A(1){Sq^1}\n");
}
