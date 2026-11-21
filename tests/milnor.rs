mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn milnor_operations_respect_coproduct_orientation_and_admissible_reduction() {
    let example = stdout(example("milnor-duality"));
    assert!(example.starts_with("A(3){Sq^3 + Sq^2 Sq^1}\nS.hi\nS.lo\n"));
    assert!(example.ends_with("S.hi\n"));

    let program = r"
        let a = cword{[3]};
        let b = cword{[0, 1]};
        let f = \theta : A(3) => milnorPair(theta, b);
        let x = mul(sq[2], sq[1]);
        let y = mul(sq[1], sq[2]);
        assert f(x) == S.hi;
        assert f(y) == S.lo;
        assert milnorPair(x, a) == S.hi;
        assert milnorPair(y, a) == S.hi;
        assert milnorPair(x, add(a, b)) == S.lo;
        assert milnorPair(add(x, y), add(a, b)) == S.hi;
        assert milnor{[0, 1, 0, 0]} == add(x, y);
        assert milnor{[0, 0]} == sq[0];
        assert milnor{[]} == sq[0];
        assert milnorPair(sq[0], cword{}) == S.hi;
        assert milnorPair(zero[A(0)], cword{}) == S.lo;
        assert milnorPair(sq[0], zero[Cobar(1, 0)]) == S.lo;
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn milnor_literals_are_exact_dual_bases_in_many_degrees() {
    for degree in 1..=15 {
        let mut monomials = Vec::new();
        for c in 0..=degree / 15 {
            for b in 0..=(degree - 15 * c) / 7 {
                for a in 0..=(degree - 15 * c - 7 * b) / 3 {
                    let first = degree - 3 * a - 7 * b - 15 * c;
                    monomials.push(format!("[{first}, {a}, {b}, {c}]"));
                }
            }
        }
        let mut program = String::new();
        for (row, operation) in monomials.iter().enumerate() {
            for (column, word) in monomials.iter().enumerate() {
                let point = if row == column { "hi" } else { "lo" };
                program.push_str(&format!(
                    "assert milnorPair(milnor{{{operation}}}, cword{{{word}}}) == S.{point};\n"
                ));
            }
        }
        assert!(stdout(interpret(&program)).is_empty(), "degree {degree}",);
    }
}

#[test]
fn coaction_transposes_the_action_on_cohomology_and_has_a_counit() {
    let calculations = r"
        let onRP = \u : H(RP2, 1) => coact[2, 0](u);
        let onWedge = \u : H(Wedge, 1) => coact[2, 0](u);
        assert onRP(alpha) == cword{[1]};
        assert onWedge(beta) == zero[Cobar(1, 1)];
        assert onRP(zero[RP2, 1]) == zero[Cobar(1, 1)];
        assert coact[3, 0](susp(alpha)) == cword{[1]};
        assert coact[3, 0](susp(beta)) == zero[Cobar(1, 1)];
        assert coact[2, 0](gen[RP2, 0, 0]) == zero[Cobar(1, 2)];
        assert coact[1, 0](alpha) == cword{};
        assert coact[1, 0](zero[RP2, 1]) == zero[Cobar(0, 0)];
        space D { a, b };
        assert coact[0, 1](gen[D, 0, 0]) == zero[Cobar(0, 0)];
        assert coact[0, 1](gen[D, 0, 1]) == cword{};
    ";
    let program = format!(
        "{}\n{calculations}",
        include_str!("../examples/projective-plane.asl")
    );
    assert!(stdout(interpret(&program)).ends_with("Stable(RP2, 1000000, 1000002){g0}\n"));
}

#[test]
fn coaction_detects_the_second_milnor_monomial_on_projective_four_space() {
    fn face(axes: u32, mut signs: u32) -> String {
        let first = 1 << axes.trailing_zeros();
        if signs & first != 0 {
            signs ^= axes;
        }
        format!("f{axes}_{signs}")
    }

    // The antipodal face poset of the five-dimensional cross-polytope
    // has order complex a subdivision of RP^4.
    let mut relations = Vec::new();
    for axes in 1u32..32 {
        let first = 1 << axes.trailing_zeros();
        for signs in 0..32 {
            if signs & !axes != 0 || signs & first != 0 {
                continue;
            }
            for axis in 0..5 {
                let lower = axes ^ (1 << axis);
                if axes & (1 << axis) != 0 && lower != 0 {
                    relations.push(format!(
                        "{} < {}",
                        face(lower, signs & lower),
                        face(axes, signs),
                    ));
                }
            }
        }
    }
    let program = format!(
        "space RP4 {{ {} }};
         let x = gen[RP4, 1, 0];
         assert act(sq[3], x) == zero[RP4, 4];
         assert act(milnor{{[0, 1]}}, x) == gen[RP4, 4, 0];
         assert coact[4, 0](x) == cword{{[0, 1]}};",
        relations.join(", ")
    );
    assert!(stdout(interpret(&program)).is_empty());
}

#[test]
fn milnor_ranges_are_checked_and_sparse_pairing_does_not_enumerate() {
    for (expression, reason) in [
        ("milnorPair(sq[3], cword{[1]})", "Cobar(1, 3)"),
        ("milnorPair(sq[3], cword{[1], [2]})", "Cobar(1, 3)"),
        ("milnorPair(cword{}, cword{})", "Steenrod operation"),
        ("coact[1, 0](zero[One, 2])", "below input degree"),
        ("coact[0, 1](zero[One, 0])", "rank 1"),
        ("coact[1, 0](zero[One, 0])", "rank 0"),
        ("coact[0, 0](shift[0](gen[One, 0, 0]))", "cohomology class"),
        ("milnor{[0, 1000000]}", "enumeration limit"),
    ] {
        let stderr = diagnostic(interpret(&format!("show {expression};")));
        assert!(stderr.contains(reason), "{expression}\n{stderr}");
    }
    let max = usize::MAX;
    let overflow = diagnostic(interpret(&format!("show milnor{{[{max}, 1]}};")));
    assert!(overflow.contains("milnor degree is out of range"));
    let high = 1usize << (usize::BITS - 1);
    let total = high + 1;
    assert!(stdout(interpret(&format!(
        "assert milnorPair(milnor{{[{max}]}}, cword{{[{max}]}}) == S.hi;\n
             assert milnorPair(mul(sq[{high}], sq[1]), cword{{[{total}]}}) == S.hi;"
    )))
    .is_empty());
}
