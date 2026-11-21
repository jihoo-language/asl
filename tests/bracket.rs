mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn a_non_affine_fourfold_locus_and_three_different_torus_loci() {
    assert_eq!(
        stdout(example("massey-locus")),
        "card = 344 / 512\n(S.hi, S.hi)\nS.lo\ncard = 0 / 2\ncard = 1 / 2\ncard = 2 / 2\n"
    );
}

#[test]
fn every_matrix_class_agrees_with_the_determinant_test() {
    let mut source = String::from(
        r"space G { p < a, q < a, p < b, q < b, p < c, q < c, p < d, q < d };
          let l = \v : G * G => fst(v);
          let r = \v : G * G => snd(v);
          let z = zero[G * G, 1];
          let image = bracket(z, z, z, z);
          let m0 = zero[G * G, 2];",
    );
    for i in 0..3 {
        for j in 0..3 {
            writeln!(
                source,
                "let e{} = cup(pull(l, gen[G, 1, {i}]), pull(r, gen[G, 1, {j}]) );",
                3 * i + j
            )
            .unwrap();
        }
    }
    for mask in 0..512u16 {
        if mask != 0 {
            let bit = mask.trailing_zeros();
            let rest = mask ^ (1 << bit);
            writeln!(source, "let m{mask} = add(m{rest}, e{bit});").unwrap();
        }
        // Each term selects one entry of each row and column. Signs vanish in F₂.
        let at = |i: u32, j: u32| (mask >> (3 * i + j)) & 1;
        let mut determinant = 0;
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    determinant ^= at(0, i) & at(1, j) & at(2, 3 - i - j);
                }
            }
        }
        let expected = if determinant == 0 { "hi" } else { "lo" };
        writeln!(source, "assert image(m{mask}) == S.{expected};").unwrap();
    }
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn mixed_degrees_remember_which_tensor_classes_are_decomposable() {
    let source = r"
        space G { p < e, q < e, p < f, q < f, p < g, q < g };
        let l = \v : G * Susp(G) => fst(v);
        let r = \v : G * Susp(G) => snd(v);
        let a = gen[G, 1, 0];
        let b = gen[G, 1, 1];
        let p = cup(pull(l, a), pull(r, susp(a)));
        let q = cup(pull(l, b), pull(r, susp(b)));
        let z = zero[G * Susp(G), 1];
        let w = zero[G * Susp(G), 2];
        let image = bracket(z, z, z, w);
        assert image == bracket(w, z, z, z);
        assert image(p) == S.hi;
        assert image(q) == S.hi;
        assert image(add(p, q)) == S.lo;
        card image;
    ";
    assert_eq!(stdout(interpret(source)), "card = 10 / 16\n");
}

fn moore_space(order: usize) -> String {
    // The mapping cylinder of a degree-order map of circles, coned off at its source.
    let mut source = String::from("space M { o < f, p < f, o < r, p < r,\n");
    for i in 0..2 * order {
        let next = (i + 1) % (2 * order);
        let (vertex, edge) = if i % 2 == 0 { ("o", "f") } else { ("p", "r") };
        writeln!(
            source,
            "v{i} < e{i}, v{next} < e{i}, v{i} < {vertex}, e{i} < {edge}, e{i} < cap,"
        )
        .unwrap();
    }
    source.push_str("}; let t = gen[M, 1, 0]; let z = zero[M, 1];\n");
    source
}

fn matrix_power(generator: &[u8; 6], exponent: usize) -> [u8; 6] {
    let mut power = [1, 2, 4, 8, 16, 32];
    for _ in 0..exponent {
        power = power.map(|row| {
            let mut result = 0;
            for (k, image) in generator.iter().enumerate() {
                if row & (1 << k) != 0 {
                    result ^= image;
                }
            }
            result
        });
    }
    power
}

#[test]
fn all_240_moore_space_brackets_agree_with_independent_unitriangular_lifts() {
    for order in [2, 4, 6, 8] {
        let mut source = moore_space(order);
        let identity = [1, 2, 4, 8, 16, 32];
        for n in 2..=5 {
            let mut entries = Vec::new();
            for i in 0..n {
                for j in i + 1..=n {
                    if (i, j) != (0, n) {
                        entries.push((i, j));
                    }
                }
            }
            let mut attainable = vec![[false; 2]; 1 << n];
            for mask in 0..1 << entries.len() {
                let mut generator = identity;
                let mut inputs = 0;
                for (k, &(i, j)) in entries.iter().enumerate() {
                    if mask & (1 << k) != 0 {
                        generator[i] |= 1 << j;
                        if j == i + 1 {
                            inputs |= 1 << i;
                        }
                    }
                }
                let mut power = matrix_power(&generator, order);
                let corner = usize::from((power[0] >> n) & 1);
                power[0] &= !(1 << n);
                if power == identity {
                    attainable[inputs][corner] = true;
                }
            }
            for (mask, answers) in attainable.iter().enumerate() {
                let inputs: Vec<_> = (0..n)
                    .map(|i| if mask & (1 << i) == 0 { "z" } else { "t" })
                    .collect();
                let name = format!("p{n}_{mask}");
                writeln!(source, "let {name} = bracket({});", inputs.join(", ")).unwrap();
                for (member, class) in answers.iter().zip(["zero[M, 2]", "gen[M, 2, 0]"]) {
                    let expected = if *member { "hi" } else { "lo" };
                    writeln!(source, "assert {name}({class}) == S.{expected};").unwrap();
                }
            }
        }
        assert!(stdout(interpret(&source)).is_empty(), "M(C{order},1)");
    }
}

#[test]
fn zero_rank_target_still_requires_every_proper_obstruction_to_vanish() {
    let source = moore_space(2)
        + "
        let w = zero[M, 2];
        card bracket(t, t, w);
        card bracket(t, t, z, w);
        card bracket(z, t, z, w);
        card bracket(z, z, w);
    ";
    assert_eq!(
        stdout(interpret(&source)),
        "card = 0 / 1\ncard = 0 / 1\ncard = 1 / 1\ncard = 1 / 1\n"
    );
}

#[test]
fn predicates_keep_the_existing_function_semantics_and_names() {
    let program = r"
        space Void {};
        card (\p : Void => S.hi);
        card (\x : S => x);
        let z = zero[One, 1];
        assert bracket(z, z, z) == (\u : H(One, 2) => S.hi);
        let card = \f : One => f;
        let bracket = card;
        assert bracket(One.pt) == One.pt;
        assert bracket(z, z)(zero[One, 2]) == S.hi;
    ";
    assert_eq!(stdout(interpret(program)), "card = 0 / 0\ncard = 1 / 2\n");
}

#[test]
fn invalid_inputs_and_exhausted_search_report_a_diagnostic() {
    for expression in [
        "bracket()",
        "bracket(zero[S, 1])",
        "bracket(zero[S, 0], zero[S, 1])",
        "bracket(zero[S, 1], zero[One, 1])",
        "bracket(bar[2]{1}, bar[2]{1})",
        "bracket(zero[Stable(S, 0, 1)], zero[S, 1])",
    ] {
        let error = diagnostic(interpret(&format!(
            "space V {{}};\nlet f = \\x : V => {expression};"
        )));
        assert!(error.contains(":2:"), "{error}");
    }
    let too_high = format!("show bracket(zero[S, {}], zero[S, 1]);", usize::MAX);
    assert!(diagnostic(interpret(&too_high)).contains("overflow"));
    assert!(diagnostic(interpret("card (\\x : S => One.pt);")).contains("card"));
    assert!(diagnostic(interpret("card S.hi;")).contains("card"));
    let wide = (0..1200).map(|_| "z").collect::<Vec<_>>().join(", ");
    let error = diagnostic(interpret(&format!(
        "let z = zero[S, 1];\ncard bracket({wide});"
    )));
    assert!(error.contains("limit"), "{error}");
}
