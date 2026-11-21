mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn the_eight_fold_obstruction_and_two_choices_for_a_triple_are_visible() {
    let expected = "H2(8){0}\n".repeat(6) + "H2(8){carry}\nH2(2){0}\nH2(2){carry}\n";
    assert_eq!(stdout(example("massey-tower")), expected);
}

#[test]
fn cyclic_cochains_and_classes_are_discrete_functional_values() {
    let source = r"
        let toggle = \p : Bar(4) * Bar(4) => add(fst(p), snd(p));
        let a = bar[4]{3, 2, 3, 1};
        assert toggle((a, bar[4]{2})) == bar[4]{1};
        assert binom[8, 3] == bar[8]{3, 7};
        assert binom[8, 8] == zero[Bar(8)];
        assert add(carry[4], carry[4]) == zero[H2(4)];
        assert carry[5] == zero[H2(5)];
        let x = binom[4, 1];
        assert massey(x, x){} == zero[H2(4)];
        show a;
        show binom[64, 64];
        betti Bar(4);
        betti H2(4);
        betti H2(5);
        homotopy (\x : One => carry[4]) ~ (\x : One => zero[H2(4)]);
        space H2 { p };
        space Bar { p };
        let massey = \v : H2 => v;
        let powerdef = \v : Bar => v;
        assert massey(H2.p) == H2.p;
        assert powerdef(Bar.p) == Bar.p;
    ";
    assert_eq!(
        stdout(interpret(source)),
        "bar[4]{1, 2}\nbar[64]{}\nbetti(Bar(4)) = [8]\nbetti(H2(4)) = [2]\nbetti(H2(5)) = [1]\nhomotopy = no\n"
    );
    assert!(diagnostic(interpret("betti Bar(64);")).contains("limit"));
}

#[test]
fn the_next_witness_fails_at_a_specific_pair_of_residues() {
    let source = (1..=8)
        .map(|k| format!("binom[8, {k}]"))
        .collect::<Vec<_>>()
        .join(", ");
    let failure = diagnostic(interpret(&format!("show powerdef({source});")));
    assert!(
        failure.contains("[0,8]") && failure.contains("(1,7)"),
        "{failure}"
    );

    for (source, interval, pair) in [
        ("massey(zero[Bar(4)], bar[4]{2, 3}){}", "[1,2]", "(1,1)"),
        ("powerdef(binom[3, 1])", "[0,1]", "(1,2)"),
        (
            "massey(binom[4, 1], binom[4, 1], binom[4, 1]){0 < 2 => binom[4, 2], 1 < 3 => bar[4]{}}",
            "[1,3]",
            "(1,1)",
        ),
    ] {
        let failure = diagnostic(interpret(&format!("show {source};")));
        assert!(failure.contains(interval) && failure.contains(pair), "{failure}");
    }
}

#[test]
fn malformed_defining_systems_are_rejected_even_under_an_empty_domain() {
    let prefix = "space V {}; let z = zero[Bar(4)];\nlet f = \\v : V => ";
    for expr in [
        "massey(z, z, z){0 < 2 => z}",
        "massey(z, z, z){0 < 2 => z, 0 < 2 => z, 1 < 3 => z}",
        "massey(z, z){0 < 2 => z}",
        "massey(z, z, z){0 < 1 => z, 1 < 3 => z}",
        "massey(z, z, z){2 < 1 => z, 1 < 3 => z}",
        "massey(z, z, z){0 < 9 => z, 1 < 3 => z}",
        "massey(z, z, z){0 < 2 => z, 1 < 3 => carry[4]}",
        "massey(z, zero[Bar(8)]){}",
        "massey(z){}",
        "powerdef()",
        "powerdef(z, bar[8]{})",
        "powerdef(carry[4])",
        "add(z, zero[H2(4)])",
        "add(carry[4], carry[8])",
        "bar[4]{0}",
        "bar[4]{4}",
        "binom[4, 0]",
        "zero[Bar(1)]",
        "carry[0]",
    ] {
        let failure = diagnostic(interpret(&format!("{prefix}{expr};")));
        assert!(failure.contains(":2:"), "{expr}\n{failure}");
    }
}

#[test]
fn the_sixty_four_fold_product_is_computed_from_literal_cochains() {
    // The bits here come from Pascal's recursion, independently of binom[m,k].
    let mut row = [false; 64];
    row[0] = true;
    let mut supports = vec![Vec::new(); 64];
    for g in 1..64 {
        for k in (1..=g).rev() {
            row[k] ^= row[k - 1];
        }
        for k in 1..=g {
            if row[k] {
                supports[k].push(g.to_string());
            }
        }
    }
    let mut source = String::new();
    let mut names = Vec::new();
    for (k, support) in supports.iter().enumerate().skip(1) {
        writeln!(source, "let b{k} = bar[64]{{{}}};", support.join(", ")).unwrap();
        names.push(format!("b{k}"));
    }
    for n in [2, 4, 8, 16, 32, 63, 64] {
        writeln!(source, "show powerdef({});", names[..n - 1].join(", ")).unwrap();
    }
    assert_eq!(
        stdout(interpret(&source)),
        "H2(64){0}\n".repeat(6) + "H2(64){carry}\n"
    );
}

#[test]
fn repeated_witnesses_match_the_full_triangle_across_word_boundaries() {
    // Use (1 + t + t² + t³)^g, whose generator has three active coefficients.
    // For both orders, (1 + t + t² + t³)^m starts with 1 + t⁴.
    for m in [68, 132] {
        let cochains = format!(
            "let x = binom[{m}, 1];\n\
             let b = add(x, binom[{m}, 2]);\n\
             let c = add(x, binom[{m}, 3]);\n"
        );
        let full = "massey(x, x, x, x){0 < 2 => b, 1 < 3 => b, 2 < 4 => b, 0 < 3 => c, 1 < 4 => c}";
        let source = format!("{cochains}assert {full} == powerdef(x, b, c);\nshow {full};");
        assert_eq!(stdout(interpret(&source)), format!("H2({m}){{carry}}\n"));

        let last = m - 1;
        let corrupt = format!("{cochains}show powerdef(x, b, add(c, bar[{m}]{{{last}}}));");
        let failure = diagnostic(interpret(&corrupt));
        assert!(
            failure.contains("[0,3]") && failure.contains(&format!("(1,{})", m - 2)),
            "{failure}"
        );
    }
}

#[test]
fn a_thousand_input_product_reaches_its_obstruction_within_the_work_budget() {
    let coefficients = (1..1024)
        .map(|k| format!("binom[1024, {k}]"))
        .collect::<Vec<_>>();
    let source = format!("show powerdef({});", coefficients.join(", "));
    assert_eq!(stdout(interpret(&source)), "H2(1024){carry}\n");
}

type Matrix = [u8; 5];
const IDENTITY: Matrix = [1, 2, 4, 8, 16];

fn product(left: Matrix, right: Matrix) -> Matrix {
    left.map(|row| {
        let mut result = 0;
        for (k, image) in right.iter().enumerate() {
            if row & (1 << k) != 0 {
                result ^= image;
            }
        }
        result
    })
}

#[test]
fn all_four_fold_defining_systems_on_small_cyclic_groups_agree_with_matrix_lifting() {
    // Each lift of the group generator has nine free entries after the
    // top-right corner is fixed at zero. Its powers give the actual cochains.
    for m in [2, 3, 4, 6, 8] {
        let mut source = String::new();
        for support in (0..1u16 << m).step_by(2) {
            let residues: Vec<_> = (1..m)
                .filter(|g| support & (1 << g) != 0)
                .map(|g| g.to_string())
                .collect();
            writeln!(
                source,
                "let c{support} = bar[{m}]{{{}}};",
                residues.join(", ")
            )
            .unwrap();
        }
        for mask in 0..512 {
            let mut generator = IDENTITY;
            let mut bit = 1;
            for (i, row) in generator.iter_mut().enumerate() {
                for j in i + 1..5 {
                    if (i, j) != (0, 4) {
                        if mask & bit != 0 {
                            *row |= 1 << j;
                        }
                        bit <<= 1;
                    }
                }
            }
            let mut power = IDENTITY;
            let mut powers = Vec::new();
            for _ in 0..m {
                powers.push(power);
                power = product(power, generator);
            }
            let obstructed = power[0] & 16 != 0 && m % 2 == 0;
            power[0] &= !16;
            if power != IDENTITY {
                continue;
            }
            let cochain = |i: usize, j: usize| {
                let mut support = 0;
                for (g, matrix) in powers.iter().enumerate() {
                    if matrix[i] & (1 << j) != 0 {
                        support |= 1 << g;
                    }
                }
                format!("c{support}")
            };
            let inputs = (0..4).map(|i| cochain(i, i + 1)).collect::<Vec<_>>();
            let entries = [(0, 2), (1, 3), (2, 4), (0, 3), (1, 4)]
                .map(|(i, j)| format!("{i} < {j} => {}", cochain(i, j)));
            let expected = if obstructed {
                format!("carry[{m}]")
            } else {
                format!("zero[H2({m})]")
            };
            writeln!(
                source,
                "assert massey({}){{{}}} == {expected};",
                inputs.join(", "),
                entries.join(", ")
            )
            .unwrap();
            if (0..4).all(|i| (i + 1..4).all(|j| cochain(i, j) == cochain(i + 1, j + 1))) {
                let diagonals = [cochain(0, 1), cochain(0, 2), cochain(0, 3)];
                writeln!(
                    source,
                    "assert powerdef({}) == {expected};",
                    diagonals.join(", ")
                )
                .unwrap();
            }
        }
        assert!(stdout(interpret(&source)).is_empty());
    }
}
