mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn odd_coherent_restriction_detects_nonzero_binary_and_ternary_components() {
    assert_eq!(
        stdout(example("odd-coherent-restriction")),
        "pextdim[3, 3, 15] = 0\nPExtAn(3, 1, 3, 15){g0}\nPExtAn(3, 1, 3, 15){g0}\n"
    );
}

#[test]
fn odd_five_coherent_restriction_uses_the_exact_incoming_correction() {
    let source = r"
        let a = pextgen[5, 5, 5, 0];
        let h = pextgen[5, 1, 40, 0];
        let lift = alift(a, h);
        assert cobarD(lift) == pscale[4](concat(extRep(a), extRep(h)));
        let cycle = resA[1](lift);
        assert cobarD(cycle) == zero[PCobarAn(5, 1, 6, 45)];
        assert extClass(cycle) == pextangen[5, 1, 5, 45, 0];
        assert resHigher[1](a, h) == extClass(cycle);
        assert resHigher[0](a, h) == zero[PExtAn(5, 0, 5, 45)];
    ";
    assert!(stdout(interpret(source)).is_empty());
}

#[test]
fn odd_finite_transfer_uses_the_same_suspension_and_coefficient_signs() {
    let source = r"
        let x = pextangen[3, 1, 1, 4, 0];
        let a = pextangen[3, 1, 1, 1, 0];
        let one = pscale[2](pextangen[3, 1, 0, 0, 0]);
        assert ainfty(x, a) == yoneda(x, a);
        assert ainfty(one, x) == pscale[2](x);
        assert ainfty(x, one) == pscale[2](x);
        assert alift(x, one) == zero[PCobarAn(3, 1, 0, 4)];
        assert alift(x, x) == resA[1](pword[3]{([2], [])});
        assert cobarD(alift(x, x)) == pscale[2](concat(extRep(x), extRep(x)));
        assert ainfty(x, x, x) == pscale[2](pextangen[3, 1, 2, 12, 0]);
        assert resHigher[1](one) == one;
        assert resHigher[1](x, x, x) == zero[PExtAn(3, 1, 1, 12)];
        let b = pextgen[3, 2, 12, 0];
        let a2 = pextgen[3, 2, 2, 0];
        let alpha = pextgen[3, 1, 1, 0];
        let period = resHigher[1](a2, b, alpha);
        assert resHigher[1](pscale[2](a2), b, alpha) == pscale[2](period);
        assert resHigher[1](a2, pscale[2](b), pscale[2](alpha)) == period;
        assert add(period, period) == resHigher[1](a2, add(b, b), alpha);
        assert resHigher[1](pextgen[3, 0, 0, 0], b) == zero[PExtAn(3, 1, 1, 12)];
    ";
    assert!(stdout(interpret(source)).is_empty());
}

#[test]
fn suspended_morphism_identities_have_the_expected_signed_third_and_fourth_components() {
    let mut source = String::from(
        "let a = pextgen[3, 1, 1, 0];\n\
         let x = pextgen[3, 1, 4, 0];\n\
         let square = pextgen[3, 2, 2, 0];\n\
         let b = pextgen[3, 2, 12, 0];\n",
    );
    let a = Expression::new("a", 1, 1);
    let x = Expression::new("x", 1, 4);
    let square = Expression::new("square", 2, 2);
    let b = Expression::new("b", 2, 12);
    for (inputs, expected) in [
        (
            vec![a.clone(), x.clone(), x],
            "pscale[2](pextangen[3, 1, 2, 9, 0])",
        ),
        (vec![a.clone(), square, b, a], "pextangen[3, 1, 4, 16, 0]"),
    ] {
        let (left, right) = identity(&inputs);
        writeln!(
            source,
            "assert {left} == {expected};\nassert {left} == {right};"
        )
        .unwrap();
    }
    assert!(stdout(interpret(&source)).is_empty());
}

#[derive(Clone)]
struct Expression {
    text: String,
    s: usize,
    t: usize,
}

impl Expression {
    fn new(text: impl Into<String>, s: usize, t: usize) -> Self {
        Self {
            text: text.into(),
            s,
            t,
        }
    }

    fn odd_stem(&self) -> bool {
        (self.t ^ self.s) & 1 != 0
    }
}

fn signed(expression: String, negative: bool) -> String {
    if negative {
        format!("pscale[2]({expression})")
    } else {
        expression
    }
}

fn suspended(operation: &str, inputs: &[Expression]) -> Expression {
    let arguments = inputs
        .iter()
        .map(|term| term.text.as_str())
        .collect::<Vec<_>>()
        .join(", ");
    let negative = inputs
        .iter()
        .rev()
        .skip(1)
        .step_by(2)
        .fold(false, |sign, term| sign ^ term.odd_stem());
    let shift = if operation == "ainfty" { 2 } else { 1 };
    let s = inputs.iter().map(|term| term.s).sum::<usize>() + shift - inputs.len();
    let t = inputs.iter().map(|term| term.t).sum();
    Expression::new(signed(format!("{operation}({arguments})"), negative), s, t)
}

fn identity(inputs: &[Expression]) -> (String, String) {
    let arity = inputs.len();
    let mut left = Vec::new();
    // b_D R on the whole tensor, with at least two consecutive output blocks.
    for mask in 1..1usize << (arity - 1) {
        let mut parts = Vec::new();
        let mut start = 0;
        for end in 1..=arity {
            if end == arity || mask & (1 << (end - 1)) != 0 {
                parts.push(suspended("resHigher[1]", &inputs[start..end]));
                start = end;
            }
        }
        left.push(suspended("ainfty", &parts).text);
    }
    let mut right = Vec::new();
    // R b_C includes the suspended-degree sign of the unchanged prefix.
    for length in 2..=arity {
        for start in 0..=arity - length {
            let end = start + length;
            let insertion = suspended("ainfty", &inputs[start..end]);
            let mut parts = inputs[..start].to_vec();
            parts.push(insertion);
            parts.extend_from_slice(&inputs[end..]);
            let negative = inputs[..start]
                .iter()
                .filter(|term| !term.odd_stem())
                .count()
                % 2
                == 1;
            right.push(signed(suspended("resHigher[1]", &parts).text, negative));
        }
    }
    let sum = |terms: Vec<String>| {
        terms
            .into_iter()
            .reduce(|a, b| format!("add({a}, {b})"))
            .unwrap()
    };
    (sum(left), sum(right))
}

#[test]
fn odd_higher_restriction_rejects_mixed_fields_profiles_and_oversized_partitions() {
    for (term, reason) in [
        (
            "resHigher[1](pextgen[3,1,1,0], pextgen[5,1,1,0])",
            "same prime",
        ),
        (
            "resHigher[1](pextgen[3,1,1,0], extgen[1,1,0])",
            "a PExt class",
        ),
        (
            "resHigher[1](pextangen[3,2,1,1,0], pextgen[3,1,1,0])",
            "same Steenrod algebra profile",
        ),
        (
            "resHigher[2](pextangen[3,1,1,1,0])",
            "smaller Steenrod algebra",
        ),
        (
            "alift(pextangen[3,1,1,1,0], pextgen[3,1,1,0])",
            "same Steenrod algebra profile",
        ),
        (
            "resHigher[1](zero[PExt(3,0,0)], zero[PExt(3,0,0)])",
            "negative cohomological degree",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
    let max = usize::MAX;
    let error = diagnostic(interpret(&format!(
        "show resHigher[1](zero[PExt(3,1,{max})], pextgen[3,1,1,0]);"
    )));
    assert!(error.contains("internal degree is out of range"), "{error}");
    let args = vec!["x"; 65].join(",");
    let error = diagnostic(interpret(&format!(
        "let x = pextgen[3,1,1,0]; show resHigher[1]({args});"
    )));
    assert!(
        error.contains("consecutive partition work limit"),
        "{error}"
    );
}
