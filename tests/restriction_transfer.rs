mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn higher_restriction_computes_nonzero_components_in_arities_two_through_four() {
    assert_eq!(
        stdout(example("ainfty-restriction")),
        "ExtAn(1, 3, 7){g0}\nExtAn(1, 4, 12){g0}\n"
    );
}

#[test]
fn transfer_over_a_one_satisfies_binary_products_and_contraction_identities() {
    let program = r"
        let a = extangen[1, 1, 1, 0];
        let b = extangen[1, 1, 2, 0];
        let unit = extangen[1, 0, 0, 0];
        assert ainfty(a, a) == yoneda(a, a);
        assert ainfty(unit, a) == a;
        assert ainfty(a, unit) == a;
        assert alift(unit, a) == zero[CobarAn(1, 0, 1)];
        assert ainfty(a, b, a) == extangen[1, 2, 4, 0];
        assert cobarD(alift(a, b)) == concat(extRep(a), extRep(b));
        let obstruction = add(concat(extRep(a), alift(b, a)), concat(alift(a, b), extRep(a)));
        assert cobarD(alift(a, b, a)) == add(obstruction, extRep(ainfty(a, b, a)));
        assert resHigher[1](a) == a;
        assert resHigher[1](a, a) == zero[ExtAn(1, 1, 2)];
        assert resHigher[1](extgen[0, 0, 0], extgen[2, 2, 0]) == zero[ExtAn(1, 1, 2)];
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn coherent_restriction_satisfies_nonzero_arity_three_and_four_morphism_equations() {
    let mut source = String::from(
        "let h0 = extgen[1, 1, 0];\n\
         let h1 = extgen[1, 2, 0];\n\
         let h0square = yoneda(h0, h0);\n",
    );
    for (inputs, generator) in [
        (&["h0", "h1", "h0"][..], "extangen[1, 2, 4, 0]"),
        (&["h1", "h0", "h1", "h0square"][..], "extangen[1, 3, 7, 0]"),
    ] {
        let (left, right) = identity(inputs);
        source.push_str(&format!(
            "assert {left} == {generator};\nassert {left} == {right};\n"
        ));
    }
    assert!(stdout(interpret(&source)).is_empty());
}

// The left side sums target operations over consecutive partitions into at
// least two blocks. The right side inserts one source operation in each slot.
fn identity(inputs: &[&str]) -> (String, String) {
    let mut left = Vec::new();
    for mask in 1..1usize << (inputs.len() - 1) {
        let mut parts = Vec::new();
        let mut start = 0;
        for end in 1..=inputs.len() {
            if end == inputs.len() || mask & (1 << (end - 1)) != 0 {
                parts.push(format!("resHigher[1]({})", inputs[start..end].join(",")));
                start = end;
            }
        }
        left.push(format!("ainfty({})", parts.join(",")));
    }
    let mut right = Vec::new();
    for length in 2..=inputs.len() {
        for start in 0..=inputs.len() - length {
            let end = start + length;
            let insertion = format!("ainfty({})", inputs[start..end].join(","));
            let mut parts = inputs[..start].to_vec();
            parts.push(&insertion);
            parts.extend_from_slice(&inputs[end..]);
            right.push(format!("resHigher[1]({})", parts.join(",")));
        }
    }
    let sum = |terms: Vec<String>| {
        terms
            .into_iter()
            .reduce(|a, b| format!("add({a},{b})"))
            .unwrap()
    };
    (sum(left), sum(right))
}

#[test]
fn higher_restriction_checks_profiles_and_bounds_partitions() {
    for (term, message) in [
        ("resHigher[1]()", "at least one Ext class"),
        ("resHigher[1](cword{[1]})", "requires an Ext class"),
        (
            "resHigher[1](extgen[1,1,0],extangen[1,1,1,0])",
            "same Steenrod algebra profile",
        ),
        (
            "resHigher[2](extangen[1,1,1,0])",
            "smaller Steenrod algebra",
        ),
        (
            "ainfty(extgen[1,1,0],extangen[1,1,1,0])",
            "same Steenrod algebra profile",
        ),
        (
            "resHigher[1](zero[Ext(0,0)],zero[Ext(0,0)])",
            "negative cohomological degree",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(message), "{error}");
    }
    let overflow = format!(
        "show resHigher[1](zero[Ext(1,{})],zero[Ext(1,1)]);",
        usize::MAX
    );
    assert!(diagnostic(interpret(&overflow)).contains("internal degree is out of range"));

    for arity in [30, 65] {
        let inputs = vec!["a"; arity].join(",");
        let program = format!("let a = extgen[1,2,0]; show resHigher[1]({inputs});");
        assert!(diagnostic(interpret(&program)).contains("consecutive partition work limit"));
        let program = format!(
            "let a = extangen[1,1,2,0]; assert resHigher[1]({inputs}) == zero[ExtAn(1,1,{})];",
            2 * arity,
        );
        assert!(stdout(interpret(&program)).is_empty());
    }
}
