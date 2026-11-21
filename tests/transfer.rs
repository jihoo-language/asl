mod common;

use common::{diagnostic, example, interpret, stdout};

const GENERATORS: &str = "
    let h0 = extgen[1, 1, 0];
    let h1 = extgen[1, 2, 0];
    let h2 = extgen[1, 4, 0];
    let unit = extgen[0, 0, 0];
";

#[test]
fn transfer_has_yoneda_as_its_binary_product_and_a_strict_unit() {
    let calculations = r"
        assert ainfty(h0, h2) == yoneda(h0, h2);
        assert ainfty(h2, h0) == yoneda(h2, h0);
        assert ainfty(h0, h1) == yoneda(h0, h1);
        assert ainfty(unit, unit) == unit;
        assert ainfty(unit, h0) == h0;
        assert ainfty(h0, unit) == h0;
        assert alift(unit, h0) == zero[Cobar(0, 1)];
        assert alift(h0, unit) == zero[Cobar(0, 1)];
        assert ainfty(unit, unit, h0) == zero[Ext(0, 1)];
        assert ainfty(h0, unit, h1) == zero[Ext(1, 3)];
        assert alift(h0, unit, h1) == zero[Cobar(0, 3)];
        assert alift(unit, unit, yoneda(h0, h0)) == zero[Cobar(0, 2)];
        assert ainfty(h0, zero[Ext(1, 2)], h0) == zero[Ext(2, 4)];
        assert alift(h0, zero[Ext(1, 2)], h0) == zero[Cobar(1, 4)];
        let middle = \x : Ext(1, 2) => ainfty(h0, x, h0);
        assert middle(h1) == yoneda(h1, h1);
        assert middle(zero[Ext(1, 2)]) == zero[Ext(2, 4)];

        assert alift(h0, h1) == add(cword{[0, 1]}, cword{[3]});
        assert cobarD(alift(h2, h0)) == add(concat(extRep(h2), extRep(h0)), extRep(ainfty(h2, h0)));
    ";
    assert!(stdout(interpret(&format!("{GENERATORS}{calculations}"))).is_empty());
}

#[test]
fn transferred_operations_satisfy_nontrivial_arity_four_identities() {
    let calculations = "
        let left = ainfty(h0, h1, h0);
        let right = ainfty(h1, h0, h1);
        assert left == yoneda(h1, h1);
        assert right == yoneda(h0, h2);
        assert ainfty(h0, h1) == zero[Ext(2, 3)];
        assert ainfty(h1, h0) == zero[Ext(2, 3)];
        assert ainfty(left, h1) == extgen[3, 6, 0];
        assert ainfty(left, h1) == ainfty(h0, right);

        let fourth = ainfty(h0, h1, h2, h0);
        assert fourth == yoneda(h2, h2);
        let f4 = alift(h0, h1, h2, h0);
        assert f4 == add(cword{[1, 0, 1]}, add(cword{[2, 2]}, cword{[5, 1]}));
        let u4a = concat(extRep(h0), alift(h1, h2, h0));
        let u4b = concat(alift(h0, h1), alift(h2, h0));
        let u4c = concat(alift(h0, h1, h2), extRep(h0));
        let u4 = add(u4a, add(u4b, u4c));
        let defect = add(concat(alift(h0, h1), cword{[1], [4]}), cword{[1], [0, 1], [4]});
        assert cobarD(u4) == defect;

        let m2left = alift(ainfty(h0, h1), h2, h0);
        let m2middle = alift(h0, ainfty(h1, h2), h0);
        let m2right = alift(h0, h1, ainfty(h2, h0));
        let m3left = alift(ainfty(h0, h1, h2), h0);
        let m3right = alift(h0, ainfty(h1, h2, h0));
        let corrections = add(add(m2left, m2middle), add(m2right, add(m3left, m3right)));
        assert cobarD(f4) == add(u4, add(extRep(fourth), corrections));
    ";
    assert!(stdout(interpret(&format!("{GENERATORS}{calculations}"))).is_empty());
}

#[test]
fn invalid_transfer_expressions_report_a_diagnostic() {
    for (expression, expected) in [
        ("ainfty(unit)", "at least two Ext classes"),
        ("alift()", "at least two Ext classes"),
        ("ainfty(h0, cword{[2]})", "an Ext class"),
        ("alift(unit, unit)", "negative cohomological degree"),
        ("ainfty(unit, unit, unit)", "negative cohomological degree"),
    ] {
        let error = diagnostic(interpret(&format!("{GENERATORS}show {expression};")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn transfer_names_keep_their_meaning_in_unary_calls_to_bound_functions() {
    let program = r"
        let ainfty = \x : One => x;
        let alift = \x : One => ainfty(x);
        assert alift(One.pt) == One.pt;
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn zero_transfer_handles_the_largest_representable_degree() {
    let max = usize::MAX;
    let inputs = format!("zero[Ext({max}, 0)], zero[Ext(1, 0)]");
    assert_eq!(
        stdout(interpret(&format!("show alift({inputs});"))),
        format!("Cobar({max}, 0){{0}}\n")
    );
    let error = diagnostic(interpret(&format!("show ainfty({inputs});")));
    assert!(
        error.contains("cohomological degree is out of range"),
        "{error}"
    );
    let error = diagnostic(interpret(&format!(
        "show alift(zero[Ext(1, {max})], zero[Ext(1, 1)]);"
    )));
    assert!(error.contains("internal degree is out of range"), "{error}");
}

#[test]
fn transfer_example_runs() {
    assert_eq!(
        stdout(example("ext-transfer")),
        concat!(
            "Ext(2, 4){g0}\n",
            "Ext(2, 8){g0}\n",
            "Cobar(1, 8){cword{[5, 1]} + cword{[2, 2]} + cword{[1, 0, 1]}}\n"
        )
    );
}
