mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn real_projective_four_space_has_a_nonzero_second_bockstein() {
    assert_eq!(
        stdout(example("vbock-rp4")),
        concat!(
            "VBock(One, RP4, [0, 1], 2, 1, 2, 2){g0}\n",
            "VBock(One, RP4, [0, 1], 2, 0, -1, 2){g0}\n",
            "vbockdim[One, RP4, [0, 1], 0, 0, 2, 1] = 1\n",
            "vbockdim[One, RP4, [0, 1], 1, 0, 2, 1] = 0\n",
            "vbockdim[One, RP4, [0, 1], 2, 0, 2, 0] = 1\n",
            "vbockdim[One, RP4, [0, 1], 3, 0, 2, 0] = 0\n",
            "vbockdim[One, RP4, [0, 1], 3, 1, 2, 2] = 0\n",
            "vextdim[One, RP4, [0, 1], 0, 0] = 1\n",
            "vextdim[One, RP4, [0, 1], 0, 2] = 0\n",
            "vextdim[One, RP4, [0, 1], 1, 2] = 0\n",
        )
    );
    let source = include_str!("../examples/vbock-rp4.asl");
    for page in ["e2", "shifted"] {
        let error = diagnostic(interpret(&format!("{source}\nshow vbNext({page});")));
        assert!(error.contains("vbNext requires a cycle for d_2"), "{error}");
    }
}

#[test]
fn incoming_boundaries_are_retained_at_the_smallest_signed_degree() {
    let min = i64::MIN;
    let t = min + 1;
    let program = format!(
        "{}
         let oneVariable = vLocal(qRes[0](mword[One, 0, 0, RP4, 1, 0](cword{{}})));
         let bottom = vbClass[0, 0](vShift[{min}](oneVariable));
         assert bottom == vbockgen[One, RP4, [0], 0, {min}, {t}, 0, 0];
         assert vbNext(bottom) == zero[VBock(One, RP4, [0], 1, {min}, {t}, 0)];",
        include_str!("../examples/vbock-rp4.asl")
    );
    assert!(stdout(interpret(&program)).ends_with("vextdim[One, RP4, [0, 1], 1, 2] = 0\n"));
}

#[test]
fn composition_adds_all_lower_exponents_and_v_shifts_preserve_filtration() {
    let program = r"
        let a = vLocal(qRes[0, 1, 2](cword{[1]}));
        let b = vLocal(qRes[0, 1, 2](cword{[0, 1]}));
        let left = vbClass[3, 1](a);
        let right = vbClass[3, 1](b);
        let product = ecomp(left, right);
        let expected = vbockgen[One, One, [0, 1, 2], 3, 2, 4, 2, 0];
        assert product == expected;
        assert vbRep(product) == ecomp(a, b);
        let translated = vShift[-6](product);
        assert translated == vbockgen[One, One, [0, 1, 2], 3, -4, -38, 2, 0];
        assert vShift[6](translated) == product;
        assert vbNext(vbPrev(translated)) == translated;
        assert vbD(translated) == zero[VBock(One, One, [0, 1, 2], 3, -3, -38, 5)];
        let shift = \x : VBock(One, One, [0, 1, 2], 3, 2, 4, 2) => vShift[-6](x);
        assert shift(product) == translated;
        assert add(product, product) == zero[VBock(One, One, [0, 1, 2], 3, 2, 4, 2)];
        betti VBock(One, One, [0, 1, 2], 3, -4, -38, 2);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(VBock(One, One, [0, 1, 2], 3, -4, -38, 2)) = [2]\n"
    );
}

#[test]
fn bockstein_operators_check_pages_indices_filtrations_and_degree_bounds() {
    for (expression, reason) in [
        ("zero[VBock(One, One, [], 0, 0, 0, 0)]", "nonempty"),
        ("vbClass[0, 0](qRes[0](cword{}))", "VKoszul"),
        (
            "vbClass[0, 1](vLocal(qRes[0](cword{})))",
            "below the requested",
        ),
        (
            "vbPrev(zero[VBock(One, One, [0], 0, 0, 0, 0)])",
            "greater than zero",
        ),
        ("vbNext(cword{})", "a VBock class"),
        (
            "ecomp(zero[VBock(One, One, [0], 0, 0, 0, 0)], zero[VBock(One, One, [0], 1, 0, 0, 0)])",
            "same page",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }
    let max = usize::MAX;
    let top = i64::MAX;
    for (expression, reason) in [
        (
            format!("vbNext(zero[VBock(One, One, [0], {max}, 0, 0, 0)])"),
            "page is out of range",
        ),
        (
            format!("vbD(zero[VBock(One, One, [0], 0, {top}, 0, 0)])"),
            "degree is out of range",
        ),
        (
            format!("vbD(zero[VBock(One, One, [0], 1, 0, 0, {max})])"),
            "filtration is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }
}
