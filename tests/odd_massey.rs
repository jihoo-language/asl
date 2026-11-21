mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn divided_powers_give_the_negative_nonzero_bockstein_classes() {
    assert_eq!(
        stdout(example("odd-massey-power")),
        "PExt(3, 2, 12){2*g0}\nPExt(5, 2, 40){4*g0}\n"
    );
}

#[test]
fn mixed_stems_use_the_left_bar_sign_and_give_closed_curvature() {
    let program = r"
        let x = pword[3]{([1], [])};
        let a = pword[3]{([], [0])};
        let T = pword[3]{([], [1])};
        let mixed = pword[3]{([1], [0])};
        let L = add(mixed, pscale[2](T));
        assert cobarD(T) == concat(x, a);
        assert cobarD(L) == pscale[2](concat(a, x));

        let result = pCurvature(x, a, x){0<2 => T, 1<3 => L};
        let expected = add(concat(x, mixed), pscale[2](add(concat(x, T), concat(T, x))));
        assert result == expected;
        assert cobarD(result) == zero[PCobar(3, 3, 9)];
        let class = pMassey(x, extClass(a), x){0<2 => T, 1<3 => L};
        assert class == extClass(result);
        show class;

        let reversed = pCurvature(a, x, a){0<2 => L, 1<3 => T};
        assert reversed == cobarD(pword[3]{([], [0, 1])});
        assert cobarD(reversed) == zero[PCobar(3, 3, 6)];
        assert pMassey(a, x, a){0<2 => L, 1<3 => T} == zero[PExt(3, 2, 6)];

        assert pCurvature(a, a){} == pscale[2](concat(a, a));
        assert pPowerDef(extClass(a)) == pscale[2](yoneda(extClass(a), extClass(a)));
        assert pMassey(zero[PCobar(3, 1, 4)], a, zero[PCobar(3, 1, 4)]){}
            == zero[PExt(3, 2, 9)];

        let ax = concat(a, x);
        let W = concat(a, T);
        let K = pword[3]{([], [0, 1])};
        let changed = add(W, cobarD(K));
        let original = pCurvature(ax, a, x){0<2 => W, 1<3 => L};
        let varied = pCurvature(ax, a, x){0<2 => changed, 1<3 => L};
        assert varied == add(original, pscale[2](cobarD(concat(K, x))));
        assert extClass(original) == extClass(varied);
    ";
    let output = stdout(interpret(program));
    assert!(output.starts_with("PExt(3, 2, 9){"), "{output}");
    assert!(!output.ends_with("{0}\n"), "{output}");
}

#[test]
fn malformed_defining_systems_report_the_failing_input_or_interval() {
    let setup = "let x = pword[3]{([1], [])};\nlet b = pword[3]{([2], [])};\n";
    for (expression, reason) in [
        ("pMassey(x, x, x){}", "equation fails on interval [0,2]"),
        (
            "pCurvature(b, x){}",
            "adjacent input on interval [0,1] is not closed",
        ),
        ("pMassey(x, pword[5]{([1], [])}){}", "common odd prime"),
        ("pMassey(x, x, x){0<2 => x}", "witness for interval [0,2]"),
        (
            "pMassey(x, x, x){0<2 => pscale[2](b), 0<2 => pscale[2](b)}",
            "duplicate witness",
        ),
        ("pMassey(x, x, x){0<1 => b}", "adjacent interval"),
        ("pMassey(x, x, x){0<3 => b}", "top interval"),
        ("pMassey(x, x, x){2<4 => b}", "invalid interval"),
        ("pPowerDef(x, b)", "equation fails on interval [0,2]"),
        ("pPowerDef()", "requires at least one"),
        ("pCurvature(x){}", "requires at least two"),
        (
            "pMassey(pword[3]{}, pword[3]{}, pword[3]{}){}",
            "negative cobar length",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }
    let max = usize::MAX;
    let error = diagnostic(interpret(&format!(
        "show pPowerDef(zero[PCobar(3, {max}, 0)]);"
    )));
    assert!(error.contains("cobar length is out of range"), "{error}");
    let error = diagnostic(interpret(&format!(
        "show pCurvature(zero[PCobar(3, 1, {max})], pword[3]{{([1], [])}}){{}};"
    )));
    assert!(error.contains("internal degree is out of range"), "{error}");
}

#[test]
fn raw_curvature_uses_sparse_terms_without_enumerating_cohomology() {
    let program = r"
        let x = pword[3]{([59049], [])};
        let z = zero[PCobar(3, 1, 236196)];
        assert pCurvature(x, x){} == concat(x, x);
        assert pCurvature(z, x){} == zero[PCobar(3, 2, 472392)];
        assert pMassey(z, x){} == zero[PExt(3, 2, 472392)];
        let u = pword[4294967291]{};
        assert pCurvature(u, u){} == pscale[4294967290](u);
    ";
    assert!(stdout(interpret(program)).is_empty());
}
