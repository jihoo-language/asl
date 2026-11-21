mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn motivic_raw_cobar_keeps_global_carries_and_realizes_the_differential() {
    assert_eq!(
        stdout(example("motivic-cobar")),
        concat!(
            "MCobar(2, 4, 1){mcword[0]{([1], [1]), ([], [1])} + mcword[1]{([1], []), ([1], [])} + mcword[0]{([], [1]), ([], [0, 1])} + mcword[0]{([], [0, 1]), ([], [1])}}\n",
            "Cobar(2, 4){cword{[3], [1]} + cword{[2], [2]} + cword{[1], [0, 1]} + cword{[0, 1], [1]}}\n",
        ),
    );
}

#[test]
fn motivic_concatenation_is_balanced_and_obeys_leibniz() {
    let source = r"
        let tau1 = mcword[0]{([],[0,1])};
        let mix = mcword[0]{([],[1,1])};
        let xi2 = mcword[0]{([0,1],[])};
        let x = mcword[0]{([1],[])};
        let x2 = mcword[0]{([2],[])};
        assert cobarD(xi2) == concat(x2,x);
        assert concat(tauShift[2](tau1),tauShift[3](mix)) == tauShift[5](concat(tau1,mix));
        assert cobarD(tauShift[5](mix)) == tauShift[5](cobarD(mix));
        assert cobarD(concat(tau1,mix)) == add(concat(cobarD(tau1),mix),concat(tau1,cobarD(mix)));
        assert realize(concat(tau1,mix)) == concat(realize(tau1),realize(mix));
        assert realize(tauShift[7](mix)) == realize(mix);
        assert concat(mcword[0]{},mix) == mix;
        assert realize(mcword[3]{}) == cword{};
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn motivic_word_spaces_have_finite_weight_specific_bases() {
    let source = r"
        card(\x : MCobar(2,4,2) => S.hi);
        card(\x : MCobar(2,4,1) => S.hi);
        card(\x : MCobar(2,4,-1) => S.hi);
        card(\x : MCobar(3,3,1) => S.hi);
        card(\x : MCobar(64,128,64) => S.hi);
        card(\x : MCobar(0,0,1) => S.hi);
        card(\x : MCobar(0,0,-3) => S.hi);
        let d = \x : MCobar(1,3,1) => cobarD(x);
        assert d(mcword[0]{([],[0,1])}) == mcword[0]{([1],[]),([],[1])};
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 2 / 2\n",
            "card = 32 / 32\n",
            "card = 32 / 32\n",
            "card = 1 / 1\n",
            "card = 2 / 2\n",
            "card = 1 / 1\n",
            "card = 2 / 2\n",
        ),
    );
}

#[test]
fn motivic_coefficients_and_operators_check_types_and_weights() {
    for (source, reason) in [
        (
            "show mcword[3]{([],[])};",
            "cobar word factors must have positive degree",
        ),
        (
            "show concat(mcword[0]{([], [1])},cword{[1]});",
            "concat requires an MCobar cochain",
        ),
        (
            "show realize(cword{[1]});",
            "realize requires an MCobar cochain",
        ),
        (
            "show tauShift[1](zero[MCobar(0,0,-9223372036854775808)]);",
            "motivic weight is out of range",
        ),
        (
            "show tauShift[1](pword[3]{([], [0])});",
            "tauShift requires an MCobar cochain",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
