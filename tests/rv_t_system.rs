mod common;

use common::{diagnostic, example, interpret, stdout};

const SETUP: &str = "
    rvideal J[[1, 2]] {
        rvmon[[1, 2]]{2, 0}, rvmon[[1, 2]]{0, 2}
    };
    let A = rvkterm[J,[1]](rvmon[[1, 2]]{1, 0});
    let B = rvkterm[J,[2]](rvmon[[1, 2]]{0, 1});
    let a = rvTorClass(A);
    let b = rvTorClass(B);
    let one = rvkterm[J,[]](rvmon[[1, 2]]{0, 0});
";

#[test]
fn four_nonzero_adjacent_classes_have_a_nonzero_fourfold_value() {
    assert_eq!(
        stdout(example("rank-massey-four")),
        concat!(
            "RVTor(Stair, 1, 2){g0}\n",
            "RVTor(Stair, 1, 2){g1}\n",
            "RVTor(Stair, 1, 2){g3}\n",
            "RVTor(Stair, 1, 2){g6}\n",
            "rvbetti[Stair, 6, 8] = 1\n",
            "RVTor(Stair, 6, 8){g0}\n"
        )
    );
}

#[test]
fn explicit_triple_agrees_with_the_complete_triple_and_twofold_is_product() {
    let flag = include_str!("../examples/rank-massey.asl");
    let source = format!(
        "{flag}
        let explicit = rvTMassey(A, B, C){{ 1 < 3 => V }};
        assert explicit == q;
        assert rvTInBracket(a, b, c, explicit) == S.hi;
        assert rvTCurvature(A, B, C){{ 1 < 3 => V }} == answer;
        "
    );
    assert!(stdout(interpret(&source)).ends_with("RVTor(Flag, 4, 6){g0}\n"));
    let twofold = format!(
        "{SETUP}
        assert rvTCurvature(A, B){{}} == rvkProd(A, B);
        assert rvTMassey(a, B){{}} == rvTorProd(a, b);
        "
    );
    assert!(stdout(interpret(&twofold)).is_empty());
}

#[test]
fn interval_errors_pinpoint_inputs_and_proper_equations() {
    for (source, expected) in [
        ("show rvTMassey(a){};", "at least two"),
        (
            "show rvTMassey(a,b,a){};",
            "equation fails on interval [0,2]",
        ),
        (
            "show rvTCurvature(rvkterm[J,[1]](rvmon[[1, 2]]{0, 0}),one){};",
            "nonzero differential on interval [0,1]",
        ),
        ("show rvTMassey(a,b,a){0<2 => A};", "RVKoszul(J, 3, 4)"),
        (
            "let z = zero[RVKoszul(J,3,4)]; show rvTMassey(a,a,a){0<2 => z, 0<2 => z};",
            "duplicate",
        ),
        ("show rvTMassey(a,a,a){0<3 => A};", "top interval"),
        ("show rvTMassey(a,a,a){2<1 => A};", "invalid interval"),
        ("show rvTMassey(a,a,a){1<2 => A};", "adjacent interval"),
        (
            "show rvTMassey(a,A,One.pt){};",
            "requires RVTor or RVKoszul",
        ),
        (
            "rvideal Other[[1, 2]] {}; show rvTMassey(a,zero[RVTor(Other,0,0)]){};",
            "same rank ideal",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{SETUP}\n{source}")));
        assert!(error.contains(expected), "{source}\n{error}");
    }
    let source = format!(
        "{SETUP} show rvTCurvature(zero[RVTor(J,{},0)],one,one){{}};",
        usize::MAX
    );
    assert!(diagnostic(interpret(&source)).contains("out of range"));
}

#[test]
fn high_arity_zero_system_does_not_build_an_ambient_basis() {
    let inputs = vec!["zero[RVTor(J,0,0)]"; 256].join(", ");
    let source = format!("{SETUP} assert rvTCurvature({inputs}){{}} == zero[RVKoszul(J,254,0)];");
    assert!(stdout(interpret(&source)).is_empty());
}
