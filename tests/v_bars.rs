mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn graded_smith_gives_a_torsion_certificate_of_exact_order_four() {
    assert_eq!(
        stdout(example("v-bars")),
        concat!(
            "vbars[One, Flash, [0, 1]] = free [], torsion [(7, 4)]\n",
            "vbardim[One, Flash, [0, 1], 0, 1] = 1\n",
            "vbardim[One, Flash, [0, 1], 0, -1] = 0\n"
        )
    );
}

#[test]
fn flash_certificate_scales_to_ninety_seven() {
    let length = 97;
    let mut basis = Vec::new();
    let mut rules = Vec::new();
    for i in 1..=length {
        basis.push(format!("a{i}:{}, b{i}:{}", 2 * i, 2 * i - 1));
        rules.push(format!("0(a{i})=>[b{i}]"));
        if i > 1 {
            rules.push(format!("1(a{i})=>[b{}]", i - 1));
        }
    }
    let degree = 2 * length - 1;
    let program = format!(
        "
        qmodule Flash[0, 1] {{ {} }} {{ {} }};
        let cycle = vbarRep[One, Flash, [0, 1], {degree}, {length}, 0];
        let null = vbarNull[One, Flash, [0, 1], {degree}, {length}, 0];
        assert qD(null) == vRatio[{length}](cycle);
        assert vRatio[{}](qClass(cycle)) == vextgen[One, Flash, [0, 1], 0, 1, 0];
        vbars[One, Flash, [0, 1]];
        ",
        basis.join(", "),
        rules.join(", "),
        length - 1
    );
    assert_eq!(
        stdout(interpret(&program)),
        "vbars[One, Flash, [0, 1]] = free [], torsion [(193, 97)]\n"
    );
}

#[test]
fn repeated_grades_free_parts_and_contractible_parts_are_counted() {
    let source = "
        space Two { a, b };
        qmodule Double[0, 1] { a:2, b:1, c:2, d:1 } { 0(a)=>[b], 0(c)=>[d] };
        qmodule Exact[0, 1] { a:3, b:0 } { 1(a)=>[b] };
        let r0 = vbarRep[One, Double, [0, 1], 1, 1, 0];
        let r1 = vbarRep[One, Double, [0, 1], 1, 1, 1];
        assert qD(vbarNull[One, Double, [0, 1], 1, 1, 0]) == vRatio[1](r0);
        assert qD(vbarNull[One, Double, [0, 1], 1, 1, 1]) == vRatio[1](r1);
        let g0 = vextgen[One, Double, [0, 1], 0, 1, 0];
        let g1 = vextgen[One, Double, [0, 1], 0, 1, 1];
        assert add(qClass(r0), qClass(r1)) == add(g0, g1);
        vbars[One, Double, [0, 1]];
        vbars[One, Two, [0, 1]];
        vbars[One, Exact, [0, 1]];
        vbardim[One, Exact, [0, 1], 0, 0];
        vbardim[One, Two, [0, 1], -1, -7];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "vbars[One, Double, [0, 1]] = free [], torsion [(1, 1), (1, 1)]\n",
            "vbars[One, Two, [0, 1]] = free [0, 0], torsion []\n",
            "vbars[One, Exact, [0, 1]] = free [], torsion []\n",
            "vbardim[One, Exact, [0, 1], 0, 0] = 0\n",
            "vbardim[One, Two, [0, 1], -1, -7] = 2\n"
        )
    );
}

#[test]
fn finite_space_bars_recover_localized_projective_four_space() {
    let program = format!(
        "{}\n{}",
        include_str!("../examples/v-koszul.asl"),
        "
        let bar = vbarRep[One, RP4, [0, 1], 3, 2, 0];
        assert qClass(bar) == z3;
        assert qD(vbarNull[One, RP4, [0, 1], 3, 2, 0]) == vRatio[2](bar);
        assert vRatio[1](qClass(bar)) == z1;
        vbars[One, RP4, [0, 1]];
        vbardim[One, RP4, [0, 1], -1, -2];
        vbardim[One, RP4, [0, 1], 2, 9];
        "
    );
    let output = stdout(interpret(&program));
    assert!(
        output.ends_with(concat!(
            "vbars[One, RP4, [0, 1]] = free [0], torsion [(3, 2)]\n",
            "vbardim[One, RP4, [0, 1], -1, -2] = 1\n",
            "vbardim[One, RP4, [0, 1], 2, 9] = 1\n"
        )),
        "{output}"
    );
}

#[test]
fn bar_degree_queries_are_wide_and_certificates_are_validated() {
    let max = i64::MAX;
    let min = i64::MIN;
    let largest = usize::BITS - 1;
    let source = format!(
        "vbardim[One, One, [0, {largest}], {min}, {max}];\n\
         vbardim[One, One, [0, 1], {max}, -{max}];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "vbardim[One, One, [0, {largest}], {min}, {max}] = 0\n\
             vbardim[One, One, [0, 1], {max}, -{max}] = 1\n"
        )
    );
    for (source, expected) in [
        ("vbars[One, One, [0]];".into(), "exactly two"),
        ("vbars[One, One, [1, 0]];".into(), "increasing"),
        ("show vbarRep[One, One, [0, 1], 0, 1, 0];".into(), "index"),
        (
            "show vbarRep[One, One, [0, 1], 0, 0, 0];".into(),
            "positive",
        ),
        (
            format!("show vbarNull[One, One, [0, 1], {min}, 1, 0];"),
            "out of range",
        ),
        (
            "show vRatio[1](zero[VExt(One, One, [0], 0, 0)]);".into(),
            "exactly two",
        ),
    ] {
        let error = diagnostic(interpret(&source));
        assert!(error.contains(expected), "{source}\n{error}");
    }
}
