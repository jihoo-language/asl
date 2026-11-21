mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn lightning_flash_has_exact_v0_torsion_order_four() {
    assert_eq!(
        stdout(example("q-module")),
        "QExt(One, Flash, [0, 1], 3, 10){g0}\nqextdim[One, Flash, [0, 1], 4, 11] = 0\n"
    );
}

#[test]
fn declared_endpoint_actions_continue_under_localization() {
    let program = format!(
        "{}\n{}",
        include_str!("../examples/q-module.asl"),
        r"
        let localized = vLocal(third);
        let shifted = vShift[-5](localized);
        assert shifted == vextgen[One, Flash, [0,1], -2, -5, 0];
        assert vShift[5](shifted) == localized;
        assert qD(vLocal(witness)) == vLocal(fourth);
        assert ecomp(vLocal(qid[Flash,[0,1]]), localized) == localized;
    "
    );
    assert_eq!(
        stdout(interpret(&program)),
        "QExt(One, Flash, [0, 1], 3, 10){g0}\nqextdim[One, Flash, [0, 1], 4, 11] = 0\n"
    );
}

#[test]
fn long_sparse_witnesses_use_polynomial_exponents() {
    let length = 32;
    let shift = 1_000_000;
    let mut generators = Vec::new();
    let mut rules = Vec::new();
    for j in 1..=length {
        generators.push(format!("a{j}:{}", 2 * j));
        generators.push(format!("b{j}:{}", 2 * j - 1));
        rules.push(format!("0(a{j})=>[b{j}]"));
        if j > 1 {
            rules.push(format!("1(a{j})=>[b{}]", j - 1));
        }
    }
    let mut program = format!(
        "qmodule M[0,1] {{ {} }} {{ {} }};\n",
        generators.join(","),
        rules.join(",")
    );
    for j in 1..=length {
        writeln!(
            program,
            "let a{j} = qterm[One,0,0,M,{},0,[0,1]]{{{},{}}};",
            2 * j,
            j - 1,
            shift + length - j
        )
        .unwrap();
        if j == 1 {
            program.push_str("let w1 = a1;\n");
        } else {
            writeln!(program, "let w{j} = add(w{}, a{j});", j - 1).unwrap();
        }
    }
    let last = 2 * length - 1;
    writeln!(
        program,
        "assert qD(w{length}) == qterm[One,0,0,M,{last},0,[0,1]]{{{length},{shift}}};"
    )
    .unwrap();
    writeln!(program, "assert qClass(qterm[One,0,0,M,{last},0,[0,1]]{{{},{shift}}}) == qextgen[One,M,[0,1],{},{},0];", length - 1, shift + length - 1, 3*shift + 3*length - 2).unwrap();
    assert_eq!(stdout(interpret(&program)), "");

    let max = usize::MAX;
    let sparse = format!(
        "qmodule Sparse[] {{ h:{max} }} {{}};\n\
         assert qClass(qterm[Sparse,{max},0,Sparse,{max},0,[]]{{}}) == qid[Sparse,[]];"
    );
    assert_eq!(stdout(interpret(&sparse)), "");
}

#[test]
fn modules_can_have_degree_zero_actions_and_commuting_nonscalar_sums() {
    let source = r"
        qmodule Edge[0] { u:1, v:0 } { 0(u)=>[v] };
        assert qD(qterm[One,0,0,Edge,1,0,[0]]{0}) == qterm[One,0,0,Edge,0,0,[0]]{1};
        assert qD(qterm[Edge,0,0,One,0,0,[0]]{0}) == qterm[Edge,1,0,One,0,0,[0]]{1};
        assert qClass(qterm[One,0,0,Edge,0,0,[0]]{1}) == zero[QExt(One,Edge,[0],1,1)];
        assert qClass(qRep(qid[Edge,[0]])) == qid[Edge,[0]];

        qmodule Square[0,1] { x:4, y:3, z:1, w:0, w2:0 } {
            0(x)=>[y], 1(x)=>[z], 1(y)=>[w2,w], 0(z)=>[w,w2]
        };
        assert qD(qD(qterm[One,0,0,Square,4,0,[0,1]]{0,0})) == zero[QKoszul(One,Square,[0,1],2,4)];
        assert qD(qterm[One,0,0,Square,1,0,[0]]{0}) == add(qterm[One,0,0,Square,0,0,[0]]{1}, qterm[One,0,0,Square,0,1,[0]]{1});

        qmodule Twice[0] { u:1, v:0 } { 0(u)=>[v,v] };
        assert qD(qterm[One,0,0,Twice,1,0,[0]]{0}) == zero[QKoszul(One,Twice,[0],1,1)];
        qmodule Empty[0] {} {};
        assert qid[Empty,[0]] == zero[QExt(Empty,Empty,[0],0,0)];
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn declarations_validate_rules_and_keep_modules_out_of_space_types() {
    for (program, expected) in [
        (
            "qmodule M[0] { a:2,b:1,c:0 } { 0(a)=>[b],0(b)=>[c] };",
            "squared is nonzero",
        ),
        (
            "qmodule M[0,1] { a:4,b:3,c:0 } { 0(a)=>[b],1(b)=>[c] };",
            "do not commute",
        ),
        (
            "qmodule M[0] { a:2,b:0 } { 0(a)=>[b] };",
            "drop degree by 1",
        ),
        (
            "qmodule M[0] { a:1 } { 0(a)=>[b] };",
            "unknown module generator",
        ),
        (
            "qmodule M[0] { a:1 } { 1(a)=>[] };",
            "not declared for this module",
        ),
        (
            "qmodule M[0] { a:1 } { 0(a)=>[],0(a)=>[] };",
            "duplicate Q_0 rule",
        ),
        ("qmodule M[0] { a:1,a:2 } {};", "duplicate module generator"),
        (
            "qmodule M[0] { a:1 } {}; show qid[M,[1]];",
            "not declared for this module",
        ),
        ("qmodule M[0] { a:1 } {}; betti M;", "unknown space"),
        (
            "qmodule M[0] { a:1 } {}; show qterm[One,0,0,M,1,1,[0]]{0};",
            "generator index",
        ),
        (
            "show qterm[One,0,0,One,0,0,[0,1]]{0};",
            "one exponent per selected primitive",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
