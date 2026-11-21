mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn explicit_p_five_and_chosen_integral_higher_values_detect_greek_beta() {
    assert_eq!(stdout(example("bp-higher")), "BPExt(5, 2, 2, 40){4*g0}\n");
    let mut program = String::new();
    let mut expected = String::new();
    for p in [3, 5, 7, 11] {
        let w = 2 * (p - 1);
        let x = format!("x{p}");
        writeln!(program, "let {x} = bpextgen[{p},2,1,{w},0];").unwrap();
        writeln!(
            program,
            "assert bpCup({x},{x}) == zero[BPExt({p},2,2,{})];",
            2 * w
        )
        .unwrap();
        for k in 3..p {
            let args = vec![x.as_str(); k].join(",");
            writeln!(
                program,
                "assert bpChosenHigher({args}) == zero[BPExt({p},2,2,{})];",
                k * w
            )
            .unwrap();
        }
        let args = vec![x.as_str(); p].join(",");
        let t = p * w;
        writeln!(
            program,
            "let b{p} = bpExtClass[{t}](bpDelta(bpDelta(bpReduce[2](bpv[{p},2,0,2]))));"
        )
        .unwrap();
        writeln!(program, "assert bpChosenHigher({args}) == bpNeg(b{p});").unwrap();
        writeln!(
            program,
            "assert bpScale[{p}](bpChosenHigher({args})) == zero[BPExt({p},2,2,{t})];"
        )
        .unwrap();
        writeln!(program, "show bpChosenHigher({args});").unwrap();
        writeln!(expected, "BPExt({p}, 2, 2, {t}){{{}*g0}}", p - 1).unwrap();
    }
    assert_eq!(stdout(interpret(&program)), expected);
}

#[test]
fn even_degree_twist_controls_the_sign_at_three_and_four_inputs() {
    let program = r"
        let beta = bpExtClass[12](bpDelta(bpDelta(bpReduce[2](bpv[3,2,0,2]))));
        let alpha = bpextgen[3,2,1,4,0];
        let expected = bpNeg(bpCup(beta,alpha));
        let tr = bpHigherSystem[3](zero[BPExt(3,2,2,8)],zero[BPExt(3,2,1,4)],alpha,bpExtRep(beta),zero[BP(3,2,1)]);
        assert tr == expected;
        let four = bpHigherSystem[4](zero[BPExt(3,2,2,4)],zero[BPExt(3,2,1,4)],zero[BPExt(3,2,1,4)],alpha,
                                     zero[BP(3,2,2)],zero[BP(3,2,1)],zero[BP(3,2,1)],bpExtRep(beta),zero[BP(3,2,1)]);
        assert four == expected;
        show four;
    ";
    assert_eq!(stdout(interpret(program)), "BPExt(3, 2, 3, 16){g0}\n");
}

#[test]
fn explicit_altered_primitives_select_both_elements_of_the_triple_coset() {
    let program = r"
        let x = bpextgen[2,2,1,2,0];
        let y = bpextgen[2,2,1,4,0];
        let t = bpt[2,2,1,1,1];
        let t2 = bpt[2,2,1,1,2];
        let vvt = bpMul(bpPow[2](bpv[2,2,1,1]),t);
        let u02 = add(add(bpNeg(vvt),bpMul(bpint[2,2,1,-2],bpPow[3](t))),bpMul(bpint[2,2,1,5],t2));
        let u13 = add(add(vvt,bpPow[3](t)),bpMul(bpint[2,2,1,-3],t2));
        let shifted = add(u02,bpExtRep(bpextgen[2,2,1,6,0]));
        let a = bpHigherSystem[3](x,y,x,u02,u13);
        let b = bpHigherSystem[3](x,y,x,shifted,u13);
        let bracket = bpTriple(x,y,x);
        bpTripleContains(bracket,a);
        bpTripleContains(bracket,b);
        assert a == bpTripleBase(bracket);
        assert b == add(a,bpextgen[2,2,2,8,0]);
        bpTripleContains(bracket,bpChosenHigher(x,y,x));
    ";
    assert_eq!(
        stdout(interpret(program)),
        "bpTripleContains = yes\nbpTripleContains = yes\nbpTripleContains = yes\n"
    );
}

#[test]
fn higher_checks_witness_order_shape_weight_and_chosen_obstructions() {
    let setup = "let x = bpextgen[3,2,1,4,0];\n";
    for (term, reason) in [
        ("bpChosenHigher(x,x)", "at least three"),
        (
            "bpHigherSystem[3](x,x,x,zero[BP(3,2,1)])",
            "proper witnesses",
        ),
        (
            "bpHigherSystem[3](x,x,x,zero[BP(3,2,2)],zero[BP(3,2,1)])",
            "BP(3, 2, 1)",
        ),
        (
            "bpHigherSystem[3](x,x,x,bpt[3,2,1,1,1],zero[BP(3,2,1)])",
            "homogeneous internal degree",
        ),
        (
            "bpHigherSystem[3](x,x,x,bpPow[2](bpv[3,2,1,1]),zero[BP(3,2,1)])",
            "normalized cochain",
        ),
        (
            "bpHigherSystem[3](x,x,x,zero[BP(3,2,1)],zero[BP(3,2,1)])",
            "wrong differential on interval [0,2)",
        ),
        (
            "bpChosenHigher(x,x,x,x)",
            "obstruction in its chosen system on interval [0,3)",
        ),
        (
            "bpChosenHigher(zero[BPExt(3,2,0,0)],x,x)",
            "positive cobar degrees",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
}
