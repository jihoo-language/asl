mod common;

use common::{diagnostic, example, interpret, stdout};

fn terms(line: &str) -> usize {
    line.split(" + ").flat_map(|sum| sum.split(" - ")).count()
}

#[test]
fn divided_greek_letter_cocycles_have_the_expected_weights_and_orders() {
    let output = stdout(example("bp-divided-greek"));
    let rows: Vec<_> = output.lines().collect();
    assert_eq!(terms(rows[0]), 5);
    assert_eq!(&rows[1..], ["bporder[8] = 2", "bporder[36] = 3"]);
    for (p, s, j, weight, count, order) in [
        (2, 4, 4, 16, 13, Some(2)),
        (3, 3, 3, 36, 15, Some(3)),
        (3, 3, 2, 40, 17, Some(3)),
        (5, 5, 5, 200, 47, None),
    ] {
        let mut program = format!(
            "let x = bpThick[{j}](bpPow[{s}](bpReduce[1](bpv[{p},2,0,2])));\n\
             assert bpD(x) == zero[BPThick({p},2,1,2,{j})];\n\
             let intermediate = bpThickDelta(x);\n\
             assert bpD(intermediate) == zero[BPI({p},2,2,1)];\n\
             let beta = bpDelta(intermediate);\n\
             assert bpD(beta) == zero[BP({p},2,3)];\nshow beta;\n"
        );
        if order.is_some() {
            program.push_str(&format!("bporder[{weight}](beta);\n"));
        }
        let text = stdout(interpret(&program));
        let mut rows = text.lines();
        assert_eq!(terms(rows.next().unwrap()), count);
        if let Some(order) = order {
            assert_eq!(rows.next().unwrap(), format!("bporder[{weight}] = {order}"));
        }
    }
}

#[test]
fn quotient_arithmetic_keeps_the_correct_nilpotent_and_arrow_terms() {
    let program = r"
        let v = bpReduce[1](bpv[3,2,0,1]);
        let a = bpThick[2](v);
        let b = bpThick[3](v);
        assert bpThickLift(a) == v;
        assert bpPow[2](a) == zero[BPThick(3,2,0,2,2)];
        assert bpThickLift(bpPow[2](b)) == bpPow[2](v);
        assert bpPow[3](b) == zero[BPThick(3,2,0,2,3)];
        assert add(a,bpNeg(a)) == zero[BPThick(3,2,0,2,2)];
        let z = bpReduce[1](bpv[3,2,0,2]);
        let t = bpReduce[1](bpt[3,2,1,1,2]);
        assert bpThickLift(bpPow[3](bpThick[2](z))) == bpPow[3](z);
        assert bpThickLift(bpPow[3](bpThick[2](t))) == bpPow[3](t);
        let x = bpThick[2](t);
        let y = bpThick[2](z);
        assert bpCup(x,y) == bpThick[2](bpCup(t,z));
        assert bpD(bpCup(x,y)) == bpThick[2](bpD(bpCup(t,z)));
        assert bpD(bpD(bpCup(x,y))) == zero[BPThick(3,2,3,2,2)];
        assert bpCoface[0](y) == bpThick[2](bpCoface[0](z));
        assert bpCoface[1](x) == bpThick[2](bpCoface[1](t));
        assert bpDegeneracy[0](x) == zero[BPThick(3,2,0,2,2)];
        assert bpRaise[3](x) == bpThick[2](bpRaise[3](t));
    ";
    assert_eq!(stdout(interpret(program)), "");
    let j = usize::MAX;
    let program = format!(
        "let v=bpThick[{j}](bpReduce[1](bpv[2,1,0,1]));\n\
         let x=bpPow[{}](v);\n\
         assert bpMul(x,x) == zero[BPThick(2,1,0,2,{j})];\n\
         assert bpPow[{j}](v) == zero[BPThick(2,1,0,2,{j})];",
        j - 1
    );
    assert_eq!(stdout(interpret(&program)), "");
}

#[test]
fn order_one_agrees_with_bpi_and_change_of_lift_changes_by_a_boundary() {
    let program = r"
        let z3 = bpReduce[1](bpv[3,2,0,2]);
        let one = bpThick[1](z3);
        assert bpThickDelta(one) == bpDelta(bpReduce[2](one));
        assert bpReduce[2](one) == bpReduce[2](z3);
        let u = bpReduce[1](bpMul(bpv[2,2,0,1],bpv[2,2,0,2]));
        let v2 = bpReduce[1](bpv[2,2,0,2]);
        let v1 = bpReduce[1](bpv[2,2,0,1]);
        let square = bpPow[2](v1);
        let x = bpPow[2](v2);
        let alternate = add(x,bpMul(square,u));
        let thick = bpThick[2](x);
        assert bpThick[2](alternate) == thick;
        let canonical_delta = bpThickDelta(thick);
        let alternate_delta = add(canonical_delta,bpD(u));
        assert bpD(alternate) == bpMul(bpCoface[1](square),alternate_delta);
        assert bpExtClass[8](bpDelta(alternate_delta)) == bpExtClass[8](bpDelta(canonical_delta));
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn height_three_divides_by_the_last_invariant_v2_power() {
    let program = r"
        let v3sq = bpPow[2](bpReduce[2](bpv[2,3,0,3]));
        let x = bpThick[2](v3sq);
        assert bpD(x) == zero[BPThick(2,3,1,3,2)];
        let t = bpReduce[2](bpt[2,3,1,1,1]);
        let v2 = bpReduce[2](bpv[2,3,1,2]);
        let y = add(bpMul(bpPow[2](v2),bpPow[2](t)),bpPow[8](t));
        assert bpThickDelta(x) == y;
        assert bpD(v3sq) == bpMul(bpPow[2](v2),y);
        assert bpD(y) == zero[BPI(2,3,2,2)];
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn exact_v1_divisibility_detects_the_first_noncycle_thickening() {
    for (p, s, j) in [
        (2, 2, 2),
        (2, 4, 4),
        (2, 6, 2),
        (3, 2, 1),
        (3, 3, 3),
        (3, 6, 3),
        (5, 5, 5),
    ] {
        let v = format!("bpPow[{s}](bpReduce[1](bpv[{p},2,0,2]))");
        let valid = format!("assert bpD(bpThick[{j}]({v})) == zero[BPThick({p},2,1,2,{j})];");
        assert_eq!(stdout(interpret(&valid)), "");
        let invalid = format!("show bpThickDelta(bpThick[{}]({v}));", j + 1);
        assert!(diagnostic(interpret(&invalid)).contains("requires a closed BPThick"));
    }
    for (term, expected) in [
        (
            "bpThick[0](zero[BPI(3,2,0,1)])",
            "exponent must be positive",
        ),
        ("zero[BPThick(3,2,0,1,2)]", "2 <= h <= n+1"),
        ("bpThick[2](zero[BPI(3,2,0,3)])", "last invariant generator"),
        ("bpThick[2](zero[BP(3,2,0)])", "requires a BPI cochain"),
        (
            "bpThickDelta(zero[BPI(3,2,0,1)])",
            "requires a BPThick cochain",
        ),
        (
            "bpThickLift(zero[BPI(3,2,0,1)])",
            "requires a BPThick cochain",
        ),
        (
            "bpReduce[1](zero[BPThick(3,2,0,2,2)])",
            "must kill the thickened generator",
        ),
        (
            "bpDelta(zero[BPThick(3,2,0,2,2)])",
            "requires a BPI cochain",
        ),
        (
            "bpCup(zero[BPThick(3,2,0,2,2)],zero[BPThick(3,2,0,2,3)])",
            "same truncation",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(expected), "{term}\n{error}");
    }
}
