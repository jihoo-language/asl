mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn chosen_triple_massey_value_is_the_negative_of_the_greek_class() {
    assert_eq!(
        stdout(example("bp-massey")),
        "BP(3, 2, 2){2*v1*t1_1*t2_1 + t1_1*t2_1^2 + 4*t1_1^2*t2_1}\n",
    );
}

#[test]
fn cup_uses_ordered_paths_and_is_a_signed_differential_algebra_product() {
    let program = r"
        let x = bpt[3, 2, 1, 1, 1];
        let y = bpt[3, 2, 1, 1, 2];
        let v = bpv[3, 2, 0, 1];
        let v1 = bpv[3, 2, 1, 1];
        let three = bpint[3, 2, 1, 3];
        assert bpCup(x,v) == bpMul(x,add(v1,bpMul(three,x)));
        assert bpCup(v,x) == bpMul(v1,x);
        assert bpCup(bpCup(x,y),v) == bpCup(x,bpCup(y,v));
        assert bpCup(bpCup(v,y),x) == bpCup(v,bpCup(y,x));
        assert bpCup(y,bpint[3, 2, 0, 1]) == y;
        assert bpCup(bpint[3, 2, 0, 1],y) == y;
        assert bpD(bpCup(y,v)) == add(bpCup(bpD(y),v),bpNeg(bpCup(y,bpD(v))));
        assert bpD(bpCup(v,y)) == add(bpCup(bpD(v),y),bpCup(v,bpD(y)));
        let yy = bpCup(y,y);
        assert bpD(bpCup(yy,v)) == add(bpCup(bpD(yy),v),bpCup(yy,bpD(v)));
        let z = bpCup(x,y);
        assert bpDegeneracy[0](z) == zero[BP(3, 2, 1)];
        assert bpDegeneracy[1](z) == zero[BP(3, 2, 1)];
        assert bpNormalize(z) == z;
        assert bpReduce[2](bpCup(y,v)) == bpCup(bpReduce[2](y),bpReduce[2](v));
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn normalization_uses_the_identity_arrow_and_commutes_with_differential() {
    let mut program = String::new();
    for (p, a) in [(2, 2), (3, 2), (5, 1), (2, 0)] {
        let mut term = format!("bpv[{p},2,{a},2]");
        for j in 1..=a {
            term = format!("add({term},bpt[{p},2,{a},{j},2])");
        }
        let f = format!("f_{p}_{a}");
        writeln!(program, "let {f} = {term};").unwrap();
        writeln!(
            program,
            "assert bpNormalize(bpNormalize({f})) == bpNormalize({f});"
        )
        .unwrap();
        writeln!(
            program,
            "assert bpD(bpNormalize({f})) == bpNormalize(bpD({f}));"
        )
        .unwrap();
        writeln!(
            program,
            "assert bpNormalize(bpReduce[2]({f})) == bpReduce[2](bpNormalize({f}));"
        )
        .unwrap();
        for j in 0..a {
            writeln!(
                program,
                "assert bpDegeneracy[{j}](bpNormalize({f})) == zero[BP({p}, 2, {})];",
                a - 1
            )
            .unwrap();
        }
        for j in 0..=a {
            for i in 0..=a + 1 {
                let right = if i < j {
                    format!("bpCoface[{i}](bpDegeneracy[{}]({f}))", j - 1)
                } else if i > j + 1 {
                    format!("bpCoface[{}](bpDegeneracy[{j}]({f}))", i - 1)
                } else {
                    f.clone()
                };
                writeln!(
                    program,
                    "assert bpDegeneracy[{j}](bpCoface[{i}]({f})) == {right};"
                )
                .unwrap();
            }
        }
        if a >= 2 {
            for j in 0..a - 1 {
                for i in 0..=j {
                    writeln!(program, "assert bpDegeneracy[{j}](bpDegeneracy[{i}]({f})) == bpDegeneracy[{i}](bpDegeneracy[{}]({f}));", j+1).unwrap();
                }
            }
        }
    }
    program.push_str("assert bpNormalize(bpt[2,2,2,2,2]) == bpD(bpt[2,2,1,1,2]);");
    assert_eq!(stdout(interpret(&program)), "");
}

#[test]
fn massey_signs_and_witnesses_work_in_integral_and_quotient_rings() {
    let program = r"
        let a = bpt[3, 2, 1, 1, 1];
        let v = bpv[3, 2, 1, 1];
        let w = add(bpMul(v,a),bpPow[2](a));
        let aa = bpCup(a,a);
        let u = bpCup(w,a);
        let even = bpMassey(aa,a,a,u,w);
        assert even == bpD(bpCup(w,w));
        let x = bpReduce[1](a);
        let y = bpReduce[1](w);
        assert bpMassey(x,x,x,y,y) == bpReduce[1](bpMassey(a,a,a,w,w));
        assert bpMassey(bpCup(x,x),x,x,bpCup(y,x),y) == bpReduce[1](even);
        show even;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "BP(3, 2, 3){t1_1*t2_1*t3_1^2 + 3*t1_1*t2_1^2*t3_1 + 2*t1_1^2*t2_1*t3_1}\n",
    );
}

#[test]
fn malformed_cup_degeneracies_and_defining_systems_get_diagnostics() {
    let setup = r"
        let a = bpt[3, 2, 1, 1, 1];
        let q = bpt[3, 2, 1, 1, 2];
        let v = bpv[3, 2, 1, 1];
        let w = add(bpMul(v,a),bpPow[2](a));
    ";
    for (term, message) in [
        (
            "bpCup(a,bpv[2,2,0,1])",
            "same prime, height, and ideal depth",
        ),
        (
            "bpCup(a,bpReduce[1](a))",
            "same prime, height, and ideal depth",
        ),
        ("bpCup(a,zero[BP(3,2,16)])", "arity exceeds the limit"),
        ("bpDegeneracy[1](a)", "degeneracy index"),
        (
            "bpDegeneracy[0](bpv[3,2,0,1])",
            "requires at least one arrow",
        ),
        ("bpMassey(a,a,a,w)", "requires five cochains"),
        ("bpMassey(q,a,a,w,w)", "input A must be closed"),
        (
            "bpMassey(a,a,a,zero[BP(3,2,1)],w)",
            "requires d(U) = bpCup(A, B)",
        ),
        (
            "bpMassey(a,a,a,w,zero[BP(3,2,1)])",
            "requires d(V) = bpCup(B, C)",
        ),
        (
            "bpMassey(a,a,a,zero[BP(3,2,2)],w)",
            "witness U must have arity",
        ),
        (
            "bpMassey(a,a,a,bpReduce[1](w),w)",
            "same prime, height, and ideal depth",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\nshow {term};")));
        assert!(error.contains(message), "{term}\n{error}");
    }
}
