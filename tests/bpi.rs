mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn invariant_ideal_reductions_commute_with_cofaces_and_preserve_arithmetic() {
    let program = r"
        let v = bpv[3, 3, 2, 1];
        let u = bpv[3, 3, 2, 2];
        let w = bpv[3, 3, 2, 3];
        let t = bpt[3, 3, 2, 2, 2];
        let f = add(bpMul(v,t), add(bpMul(w,t), u));
        let q = bpReduce[2](f);
        assert q == bpReduce[2](add(bpMul(w,t),u));
        assert bpReduce[2](q) == q;
        assert bpReduce[3](q) == bpReduce[3](bpMul(w,t));
        assert bpReduce[2](bpLift(q)) == q;
        assert bpReduce[2](bpCoface[0](f)) == bpCoface[0](q);
        assert bpReduce[2](bpCoface[1](f)) == bpCoface[1](q);
        assert bpReduce[2](bpCoface[2](f)) == bpCoface[2](q);
        assert bpReduce[2](bpCoface[3](f)) == bpCoface[3](q);
        assert bpReduce[2](bpD(f)) == bpD(q);
        assert bpD(bpD(q)) == zero[BPI(3, 3, 4, 2)];
        assert add(q, bpNeg(q)) == zero[BPI(3, 3, 2, 2)];
        assert bpMul(q, q) == bpPow[2](q);
        let two = bpReduce[1](bpint[3, 0, 0, -1]);
        assert bpPow[1000](two) == bpReduce[1](bpint[3, 0, 0, 1]);
        let top = bpReduce[4](bpv[4294967291, 4, 0, 4]);
        assert bpD(top) == zero[BPI(4294967291, 4, 1, 4)];
        show two;
        show bpLift(two);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "BPI(3, 0, 0, 1){2}\nBP(3, 0, 0){2}\n"
    );
}

#[test]
fn beta_and_alpha_have_the_specified_integral_lifts() {
    assert_eq!(
        stdout(example("bp-greek")),
        concat!(
            "BPI(3, 2, 1, 1){2*v1^2*t1_1 + t1_1^3}\n",
            "BP(3, 2, 2){4*v1*t1_1*t2_1 - t1_1*t2_1^2 + 5*t1_1^2*t2_1}\n",
            "BP(2, 1, 1){4*v1*t1_1^3 + 3*v1^2*t1_1^2 + v1^3*t1_1 + 2*t1_1^4}\n",
        )
    );
}

#[test]
fn beta_is_not_a_weight_twelve_boundary() {
    // These monomials span normalized one-cochains in weight twelve.
    // In their images the coefficient of t1_1^2*t2_1 is divisible by 3.
    let program = r"
        let v = bpv[3, 2, 1, 1];
        let t = bpt[3, 2, 1, 1, 1];
        show bpD(bpMul(bpPow[2](v),t));
        show bpD(bpMul(v,bpPow[2](t)));
        show bpD(bpPow[3](t));
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "BP(3, 2, 2){6*v1*t1_1*t2_1 + 9*t1_1^2*t2_1}\n",
            "BP(3, 2, 2){-2*v1*t1_1*t2_1 + 3*t1_1*t2_1^2}\n",
            "BP(3, 2, 2){-3*t1_1*t2_1^2 - 3*t1_1^2*t2_1}\n",
        )
    );
}

#[test]
fn canonical_higher_lifts_differ_between_two_and_three() {
    let program = r"
        assert bpDelta(bpDelta(bpReduce[2](bpv[2, 2, 0, 2]))) == zero[BP(2, 2, 2)];
        let gamma2 = bpDelta(bpDelta(bpDelta(bpReduce[3](bpv[2, 3, 0, 3]))));
        assert gamma2 == zero[BP(2, 3, 3)];
        let gamma3 = bpDelta(bpDelta(bpDelta(bpReduce[3](bpv[3, 3, 0, 3]))));
        assert bpD(gamma3) == zero[BP(3, 3, 4)];
        show gamma3;
    ";
    let output = stdout(interpret(program));
    let terms = output.strip_prefix("BP(3, 3, 3){").unwrap();
    assert_eq!(
        terms.split(" + ").flat_map(|sum| sum.split(" - ")).count(),
        27
    );
    assert!(terms.starts_with("4*v1*v2*t1_1*t2_1*t3_1"));
}

#[test]
fn connecting_maps_check_the_source_and_exact_divisibility() {
    for (program, message) in [
        ("show zero[BPI(3, 2, 0, 0)];", "ideal depth"),
        ("show zero[BPI(3, 2, 0, 4)];", "ideal depth"),
        (
            "show bpReduce[1](zero[BPI(3, 2, 0, 2)]);",
            "cannot decrease",
        ),
        ("show bpLift(zero[BP(3, 2, 0)]);", "requires a BPI cochain"),
        ("show bpDelta(zero[BP(3, 2, 0)]);", "requires a BPI cochain"),
        (
            "show bpDelta(bpReduce[2](bpt[3, 2, 1, 1, 2]));",
            "requires a closed BPI",
        ),
        (
            "show bpDelta(bpReduce[1](bpt[3, 2, 1, 1, 2]));",
            "requires a closed BPI",
        ),
        (
            "show bpBockP[1](zero[BPI(3, 2, 0, 1)]);",
            "requires an integral BP",
        ),
        (
            "show bpBockP[4](bpPow[4](bpv[2, 1, 0, 1]));",
            "divisible by p^4",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(message), "{program}\n{error}");
    }
    let exponent = usize::MAX;
    let source = format!("assert bpBockP[{exponent}](zero[BP(2, 0, 0)]) == zero[BP(2, 0, 1)];");
    assert_eq!(stdout(interpret(&source)), "");
    let program = format!("show bpBockP[{exponent}](bpv[2, 1, 0, 1]);");
    assert!(diagnostic(interpret(&program)).contains("differential to be divisible"));
    let coefficient = i128::MIN;
    let program =
        format!("assert bpBockP[127](bpint[2, 0, 1, {coefficient}]) == bpint[2, 0, 2, -1];");
    assert_eq!(stdout(interpret(&program)), "");
}
