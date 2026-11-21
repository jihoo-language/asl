mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn hazewinkel_right_unit_and_internal_composition_at_two() {
    assert_eq!(
        stdout(example("bp-cobar")),
        concat!(
            "BP(2, 2, 1){v1 + 2*t1_1}\n",
            "BP(2, 2, 1){-5*v1*t1_1^2 - 3*v1^2*t1_1 + v2 - 4*t1_1^3 + 2*t1_2}\n",
            "BP(2, 2, 2){-v1*t1_1*t2_1 + t1_1*t2_1^2 + t1_2 + t2_2}\n",
            "BP(2, 2, 2){v1*t1_1*t2_1 - t1_1*t2_1^2}\n",
        )
    );
}

#[test]
fn hazewinkel_cofaces_transport_the_base_of_later_arrows() {
    let program = r"
        let v = bpv[3, 2, 2, 1];
        let a = bpt[3, 2, 2, 1, 1];
        let b = bpt[3, 2, 2, 2, 1];
        let three = bpint[3, 2, 2, 3];
        assert bpCoface[1](bpt[3, 2, 1, 1, 1]) == add(a, b);
        assert bpCoface[1](bpCoface[0](bpv[3, 2, 0, 1])) == add(v, bpMul(three, add(a, b)));
        show bpCoface[0](bpv[3, 2, 0, 2]);
        show bpCoface[2](bpt[3, 2, 2, 2, 2]);
        show bpCoface[2](bpt[2, 2, 2, 2, 2]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "BP(3, 2, 1){-35*v1*t1_1^3 - 18*v1^2*t1_1^2 - 4*v1^3*t1_1 + v2 - 27*t1_1^4 + 3*t1_2}\n",
            "BP(3, 2, 3){-v1*t2_1*t3_1^2 - v1*t2_1^2*t3_1 - 3*t1_1*t2_1*t3_1^2 - 3*t1_1*t2_1^2*t3_1 + t2_1*t3_1^3 + t2_2 + t3_2}\n",
            "BP(2, 2, 3){-v1*t2_1*t3_1 - 2*t1_1*t2_1*t3_1 + t2_1*t3_1^2 + t2_2 + t3_2}\n",
        )
    );
}

#[test]
fn cofaces_are_ring_maps_and_satisfy_cosimplicial_identities() {
    let mut program = String::new();
    use std::fmt::Write;
    for (p, n, a) in [(2, 3, 0), (2, 3, 1), (2, 2, 2), (3, 2, 2), (5, 2, 1)] {
        let next = a + 2;
        let mut term = format!("bpv[{p},{n},{a},{n}]");
        for slot in 1..=a {
            term = format!("add({term},bpt[{p},{n},{a},{slot},{n}])");
        }
        writeln!(program, "let f_{p}_{n}_{a} = {term};").unwrap();
        let f = format!("f_{p}_{n}_{a}");
        for j in 1..=a + 2 {
            for i in 0..j {
                let k = j - 1;
                writeln!(program, "assert bpCoface[{j}](bpCoface[{i}]({f})) == bpCoface[{i}](bpCoface[{k}]({f}));").unwrap();
            }
        }
        for i in 0..=a + 1 {
            writeln!(
                program,
                "assert bpCoface[{i}](bpPow[2]({f})) == bpPow[2](bpCoface[{i}]({f}));"
            )
            .unwrap();
        }
        writeln!(program, "assert bpD(bpD({f})) == zero[BP({p},{n},{next})];").unwrap();
    }
    assert_eq!(stdout(interpret(&program)), "");
}

#[test]
fn bp_arithmetic_is_exact_and_checks_bounds() {
    assert_eq!(
        stdout(interpret(
            r"
            let v = bpv[3, 1, 0, 1];
            let a = add(v, bpint[3, 1, 0, -2]);
            assert add(a, bpNeg(a)) == zero[BP(3, 1, 0)];
            assert bpPow[2](a) == bpMul(a, a);
            assert bpPow[0](zero[BP(3, 1, 0)]) == bpint[3, 1, 0, 1];
            assert bpD(bpint[2, 0, 0, -170141183460469231731687303715884105728]) == zero[BP(2, 0, 1)];
            show bpPow[2](a);
        "
        )),
        "BP(3, 1, 0){4 - 4*v1 + v1^2}\n"
    );
    for (program, message) in [
        ("show bpv[4, 1, 0, 1];", "requires a prime"),
        ("show bpv[9, 1, 0, 1];", "is not prime"),
        ("show bpv[2, 9, 0, 1];", "height exceeds the limit"),
        ("show bpv[2, 1, 0, 0];", "generator index"),
        ("show bpt[2, 1, 1, 0, 1];", "arrow slots start"),
        ("show bpt[2, 1, 1, 2, 1];", "arrow slot must be"),
        ("show bpCoface[2](bpv[2, 1, 0, 1]);", "coface index"),
        (
            "show bpD(bpv[2, 1, 16, 1]);",
            "arrow arity exceeds the limit",
        ),
        ("show bpCoface[0](zero[Ext(0, 0)]);", "a BP polynomial"),
        (
            "show bpMul(bpv[2, 1, 0, 1], bpv[3, 1, 0, 1]);",
            "BP(2, 1, 0)",
        ),
        ("show bpPow[127](bpint[2, 0, 0, 2]);", "coefficient exceeds"),
        (
            "show bpNeg(bpint[2, 0, 0, -170141183460469231731687303715884105728]);",
            "coefficient exceeds",
        ),
        (
            "show bpint[2, 0, 0, 170141183460469231731687303715884105728];",
            "coefficient exceeds",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(message), "{program}\n{error}");
    }
    let exponent = usize::MAX;
    let source = format!("let v = bpv[2, 1, 0, 1];\nshow bpMul(bpPow[{exponent}](v), v);");
    assert!(diagnostic(interpret(&source)).contains("BP monomial exponent overflow"));
}
