mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn apery_franel_and_alternating_telescopers() {
    let line = "rvWZCheck = certified (2 initial coefficients checked)\n";
    assert_eq!(stdout(example("rv-wz-period")), line.repeat(3));
}

#[test]
fn literal_zero_denominator_count_and_negative_base() {
    assert_eq!(
        stdout(interpret(
            r"
            let z = add(rvOre[0,1],add(rvOre[1,1],rvOre[1,0]));
            rvWZCheck[1,0,-2,1](z) {{(0,1,-1)},{}};
            rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(0,1,-1),(6,0,4),(6,0,-4)},{}};
            "
        )),
        "rvWZCheck = certified (1 initial coefficient checked)\n".repeat(2)
    );
}

#[test]
fn a_true_tail_identity_does_not_hide_an_inhomogeneous_initial_term() {
    let output = diagnostic(interpret(
        "rvWZCheck[1,0,-1,1](rvOre[0,0]) {{(0,1,-1)},{(1,0,0,1)}};",
    ));
    assert!(
        output.contains("initial coefficient of L(sum) at q^0 is 1"),
        "{output}"
    );
}

#[test]
fn malformed_or_false_certificates_fail_even_with_the_right_low_coefficients() {
    for (input, reason) in [
        (
            "rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(0,1,-1),(6,13,1)}, {}};",
            "cleared bivariate identity fails",
        ),
        (
            "rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(0,0,1),(0,1,-1)}, {}};",
            "R(n,0) = 0",
        ),
        (
            "rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(0,1,-1)},{(1,-1,1,1)}};",
            "positive on n>=n0, 0<=k<=n+1",
        ),
        (
            "rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(0,1,-1)},{(1,-2,100,1)}};",
            "positive on n>=n0, 0<=k<=n+1",
        ),
        (
            "rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(0,1,-1)},{(-1,2,100,1)}};",
            "positive on n>=n0, 0<=k<=n+1",
        ),
        (
            "rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(0,1,-1)},{(1,1,100,0)}};",
            "denominator multiplicity must be positive",
        ),
        (
            "rvWZCheck[1,0,-1,1](rvOre[2,1]) {{},{}};",
            "largest q-shift 2",
        ),
        (
            "rvWZCheck[1,0,-1,0](rvOre[0,1]) {{},{}};",
            "start must be positive",
        ),
        (
            "rvWZCheck[0,1,1,1](rvOre[0,1]) {{},{}};",
            "requires u >= 1 and nonzero z",
        ),
        (
            "rvWZCheck[1,0,0,1](rvOre[0,1]) {{},{}};",
            "requires u >= 1 and nonzero z",
        ),
        (
            "rvWZCheck[1,0,1,1](S.hi) {{},{}};",
            "requires an Ore operator",
        ),
    ] {
        let output = diagnostic(interpret(input));
        assert!(output.contains(reason), "{input}\n{output}");
    }
}

#[test]
fn a_positive_negative_slope_factor_is_allowed_when_both_endpoints_are_safe() {
    // Multiply R=-k by (n-k+2)/(n-k+2) to exercise a nonconstant,
    // decreasing-in-k denominator that remains positive at k=n+1.
    assert_eq!(
        stdout(interpret(
            "rvWZCheck[1,0,-1,1](rvOre[0,1]) {{(1,1,-1),(0,2,1),(0,1,-2)},{(1,-1,2,1)}};"
        )),
        "rvWZCheck = certified (1 initial coefficient checked)\n",
    );
}

#[test]
fn all_initial_coefficients_are_checked_at_a_later_start() {
    let source = include_str!("../examples/rv-wz-period.asl")
        .replace("rvWZCheck[2,1,1,2]", "rvWZCheck[2,1,1,8]")
        .replace("rvWZCheck[2,0,-1,2]", "rvWZCheck[2,0,-1,7]");
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvWZCheck = certified (8 initial coefficients checked)\n",
            "rvWZCheck = certified (2 initial coefficients checked)\n",
            "rvWZCheck = certified (7 initial coefficients checked)\n",
        )
    );
}

#[test]
fn an_excessive_initial_range_reports_the_computation_limit() {
    let output = diagnostic(interpret(
        "rvWZCheck[1,0,-1,1000000000](rvOre[0,1]) {{(0,1,-1)},{}};",
    ));
    assert!(
        output.contains("statement exceeds the computation limit"),
        "{output}"
    );
}

#[test]
fn larger_rising_binomial_exponent_and_signed_apery_variant() {
    let input = r"
        let signed = add(rvOre[0,2],add(rvOre[1,2,11],add(rvOre[1,1,11],
            add(rvOre[1,0,3],add(rvOre[2,2,-1],add(rvOre[2,1,-2],rvOre[2,0,-1]))))));
        rvWZCheck[1,2,-1,6](signed) {
            {(3,3,-11),(2,4,-6),(2,3,15),(1,5,1),(1,4,1),(1,3,-4)},
            {(1,1,0,2),(1,1,-1,2)}
        };
    ";
    assert_eq!(
        stdout(interpret(input)),
        "rvWZCheck = certified (6 initial coefficients checked)\n"
    );
}
