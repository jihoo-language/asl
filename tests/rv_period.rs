mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const APERY: &str = r"
    rvperiod A[2] {
        ([-1,-1],1), ([-1,0],2), ([-1,1],1), ([0,-1],2),
        ([0,0],3), ([0,1],1), ([1,-1],1), ([1,0],1)
    };
";

#[test]
fn actual_section_states_retain_the_carries_of_the_laurent_period() {
    assert_eq!(
        stdout(example("rv-laurent-period")),
        concat!(
            "rvCTExact[Apery, 0] = 1\n",
            "rvCTExact[Apery, 1] = 3\n",
            "rvCTExact[Apery, 2] = 19\n",
            "rvCTExact[Apery, 3] = 147\n",
            "RVCTResidue(3, 2){6}\n",
            "RVCTResidue(3, 3){19}\n",
            "RVCTResidue(2, 4){3}\n",
            "RVCTState(Carry, 2, 1){1*z^0[0] + 1*z^0[1]}\n",
            "rvCTExact[Carry, 3] = 4\n",
        )
    );
}

fn binomial(n: u128, k: u128) -> u128 {
    (0..k).fold(1, |x, j| x * (n - j) / (j + 1))
}

fn apery(n: u128) -> u128 {
    (0..=n)
        .map(|k| binomial(n, k).pow(2) * binomial(n + k, k))
        .sum()
}

#[test]
fn multivariate_and_signed_periods_match_exact_powers_and_prime_power_reductions() {
    let mut source = APERY.to_string();
    source.push_str(
        r"
        rvperiod B[1] {([0],2),([1],1),([-1],2),([-1],-1)};
        rvperiod Signed[1] {([0],-2),([1],1),([-1],-1)};
        rvperiod Empty[2] {};
        rvCTExact[Signed,1];
        rvCTExact[Signed,2];
        rvCTExact[Signed,3];
        rvCTExact[Empty,0];
        rvCTExact[Empty,4];
        assert rvCTMod[Empty,3,3,99] == rvCTResidue[3,3,0];
    ",
    );
    for n in 0..=10 {
        writeln!(source, "rvCTExact[A,{n}];").unwrap();
        for (p, s) in [(2u128, 1), (2, 3), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2)] {
            let modulus = p.pow(s);
            let expected = apery(n) % modulus;
            let central = binomial(2 * n, n) % modulus;
            writeln!(
                source,
                "assert rvCTMod[A,{p},{s},{n}] == rvCTResidue[{p},{s},{expected}];"
            )
            .unwrap();
            writeln!(
                source,
                "assert rvCTMod[B,{p},{s},{n}] == rvCTResidue[{p},{s},{central}];"
            )
            .unwrap();
        }
    }
    let mut expected = "rvCTExact[Signed, 1] = -2\nrvCTExact[Signed, 2] = 2\nrvCTExact[Signed, 3] = 4\nrvCTExact[Empty, 0] = 1\nrvCTExact[Empty, 4] = 0\n".to_string();
    for n in 0..=10 {
        writeln!(expected, "rvCTExact[A, {n}] = {}", apery(n)).unwrap();
    }
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn section_composition_and_linearity_are_native_first_class_operations() {
    let source = format!(
        "{APERY}{}",
        r"
        let init = rvCTInitial[A,3,2];
        let one = rvCTSection[1](init);
        let two = rvCTSection[2](init);
        let section12 = rvCTSection[2](one);
        let section21 = rvCTSection[1](two);
        assert rvCTRead(section12) == rvCTMod[A,3,2,7];
        assert rvCTRead(section21) == rvCTMod[A,3,2,5];
        assert rvCTSection[2](add(init,init)) == add(two,two);
        assert rvCTRead(add(one,two)) == add(rvCTRead(one),rvCTRead(two));
        assert rvCTSection[1](zero[RVCTState(A,3,2)]) == zero[RVCTState(A,3,2)];
        assert rvCTRead(zero[RVCTState(A,3,2)]) == zero[RVCTResidue(3,2)];
        assert rvCTMod[A,3,2,1350851717672992088] == rvCTResidue[3,2,1];
        rvperiod Unit[0] {([],1)};
        let step = \x : RVCTState(Unit,2,1) => rvCTSection[1](x);
        let read = \x : RVCTState(Unit,2,1) => rvCTRead(x);
        assert read(step(step(rvCTInitial[Unit,2,1]))) == rvCTResidue[2,1,1];
        card (\x : RVCTState(Unit,2,1) => S.hi);
        card (\x : RVCTResidue(2,2) => S.hi);
        assert rvCTMod[Unit,2,1,340282366920938463463374607431768211455] == rvCTResidue[2,1,1];
    "
    );
    assert_eq!(stdout(interpret(&source)), "card = 2 / 2\ncard = 4 / 4\n");
}

#[test]
fn a_wide_prime_power_uses_checked_modular_arithmetic() {
    let large = 1u128 << 126;
    let modulus = 2 * large;
    let source = format!(
        r"
        rvperiod Nil[0] {{([],{large})}};
        assert rvCTMod[Nil,2,127,1] == rvCTResidue[2,127,{large}];
        assert rvCTMod[Nil,2,127,2] == rvCTResidue[2,127,0];
        assert add(rvCTResidue[2,127,{large}],rvCTResidue[2,127,{large}]) == zero[RVCTResidue(2,127)];
        assert rvCTResidue[2,127,{}] == rvCTResidue[2,127,{}];
        rvperiod Min[0] {{([],{})}};
        rvCTExact[Min,1];
        assert rvCTMod[Min,3,1,1] == rvCTResidue[3,1,1];
        ",
        u128::MAX,
        modulus - 1,
        i128::MIN
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!("rvCTExact[Min, 1] = {}\n", i128::MIN)
    );
}

#[test]
fn malformed_periods_digits_and_storage_have_distinct_diagnostics() {
    for (statement, expected) in [
        ("rvperiod Bad[2]{([1],2)};", "requires 2 spatial exponents"),
        ("rvperiod A[0]{};", "already defined"),
        ("rvCTExact[Unknown,0];", "unknown Laurent period"),
        ("show rvCTInitial[A,3,0];", "precision must be positive"),
        ("show rvCTMod[A,9,1,1];", "9 is not prime"),
        (
            "show rvCTInitial[A,2,128];",
            "prime-power modulus exceeds u128",
        ),
        (
            "show rvCTSection[3](rvCTInitial[A,3,1]);",
            "digit 3 is outside 0..3",
        ),
        (
            "show rvCTRead(rvCTResidue[3,1,0]);",
            "requires an RVCTState",
        ),
        (
            "show rvCTMod[A,3,1,340282366920938463463374607431768211456];",
            "period index exceeds u128",
        ),
        (
            "card (\\k : RVCTState(A,3,2) => S.hi);",
            "Cartier state exceeds the finite type enumeration limit",
        ),
        (
            "rvperiod Big[0]{([],170141183460469231731687303715884105727)}; rvCTExact[Big,2];",
            "exact Laurent period coefficient exceeds i128",
        ),
        (
            "rvperiod Big[1]{([-170141183460469231731687303715884105728],1)}; rvCTExact[Big,2];",
            "Laurent period exponent is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{APERY}{statement}")));
        assert!(error.contains(expected), "{statement}\n{error}");
    }
}

#[test]
fn exact_apery_terms_agree_with_the_grassmannian_normalization() {
    let mut source = APERY.to_string();
    source.push_str("rvQGrPeriod[2,5,7];\n");
    for n in 0..=7 {
        writeln!(source, "rvCTExact[A,{n}];").unwrap();
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    let fractions = lines
        .next()
        .unwrap()
        .strip_prefix("rvQGrPeriod = [")
        .unwrap()
        .strip_suffix(']')
        .unwrap();
    for (n, (fraction, period)) in fractions.split(", ").zip(lines).enumerate() {
        let (p, q) = fraction.split_once('/').unwrap_or((fraction, "1"));
        let (p, q): (u128, u128) = (p.parse().unwrap(), q.parse().unwrap());
        let value: u128 = period.rsplit_once(" = ").unwrap().1.parse().unwrap();
        let factorial = (1..=n as u128).product::<u128>();
        assert_eq!(value * q, p * factorial.pow(5));
    }
}
