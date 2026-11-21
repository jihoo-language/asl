mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

const BASIS: &str = "
    let x = extClass(pword[3]{([1], [])});
    let a = extClass(pword[3]{([], [0])});
";

#[test]
fn the_first_nonzero_power_has_the_suspension_sign_at_each_odd_prime() {
    assert_eq!(
        stdout(example("odd-transfer")),
        "PExt(3, 2, 12){2*g0}\nPExt(5, 2, 40){g0}\n",
    );
    for p in [3, 5, 7] {
        let degree = 2 * (p - 1);
        let mut source = format!("let x = extClass(pword[{p}]{{([1], [])}});\n");
        let mut factorial = 1;
        for n in 2..=p {
            let t = n * degree;
            let args = vec!["x"; n].join(", ");
            let negative = n * (n - 1) / 2 % 2 == 1;
            if n == p {
                let sign = if negative { p - 1 } else { 1 };
                writeln!(
                    source,
                    "assert ainfty({args}) == pscale[{sign}](pextgen[{p}, 2, {t}, 0]);"
                )
                .unwrap();
            } else {
                factorial = factorial * n % p;
                let sign = if negative { p - 1 } else { 1 };
                let coefficient = (1..p).find(|k| k * factorial % p == sign).unwrap();
                writeln!(
                    source,
                    "assert ainfty({args}) == zero[PExt({p}, 2, {t})];
                     assert alift({args}) == pscale[{coefficient}](pword[{p}]{{([{n}], [])}});"
                )
                .unwrap();
            }
        }
        assert!(stdout(interpret(&source)).is_empty());
    }
}

#[test]
fn binary_multiplication_scalar_units_and_multilinearity_keep_field_coefficients() {
    let source = format!(
        "{BASIS}
        let y = extClass(pword[3]{{([9], [])}});
        assert ainfty(x, y) == yoneda(x, y);
        assert ainfty(y, x) == pscale[2](ainfty(x, y));
        let u = pscale[2](pextgen[3, 0, 0, 0]);
        assert ainfty(u, x) == pscale[2](x);
        assert ainfty(x, u) == pscale[2](x);
        assert alift(u, x) == zero[PCobar(3, 0, 4)];
        assert alift(x, u, x) == zero[PCobar(3, 0, 8)];
        assert ainfty(x, u, x) == zero[PExt(3, 1, 8)];
        let m = ainfty(x, a, x);
        assert m == pscale[2](pextgen[3, 2, 9, 0]);
        assert ainfty(pscale[2](x), pscale[2](a), pscale[2](x)) == pscale[2](m);
        assert ainfty(add(x, x), a, x) == add(m, m);
        assert alift(x, x, a) == pword[3]{{([1], [1])}};
        assert alift(pscale[2](x), x, a) == pscale[2](alift(x, x, a));
        assert cobarD(alift(x, y)) == add(extRep(ainfty(x, y)),
                                         pscale[2](concat(extRep(x), extRep(y))));
    "
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn signed_stasheff_and_morphism_equations_hold_through_arity_five() {
    let x = Expression::new("x", 1, 4);
    let a = Expression::new("a", 1, 1);
    let mut program = format!(
        "{BASIS}
        assert ainfty(a, ainfty(x, x, x)) == pscale[2](pextgen[3, 3, 13, 0]);
        assert ainfty(ainfty(a, x, x), x) == pextgen[3, 3, 13, 0];
        assert ainfty(x, x, a, x) == pscale[2](pextgen[3, 2, 13, 0]);
        assert alift(x, x, a, x)
            == add(pword[3]{{([3], [0])}}, pscale[2](pword[3]{{([2], [1])}}));
        assert alift(x, x, a, x, x) == pscale[2](pword[3]{{([0, 1], [0])}});
    "
    );
    identities(
        &mut program,
        3,
        &[a.clone(), x.clone(), x.clone(), x.clone()],
    );
    identities(
        &mut program,
        3,
        &[x.clone(), x.clone(), a.clone(), x.clone()],
    );
    identities(&mut program, 3, &[x.clone(), x.clone(), a, x.clone(), x]);
    assert!(stdout(interpret(&program)).is_empty());
}

#[derive(Clone)]
struct Expression {
    text: String,
    s: usize,
    t: usize,
}

impl Expression {
    fn new(text: impl Into<String>, s: usize, t: usize) -> Self {
        Self {
            text: text.into(),
            s,
            t,
        }
    }

    fn odd(&self) -> bool {
        (self.t - self.s) % 2 == 1
    }
}

fn signed(text: String, p: usize, negative: bool) -> String {
    if negative {
        format!("pscale[{}]({text})", p - 1)
    } else {
        text
    }
}

fn suspended(operation: &str, p: usize, inputs: &[Expression]) -> Expression {
    let n = inputs.len();
    if n == 1 {
        assert_eq!(operation, "alift");
        let input = &inputs[0];
        return Expression::new(format!("extRep({})", input.text), input.s, input.t);
    }
    let text = inputs
        .iter()
        .map(|input| input.text.as_str())
        .collect::<Vec<_>>()
        .join(", ");
    let sign: usize = inputs
        .iter()
        .enumerate()
        .map(|(i, input)| (n - i - 1) * usize::from(input.odd()))
        .sum();
    let s = inputs.iter().map(|input| input.s).sum::<usize>() - n
        + if operation == "alift" { 1 } else { 2 };
    let t = inputs.iter().map(|input| input.t).sum();
    Expression::new(
        signed(format!("{operation}({text})"), p, sign % 2 == 1),
        s,
        t,
    )
}

fn sum(terms: Vec<String>) -> String {
    terms
        .into_iter()
        .reduce(|left, right| format!("add({left}, {right})"))
        .unwrap()
}

fn identities(source: &mut String, p: usize, inputs: &[Expression]) {
    let n = inputs.len();
    let mut algebra = Vec::new();
    let mut comparison = Vec::new();
    for k in 2..=n {
        for start in 0..=n - k {
            let inserted = suspended("ainfty", p, &inputs[start..start + k]);
            let args: Vec<_> = inputs[..start]
                .iter()
                .cloned()
                .chain([inserted])
                .chain(inputs[start + k..].iter().cloned())
                .collect();
            let negative = inputs[..start].iter().filter(|input| !input.odd()).count() % 2 == 1;
            comparison.push(signed(suspended("alift", p, &args).text, p, negative));
            if k < n {
                algebra.push(signed(suspended("ainfty", p, &args).text, p, negative));
            }
        }
    }
    let full = suspended("alift", p, inputs);
    let mut direct = vec![format!("cobarD({})", full.text)];
    for cut in 1..n {
        let left = suspended("alift", p, &inputs[..cut]);
        let right = suspended("alift", p, &inputs[cut..]);
        direct.push(signed(
            format!("concat({}, {})", left.text, right.text),
            p,
            left.odd(),
        ));
    }
    writeln!(
        source,
        "assert {} == zero[PExt({p}, {}, {})];\nassert {} == {};",
        sum(algebra),
        full.s + 2,
        full.t,
        sum(direct),
        sum(comparison)
    )
    .unwrap();
}

#[test]
fn odd_transfer_checks_primes_grades_and_lift_table_size() {
    for (expression, expected) in [
        ("ainfty(x)".into(), "at least two"),
        ("ainfty(x, pextgen[5, 1, 8, 0])".into(), "same prime"),
        ("alift(x, extgen[1, 1, 0])".into(), "a PExt class"),
        ("alift(x, pword[3]{([1], [])})".into(), "a PExt class"),
        ("ainfty(extgen[1, 1, 0], x)".into(), "an Ext class"),
        (
            "alift(pextgen[3, 0, 0, 0], pextgen[3, 0, 0, 0])".into(),
            "negative",
        ),
        (
            format!("ainfty(zero[PExt(3, 1, {})], x)", usize::MAX),
            "internal degree",
        ),
        (
            format!("ainfty(zero[PExt(3, {}, 0)], x)", usize::MAX),
            "cohomological degree",
        ),
        (
            format!("ainfty({})", vec!["x"; 1000].join(",")),
            "lift table exceeds",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{BASIS}show {expression};")));
        assert!(error.contains(expected), "{expression}\n{error}");
    }
    let source = format!(
        "{BASIS}assert alift(zero[PExt(3, {0}, 0)], x) == zero[PCobar(3, {0}, 4)];",
        usize::MAX,
    );
    assert!(stdout(interpret(&source)).is_empty());
}
