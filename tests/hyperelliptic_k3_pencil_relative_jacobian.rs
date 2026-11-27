mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const BASIS: &str = r"
    let r0=rvJacPoint[C](z,z);
    let r1=rvJacPoint[C](o,z);
    let r2=rvJacPoint[C](two,z);
    let r3=rvJacPoint[C](three,z);
    rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
    rvk3pencil P[B,0,1];
";

fn split(prime: i64) -> String {
    format!(
        r"
        rvunram F[{prime},1]{{0,1}};
        let z=rvUnramElem[F]{{}};
        let o=rvUnramElem[F]{{1}};
        let two=rvUnramElem[F]{{2}};
        let three=rvUnramElem[F]{{3}};
        let c1=rvUnramElem[F]{{24}};
        let c2=rvUnramElem[F]{{-50}};
        let c3=rvUnramElem[F]{{35}};
        let c4=rvUnramElem[F]{{-10}};
        rvhyperelliptic C[F]{{z,c1,c2,c3,c4,o}};
        {BASIS}"
    )
}

#[test]
fn global_model_evaluates_to_each_individual_fiber_including_infinity() {
    for prime in [5, 13] {
        let mut source = split(prime);
        source.push_str("rvK3PencilWeierstrass[P]; rvK3PencilDiscriminant[P];");
        for value in 0..prime {
            writeln!(
                source,
                "rvK3FiberJacobian[P](rvK3PencilParam[P](o,rvUnramElem[F]{{{value}}}));"
            )
            .unwrap();
        }
        source.push_str("rvK3FiberJacobian[P](rvK3PencilParam[P](z,o));");
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        let relative = lines.next().unwrap();
        let discriminant_report = lines.next().unwrap();
        let forms = [
            coefficients(relative, "a2="),
            coefficients(relative, "a4="),
            coefficients(relative, "a6="),
            coefficients(relative, "c4="),
            coefficients(relative, "c6="),
            coefficients(discriminant_report, "coefficients="),
        ];
        assert_eq!(forms.each_ref().map(Vec::len), [5, 9, 13, 9, 13, 25]);
        let infinity = forms.each_ref().map(|form| *form.last().unwrap());
        assert_eq!(coefficients(relative, "discriminant)="), infinity);
        assert_eq!(forms[5], expected_discriminant(prime));
        assert!(discriminant_report.contains("degree=24"));
        assert!(discriminant_report.contains("infinity order=6"));

        let [_, _, _, c4, c6, delta] = &forms;
        let cube = multiply(prime, &multiply(prime, c4, c4), c4);
        let square = multiply(prime, c6, c6);
        let difference: Vec<_> = cube
            .iter()
            .zip(square)
            .map(|(cube, square)| (cube - square).rem_euclid(prime))
            .collect();
        assert_eq!(
            difference,
            delta.iter().map(|c| c * 1728 % prime).collect::<Vec<_>>()
        );
        for (index, report) in lines.enumerate() {
            let actual = ["a2=", "a4=", "a6=", "discriminant="].map(|label| {
                let entry = report
                    .split_once(label)
                    .unwrap()
                    .1
                    .split(']')
                    .next()
                    .unwrap();
                entry.strip_prefix('[').unwrap().parse::<i64>().unwrap()
            });
            let expected = [&forms[0], &forms[1], &forms[2], &forms[5]].map(|form| {
                if index as i64 == prime {
                    *form.last().unwrap()
                } else {
                    evaluate(prime, form, index as i64)
                }
            });
            assert_eq!(actual, expected);
        }
    }
}

#[test]
fn tame_local_orders_distinguish_the_star_and_colliding_multiplicative_fibers() {
    let setup = split(5);
    let output = stdout(interpret(&format!(
        r"{setup}
        rvK3PencilDiscriminant[P];
        rvK3FiberReduction[P](rvK3PencilParam[P](o,z));
        rvK3FiberReduction[P](rvK3PencilParam[P](z,o));
        rvK3FiberReduction[P](rvK3PencilParam[P](o,o));
        rvK3FiberReduction[P](rvK3PencilParam[P](o,two));
        rvK3FiberReduction[P](rvK3PencilParam[P](two,rvUnramElem[F]{{4}}));
        rvK3FiberReduction[P](rvK3PencilParam[P](o,three));
        "
    )));
    let mut lines = output.lines();
    let discriminant = lines.next().unwrap();
    assert!(discriminant.contains("degree=24"));
    assert!(discriminant.ends_with(
        "[(1, 6, [[0], [1]]), (1, 2, [[1], [1]]), (1, 4, [[2], [1]]), (1, 4, [[3], [1]]), (1, 2, [[4], [1]])]"
    ));
    let zero = lines.next().unwrap();
    let infinity = lines.next().unwrap();
    assert!(zero.contains("uniformizer=v/u-[0]"));
    assert!(infinity.contains("uniformizer=u/v"));
    for place in [zero, infinity] {
        assert!(place.contains("scale=0"));
        assert!(place.ends_with("geometric Kodaira=I0*"));
    }
    let first = lines.next().unwrap();
    assert!(first.ends_with("j order=-2, geometric Kodaira=I2"));
    let second = lines.next().unwrap();
    assert!(second.ends_with("j order=-4, geometric Kodaira=I4"));
    assert_eq!(lines.next().unwrap(), second);
    assert!(lines.next().unwrap().ends_with("geometric Kodaira=I4"));
    assert!(lines.next().is_none());
}

#[test]
fn characteristic_three_preserves_orders_divisible_by_three() {
    let source = format!(
        r"
        rvunram F[3,1]{{1,0,1}};
        let z=rvUnramElem[F]{{}};
        let o=rvUnramElem[F]{{1}};
        let two=neg(o);
        let three=rvUnramElem[F]{{0,1}};
        rvhyperelliptic C[F]{{z,two,z,z,z,o}};
        {BASIS}
        rvK3PencilDiscriminant[P];
        rvK3FiberReduction[P](rvK3PencilParam[P](o,z));
        rvK3FiberReduction[P](rvK3PencilParam[P](z,o));
        rvK3FiberReduction[P](rvK3PencilParam[P](o,three));
        rvK3FiberReduction[P](rvK3PencilParam[P](o,add(o,three)));
        "
    );
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    let discriminant = lines.next().unwrap();
    assert!(discriminant.contains("degree=24"));
    assert!(discriminant.contains("infinity order=6"));
    assert!(discriminant.contains("(1, 6, [[0, 0], [1, 0]])"));
    for order in [6, 6, 4, 0] {
        let line = lines.next().unwrap();
        assert!(
            line.ends_with(&format!(
                ", {order}], Kodaira=unclassified in characteristic 3"
            )),
            "{line}"
        );
        assert!(!line.contains("minimal"));
    }
    assert!(lines.next().is_none());
}

#[test]
fn example_and_nominal_type_errors() {
    let output = stdout(example("hyperelliptic-k3-pencil-relative-jacobian"));
    assert!(output.contains("rvK3PencilWeierstrass[P]"));
    assert!(output.contains("geometric Kodaira=I0\n"));
    assert!(output.contains("geometric Kodaira=I2\n"));
    assert!(output.contains("geometric Kodaira=I0*\n"));
    let setup = split(13);
    for (input, message) in [
        ("rvK3PencilWeierstrass[NoSuchPencil];", "unknown K3 pencil"),
        (
            "rvK3FiberReduction[P](z);",
            "expected type RVK3PencilBase(P)",
        ),
        (
            "rvk3pencil P2[B,1,0]; rvK3FiberReduction[P](rvK3PencilParam[P2](o,z));",
            "expected type RVK3PencilBase(P)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {input}")));
        assert!(error.contains(message), "{error}");
    }
}

fn coefficients(report: &str, label: &str) -> Vec<i64> {
    report
        .split_once(label)
        .unwrap()
        .1
        .split_once("]]")
        .unwrap()
        .0
        .split(|c: char| !c.is_ascii_digit())
        .filter(|token| !token.is_empty())
        .map(|token| token.parse().unwrap())
        .collect()
}

fn evaluate(prime: i64, polynomial: &[i64], point: i64) -> i64 {
    polynomial
        .iter()
        .rev()
        .fold(0, |acc, next| (acc * point + next) % prime)
}

fn multiply(prime: i64, a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut out = vec![0; a.len() + b.len() - 1];
    for (i, a) in a.iter().enumerate() {
        for (j, b) in b.iter().enumerate() {
            out[i + j] = (out[i + j] + a * b).rem_euclid(prime);
        }
    }
    out
}

fn expected_discriminant(prime: i64) -> Vec<i64> {
    // For this secant the four branch factors have the form
    // −alpha X + v + alpha²u. Take all six squared resultants.
    let mut polynomial = vec![0; 7];
    polynomial[6] = 16i64.pow(7) % prime;
    for alpha in 1..=4 {
        for beta in alpha + 1..=4 {
            let difference = beta - alpha;
            let factor = [-alpha * beta * difference, difference];
            polynomial = multiply(prime, &polynomial, &multiply(prime, &factor, &factor));
        }
    }
    polynomial.resize(25, 0);
    polynomial
}
