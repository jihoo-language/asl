mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup(prime: i64, extension: bool) -> String {
    let mut source = format!(
        r"
        rvunram K[{prime},1]{{0,1}};
        let zk=rvUnramElem[K]{{}};
        let ok=rvUnramElem[K]{{1}};
        let c1=rvUnramElem[K]{{24}};
        let c2=rvUnramElem[K]{{-50}};
        let c3=rvUnramElem[K]{{35}};
        let c4=rvUnramElem[K]{{-10}};
        rvhyperelliptic C[K]{{zk,c1,c2,c3,c4,ok}};
        "
    );
    let (curve, field) = if extension {
        source.push_str("rvunram F[5,1]{2,0,1}; rvjacbase L[C,F];");
        ("L", "F")
    } else {
        ("C", "K")
    };
    write!(
        source,
        r"
        let z=rvUnramElem[{field}]{{}};
        let o=rvUnramElem[{field}]{{1}};
        let two=rvUnramElem[{field}]{{2}};
        let three=rvUnramElem[{field}]{{3}};
        let r0=rvJacPoint[{curve}](z,z);
        let r1=rvJacPoint[{curve}](o,z);
        let r2=rvJacPoint[{curve}](two,z);
        let r3=rvJacPoint[{curve}](three,z);
        rvjactorsionbasis B[{curve},2](neg(o)){{r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)}};
        rvk3pencil P[B,0,1];
        "
    )
    .unwrap();
    source
}

#[test]
fn independent_euler_product_matches_over_prime_and_quadratic_fields() {
    let base = euler_product(5, false, 8);
    let cover = euler_product(5, true, 8);
    let mut even_product = vec![0; 2 * base.len() - 1];
    for (i, &a) in base.iter().enumerate() {
        for (j, &b) in base.iter().enumerate() {
            even_product[i + j] += if j % 2 == 0 { a * b } else { -a * b };
        }
    }
    for (j, coefficient) in even_product.into_iter().enumerate() {
        let expected = if j % 2 == 0 { cover[j / 2] } else { 0 };
        assert_eq!(coefficient, expected);
    }
    euler_product(13, false, 10);
}

fn euler_product(prime: i64, extension: bool, conductor: usize) -> Vec<i128> {
    let mut source = setup(prime, extension);
    source.push_str(
        "rvK3PencilLFunction[P];
         let inf=rvK3PencilParam[P](z,o);
         rvK3FiberEuler[P](inf);
         rvK3FiberJacobian[P](inf);
         rvK3PencilFiber[P](inf);",
    );
    let field = if extension { "F" } else { "K" };
    let elements = elements(prime, extension);
    for (index, [constant, linear]) in elements.iter().enumerate() {
        let input = if extension {
            format!("{constant},{linear}")
        } else {
            constant.to_string()
        };
        writeln!(
            source,
            "let p{index}=rvK3PencilParam[P](o,rvUnramElem[{field}]{{{input}}});
             rvK3FiberEuler[P](p{index});
             rvK3FiberJacobian[P](p{index});
             rvK3PencilFiber[P](p{index});"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    let report = lines.next().unwrap();
    let polynomial = list(report, "coefficients=[");
    assert_eq!(polynomial.len() - 1, conductor - 4);
    assert_eq!(value(report, "conductor degree="), conductor as i128);
    let q = elements.len() as i128;
    let (mut sum, mut trivial, mut seen, mut observed_conductor) = (0, 2, 0, 0);
    while let Some(euler) = lines.next() {
        seen += 1;
        let factor = list(euler, "factor=[");
        let cycles = list(euler, "component cycles=[");
        let actual_trace = -factor.get(1).copied().unwrap_or(0);
        let jacobian = lines.next().unwrap();
        let coefficients = ["a2=[", "a4=[", "a6=["].map(|label| {
            let coefficients = list(jacobian, label);
            [
                coefficients[0] as i64,
                coefficients.get(1).copied().unwrap_or(0) as i64,
            ]
        });
        let trace = trace_cubic(prime, &elements, coefficients);
        assert_eq!(trace, actual_trace, "{euler}");
        let points = lines
            .next()
            .unwrap()
            .split_once(" = ")
            .unwrap()
            .1
            .parse::<i128>()
            .unwrap();
        let fixed = cycles.iter().filter(|&&cycle| cycle == 1).count() as i128;
        assert_eq!(points, 1 - trace + q * (1 + fixed), "{euler}");
        sum += trace;
        trivial += cycles.len() as i128;
        observed_conductor += value(euler, "conductor order=");
    }
    assert_eq!(seen, q + 1);
    assert_eq!(observed_conductor, conductor as i128);
    assert_eq!(sum, polynomial[1]);
    let rank = value(report, "Mordell-Weil rank=");
    assert_eq!(value(report, "NS rank="), rank + trivial);
    assert_eq!(value(report, "trivial rank="), trivial);
    polynomial
}

fn value(report: &str, label: &str) -> i128 {
    report
        .split_once(label)
        .unwrap()
        .1
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

fn list(report: &str, label: &str) -> Vec<i128> {
    let text = report
        .split_once(label)
        .unwrap()
        .1
        .split(']')
        .next()
        .unwrap();
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(',')
            .map(|part| part.trim().parse().unwrap())
            .collect()
    }
}

fn elements(prime: i64, extension: bool) -> Vec<[i64; 2]> {
    let mut values = Vec::new();
    for constant in 0..prime {
        for linear in 0..if extension { prime } else { 1 } {
            values.push([constant, linear]);
        }
    }
    values
}

fn trace_cubic(prime: i64, elements: &[[i64; 2]], [a2, a4, a6]: [[i64; 2]; 3]) -> i128 {
    let mult = |[a, b]: [i64; 2], [c, d]: [i64; 2]| {
        [
            (a * c - 2 * b * d).rem_euclid(prime),
            (a * d + b * c).rem_euclid(prime),
        ]
    };
    let add =
        |[a, b]: [i64; 2], [c, d]: [i64; 2]| [(a + c).rem_euclid(prime), (b + d).rem_euclid(prime)];
    let squares: Vec<_> = elements.iter().map(|&y| mult(y, y)).collect();
    let affine: usize = elements
        .iter()
        .map(|&x| {
            let right = add(mult(add(mult(add(x, a2), x), a4), x), a6);
            squares.iter().filter(|&&square| square == right).count()
        })
        .sum();
    elements.len() as i128 - affine as i128
}

#[test]
fn split_and_nonsplit_multiplicative_components_have_different_actions() {
    let setup = setup(5, false);
    let output = stdout(interpret(&format!(
        r"{setup}
        rvK3FiberEuler[P](rvK3PencilParam[P](o,o));
        rvK3FiberEuler[P](rvK3PencilParam[P](o,two));
        rvK3FiberEuler[P](rvK3PencilParam[P](z,o));
        rvK3PencilLFunction[P];
        "
    )));
    assert_eq!(
        output,
        concat!(
            "rvK3FiberEuler[P] = q=5, geometric Kodaira=I2, factor=[1, -1], conductor order=1, split=yes, component cycles=[1], Tamagawa=2\n",
            "rvK3FiberEuler[P] = q=5, geometric Kodaira=I4, factor=[1, 1], conductor order=1, split=no, component cycles=[1, 2], Tamagawa=2\n",
            "rvK3FiberEuler[P] = q=5, geometric Kodaira=I0*, factor=[1], conductor order=2, component cycles=[1, 1, 1, 1], Tamagawa=4\n",
            "rvK3PencilLFunction[P] = q=5, coefficients=[1, 0, -50, 0, 625], degree=4, conductor degree=8, NS rank=18, trivial rank=16, Mordell-Weil rank=2\n",
        )
    );
}

#[test]
fn example_and_scope_diagnostics() {
    let output = stdout(example("hyperelliptic-k3-pencil-l-function"));
    assert!(output.contains("factor=[1, -2, 13]"));
    assert!(output.contains("coefficients=[1, 10, -169, -3380, -28561, 285610, 4826809]"));
    assert!(
        output.ends_with("conductor degree=10, NS rank=18, trivial rank=16, Mordell-Weil rank=2\n")
    );
    let p5 = setup(5, false);
    let error = diagnostic(interpret(&format!("{p5} rvK3FiberEuler[P](o);")));
    assert!(error.contains("expected type RVK3PencilBase(P)"), "{error}");

    let p3 = r"
        rvunram F[3,1]{1,0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let i=rvUnramElem[F]{0,1};
        rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
        let r0=rvJacPoint[C](z,z);
        let r1=rvJacPoint[C](o,z);
        let r2=rvJacPoint[C](neg(o),z);
        let r3=rvJacPoint[C](i,z);
        rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
        rvk3pencil P[B,0,1];
    ";
    for command in [
        "rvK3PencilLFunction[P];",
        "rvK3FiberEuler[P](rvK3PencilParam[P](o,z));",
    ] {
        let error = diagnostic(interpret(&format!("{p3}{command}")));
        assert!(
            error.contains("require characteristic at least five"),
            "{error}"
        );
    }
}
