mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

type Weight = (i32, i32, i32);
type Fourier = BTreeMap<Weight, i128>;

fn coefficients(name: &str, weights: impl IntoIterator<Item = Weight>) -> Fourier {
    let weights = weights.into_iter().collect::<Vec<_>>();
    let mut source = String::new();
    for &(a, b, r) in &weights {
        writeln!(source, "{name}[{a}, {b}, {r}];").unwrap();
    }
    let output = stdout(interpret(&source));
    assert_eq!(output.lines().count(), weights.len());
    weights
        .into_iter()
        .zip(output.lines())
        .map(|((a, b, r), line)| {
            let prefix = format!("{name}[{a}, {b}, {r}] = ");
            (
                (a, b, r),
                line.strip_prefix(&prefix).unwrap().parse().unwrap(),
            )
        })
        .collect()
}

#[test]
fn weights_add_and_the_three_simple_reflections_are_isometries() {
    let source = r"
        let e0=rvSiegelWeight[0,0,-2];
        let e1=rvSiegelWeight[2,0,2];
        let e2=rvSiegelWeight[0,2,2];
        let rho=rvSiegelWeight[1,1,1];
        assert add(rho,neg(rho)) == zero[RVSiegelWeight];
        assert add(add(e0,e1),e2) == add(rho,rho);
        assert rvSiegelReflect[0](e0) == neg(e0);
        assert rvSiegelReflect[1](e1) == neg(e1);
        assert rvSiegelReflect[2](e2) == neg(e2);
        assert rvSiegelReflect[0](rho) == rvSiegelWeight[1,1,-1];
        assert rvSiegelReflect[1](rho) == rvSiegelWeight[3,1,3];
        assert rvSiegelReflect[2](rho) == rvSiegelWeight[1,3,3];
        assert rvSiegelReflect[0](rvSiegelReflect[0](rho)) == rho;
        assert rvSiegelReflect[1](rvSiegelReflect[1](rho)) == rho;
        assert rvSiegelReflect[2](rvSiegelReflect[2](rho)) == rho;
        rvSiegelPair(e0,e0); rvSiegelPair(e1,e1); rvSiegelPair(e2,e2);
        rvSiegelPair(e0,e1); rvSiegelPair(e0,e2); rvSiegelPair(e1,e2);
        rvSiegelPair(rho,e0); rvSiegelPair(rho,e1); rvSiegelPair(rho,e2);
        rvSiegelPair(rho,rho);
        rvSiegelPair(rvSiegelReflect[0](rho),rvSiegelReflect[0](rho));
        rvSiegelPair(rvSiegelReflect[1](rho),rvSiegelReflect[1](rho));
        rvSiegelPair(rvSiegelReflect[2](rho),rvSiegelReflect[2](rho));
    ";
    let expected = [2, 2, 2, -2, -2, -2, -1, -1, -1]
        .map(|n| format!("rvSiegelPair = {n}\n"))
        .concat()
        + &"rvSiegelPair = -3/2\n".repeat(4);
    assert_eq!(stdout(interpret(source)), expected);
}

#[test]
fn roots_have_signed_multiplicities_and_reduce_by_discriminant() {
    let source = r"
        rvSiegelRoot(rvSiegelWeight[0,0,-2]);
        rvSiegelRoot(rvSiegelWeight[2,4,6]);
        rvSiegelRoot(neg(rvSiegelWeight[2,4,6]));
        rvSiegelRoot(rvSiegelWeight[2,8,8]);
        rvSiegelRoot(rvSiegelWeight[2,2,2]);
        rvSiegelRoot(rvSiegelWeight[2,2,0]);
        rvSiegelRoot(neg(rvSiegelWeight[2,2,2]));
        rvSiegelRoot(rvSiegelWeight[0,0,-4]);
        rvSiegelRoot(rvSiegelWeight[2,-2,0]);
        rvSiegelRoot(zero[RVSiegelWeight]);
    ";
    assert_eq!(stdout(interpret(source)), concat!(
        "rvSiegelRoot = real, orientation=positive, simple=(1, 0, 0), norm=2, supermultiplicity=1\n",
        "rvSiegelRoot = real, orientation=positive, simple=(0, 1, 2), norm=2, supermultiplicity=1\n",
        "rvSiegelRoot = real, orientation=negative, simple=(0, -1, -2), norm=2, supermultiplicity=1\n",
        "rvSiegelRoot = null, orientation=positive, simple=(1, 1, 4), norm=0, supermultiplicity=10\n",
        "rvSiegelRoot = imaginary, orientation=positive, simple=(1, 1, 1), norm=-6, supermultiplicity=-64\n",
        "rvSiegelRoot = imaginary, orientation=positive, simple=(2, 1, 1), norm=-8, supermultiplicity=108\n",
        "rvSiegelRoot = imaginary, orientation=negative, simple=(-1, -1, -1), norm=-6, supermultiplicity=-64\n",
        "rvSiegelRoot = absent, orientation=positive, simple=(2, 0, 0), norm=8, supermultiplicity=0\n",
        "rvSiegelRoot = absent, orientation=neither, simple=(0, 1, -1), norm=8, supermultiplicity=0\n",
        "rvSiegelRoot = absent, orientation=neither, simple=(0, 0, 0), norm=0, supermultiplicity=0\n",
    ));

    // Unimodular substitutions of uv, u², and u²+v² produce enormous
    // coordinates with the original discriminants.
    let k = 100_000_000_000_000_000_i128;
    let (a, b, r) = (k * (k - 1), k * (k + 1), 2 * k * k - 1);
    let (c, d, s) = (k * k, (k + 1) * (k + 1), 2 * k * (k + 1));
    let (e, f, t) = (k * k + d, (k - 1) * (k - 1) + c, 4 * k * k);
    let output = stdout(interpret(&format!(
        "rvSiegelRoot(rvSiegelWeight[{},{},{}]); rvSiegelRoot(rvSiegelWeight[{},{},{}]);\
         rvSiegelRoot(rvSiegelWeight[{},{},{}]);",
        2 * a,
        2 * b,
        2 * r,
        2 * c,
        2 * d,
        2 * s,
        2 * e,
        2 * f,
        2 * t
    )));
    assert_eq!(output, format!(
        "rvSiegelRoot = real, orientation=positive, simple=(1, {a}, {b}), norm=2, supermultiplicity=1\n\
         rvSiegelRoot = null, orientation=positive, simple=(1, {c}, {d}), norm=0, supermultiplicity=10\n\
         rvSiegelRoot = imaginary, orientation=positive, simple=(2, {e}, {f}), norm=-8, supermultiplicity=108\n"
    ));
}

#[test]
fn the_denominator_is_anti_invariant_under_actual_reflections() {
    let probes = [
        (1, 1, 1),
        (1, 3, 1),
        (3, 1, 3),
        (3, 3, 1),
        (3, 3, -3),
        (3, 5, 1),
        (1, 1, 3),
    ];
    let images = probes.map(|(a, b, r)| {
        [
            (a, b, -r),
            (a + 4 * b - 2 * r, b, 4 * b - r),
            (a, b + 4 * a - 2 * r, 4 * a - r),
        ]
    });
    let values = coefficients(
        "rvSiegelDelta5",
        probes.into_iter().chain(images.into_iter().flatten()),
    );
    let mut source = String::new();
    let mut expected = String::new();
    for ((a, b, r), reflected) in probes.into_iter().zip(images) {
        let value = values[&(a, b, r)];
        for (i, (x, y, z)) in reflected.into_iter().enumerate() {
            assert_eq!(values[&(x, y, z)], -value);
            writeln!(
                source,
                "rvSiegelDelta5At(rvSiegelReflect[{i}](rvSiegelWeight[{a},{b},{r}]));"
            )
            .unwrap();
            writeln!(expected, "rvSiegelDelta5[{x}, {y}, {z}] = {}", -value).unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), expected);
    assert_eq!(values[&(1, 1, 1)], 1);
    assert_eq!(values[&(1, 3, 1)], -9);
    assert_eq!(values[&(3, 3, 1)], -90);
}

#[test]
fn squaring_all_low_degree_coefficients_recovers_phi10() {
    // Each contributing root factor has |r| <= a+b. The Weyl prefactor
    // changes the doubled Laurent bound to |R| <= A+B-1 for Δ₅.
    let delta = coefficients(
        "rvSiegelDelta5",
        (1..=5).step_by(2).flat_map(|a| {
            (1..=5)
                .step_by(2)
                .flat_map(move |b| (1 - a - b..=a + b - 1).step_by(2).map(move |r| (a, b, r)))
        }),
    );
    let phi = coefficients(
        "rvIgusaPhi10",
        (0..=3).flat_map(|a| (0..=3).flat_map(move |b| (-a - b..=a + b).map(move |r| (a, b, r)))),
    );
    for ((a, b, r), value) in phi {
        let square: i128 = delta
            .iter()
            .map(|(&(x, y, z), coefficient)| {
                coefficient * delta.get(&(2 * a - x, 2 * b - y, 2 * r - z)).unwrap_or(&0)
            })
            .sum();
        assert_eq!(square, value, "P^{a} q^{b} y^{r}");
    }
    let boundaries = coefficients(
        "rvSiegelDelta5",
        [
            (0, 1, 1),
            (1, 0, -1),
            (-1, 1, 1),
            (2, 1, 1),
            (1, 2, 1),
            (1, 1, 0),
            (1, 1, -1),
            (3, 1, 1),
            (1, 3, 3),
        ],
    );
    for weight in [
        (0, 1, 1),
        (1, 0, -1),
        (-1, 1, 1),
        (2, 1, 1),
        (1, 2, 1),
        (1, 1, 0),
    ] {
        assert_eq!(boundaries[&weight], 0);
    }
    assert_eq!(boundaries[&(1, 1, -1)], -1);
    assert_eq!(boundaries[&(3, 1, 1)], -9);
    assert_eq!(boundaries[&(1, 3, 3)], -1);
}

#[test]
fn examples_and_diagnostics_cover_exact_ranges_and_nominal_types() {
    assert_eq!(stdout(example("siegel-borcherds-denominator")), concat!(
        "rvSiegelPair = -1\n",
        "rvSiegelRoot = real, orientation=positive, simple=(0, 1, 0), norm=2, supermultiplicity=1\n",
        "rvSiegelRoot = imaginary, orientation=positive, simple=(1, 1, 1), norm=-6, supermultiplicity=-64\n",
        "rvSiegelDelta5[1, 1, 1] = 1\n",
        "rvSiegelDelta5[3, 1, 3] = -1\n",
        "rvIgusaPhi10[1, 1, 1] = 1\n",
    ));
    for (source, expected) in [
        (
            "rvSiegelRoot(rvSiegelWeight[1,1,1]);",
            "integral weight coordinates",
        ),
        (
            "show rvSiegelReflect[3](zero[RVSiegelWeight]);",
            "index must be between 0 and 2",
        ),
        (
            "rvunram F[3,1]{0,1}; rvSiegelRoot(rvUnramElem[F]{1});",
            "an RVSiegelWeight",
        ),
        (
            "rvSiegelDelta5[2000001,2000001,1];",
            "coefficient storage limit",
        ),
        (
            "rvSiegelRoot(rvSiegelWeight[2000000,2000000,0]);",
            "coefficient storage limit",
        ),
        (
            "show neg(rvSiegelWeight[-170141183460469231731687303715884105728,0,0]);",
            "doubled coordinate exceeds signed 128-bit range",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
    // 4B overflows a signed coordinate, but 4B−R and the pairing fit.
    let r = 1i128 << 126;
    let b = (r / 2) + 1;
    assert_eq!(
        stdout(interpret(&format!(
            "assert rvSiegelReflect[1](rvSiegelWeight[0,{b},{r}]) == rvSiegelWeight[4,{b},{}];",
            r + 4
        ))),
        ""
    );
}
