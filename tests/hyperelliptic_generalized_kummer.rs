mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,o,z,z,z,o};
";

#[test]
fn compares_resolved_k3_counts_with_ordinary_and_certified_level_fibers() {
    assert_eq!(
        stdout(example("hyperelliptic-generalized-kummer")),
        concat!(
            "rvJacHilbPoints[C, 0, 1] = 1\n",
            "rvJacHilbPoints[C, 3, 1] = 1920\n",
            "rvGenKummerPoints[C, 2] = 24, smooth generalized Kummer fiber\n",
            "rvKummerK3Points[C, 1] = singular=12, nodes=4, resolved=24\n",
            "rvGenKummerPoints[C, 3] = 178, summation fiber\n",
            "rvGenKummerFiber[C, 3] = 151, summation fiber\n",
            "rvGenKummerPoints[C, 4] = 2220, smooth generalized Kummer fiber\n",
            "rvGenKummerFiber[C, 4] = 1680, smooth generalized Kummer fiber\n",
            "rvGenKummerLevel[B, C, 1] = 24, smooth generalized Kummer fiber, length=2, rational torsion order=4\n",
            "rvGenKummerLevel[B, C, 2] = 248, smooth generalized Kummer fiber, length=2, rational torsion order=16\n",
            "rvGenKummerPoints[L, 2] = 248, smooth generalized Kummer fiber\n",
            "rvKummerK3Points[C, 2] = singular=104, nodes=16, resolved=248\n",
        )
    );
}

#[test]
fn all_fibers_agree_with_an_independent_closed_point_product() {
    let mut source = format!(
        r"{F3}
        let p=rvJacPoint[C](rvUnramElem[F]{{2}},o);
        let t=rvJacPoint[C](z,z);
        assert rvJacScale[6](p)==zero[RVJac(C)];
        assert rvJacU[2](rvJacScale[3](p))==o;
        assert rvJacScale[2](t)==zero[RVJac(C)];
        assert rvJacU[1](t)==o;
    "
    );
    for i in 0..12 {
        writeln!(
            source,
            "let d{i}=add(rvJacScale[{}](p),rvJacScale[{}](t));",
            i % 6,
            i / 6
        )
        .unwrap();
    }
    for n in 1..=4 {
        writeln!(source, "rvJacHilbPoints[C,{n},1];").unwrap();
        for i in 0..12 {
            writeln!(source, "rvGenKummerFiber[C,{n}](d{i});").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let reported: Vec<_> = output.lines().map(count).collect();
    let direct = literal_cycle_product();
    for n in 1..=4 {
        let group = &reported[13 * (n - 1)..13 * n];
        assert_eq!(group[0], group[1..].iter().sum());
        assert_eq!(group[1..], direct[n]);
        for d in 0..12 {
            for x in 0..12 {
                assert_eq!(direct[n][d], direct[n][add(d, times(n, x))]);
            }
        }
    }
    assert_eq!(direct[1], [1; 12]);
    assert_ne!(direct[2][0], direct[2][1]);
    assert_ne!(direct[3][0], direct[3][1]);
    assert_ne!(direct[4][0], direct[4][1]);

    // Unordered rational pairs are one part of Hilb². Quadratic closed
    // orbits contribute the other reduced cycles, and a double rational
    // point contributes a projective line of length-two subschemes.
    let pairs = (0..12).flat_map(|a| (a + 1..12).map(move |b| add(a, b)));
    let orbit_two = closed_orbits()[2];
    let mut at_two = orbit_two;
    for sum in pairs {
        at_two[sum] += 1;
    }
    for d in 0..12 {
        at_two[times(2, d)] += 4;
    }
    assert_eq!(direct[2], at_two);
}

fn add(a: usize, b: usize) -> usize {
    ((a % 6 + b % 6) % 6) + 6 * ((a / 6 + b / 6) % 2)
}

fn times(n: usize, a: usize) -> usize {
    (0..n).fold(0, |sum, _| add(sum, a))
}

fn closed_orbits() -> [[u128; 12]; 5] {
    // For t^4+2t^2+9 the four point counts follow from the roots of
    // y^2+2y+9 by taking y=x^2. The trace has fibers of equal size.
    let points = [0, 12, 144, 684, 9216];
    let mut exact = [[0; 12]; 5];
    for n in 1..=4 {
        let mut by_trace = [points[n] / 12; 12];
        for d in 1..n {
            if n % d == 0 {
                for trace in 0..12 {
                    by_trace[times(n / d, trace)] -= exact[d][trace] * d as u128;
                }
            }
        }
        for (trace, count) in by_trace.into_iter().enumerate() {
            assert_eq!(count % n as u128, 0);
            exact[n][trace] = count / n as u128;
        }
    }
    exact
}

fn literal_cycle_product() -> [[u128; 12]; 5] {
    let mut choices = [[0; 12]; 5];
    choices[0][0] = 1;
    for (degree, by_trace) in closed_orbits().into_iter().enumerate().skip(1) {
        let q = 3u128.pow(degree as u32);
        let punctual = [1, 1, 1 + q, 1 + q + q * q, 1 + q + 2 * q * q + q * q * q];
        for (trace, orbits) in by_trace.into_iter().enumerate() {
            for _ in 0..orbits {
                for length in (1..=4).rev() {
                    for multiplicity in 1..=length / degree {
                        let shift = times(multiplicity, trace);
                        for d in 0..12 {
                            choices[length][add(d, shift)] +=
                                choices[length - multiplicity * degree][d] * punctual[multiplicity];
                        }
                    }
                }
            }
        }
    }
    choices
}

fn count(line: &str) -> u128 {
    line.split_once(" = ")
        .unwrap()
        .1
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

#[test]
fn full_odd_prime_level_uses_frobenius_without_a_group_scan() {
    let source = r"
        rvunram W[5,1]{3,0,1};
        let z=rvUnramElem[W]{};
        let o=rvUnramElem[W]{1};
        let a=rvUnramElem[W]{0,1};
        rvhyperelliptic C[W]{z,o,z,z,z,o};
        let e1=rvJacScale[2](rvJacPoint[C](o,a));
        let e2=add(e1,rvJacScale[2](rvJacPoint[C](rvUnramElem[W]{1,2},a)));
        let f1=rvJacScale[4](rvJacPoint[C](rvUnramElem[W]{1,1},a));
        let f2=rvJacScale[4](rvJacPoint[C](rvUnramElem[W]{1,3},a));
        rvjactorsionbasis B[C,3](rvUnramElem[W]{2,2}){e1,e2,f1,f2};
        rvGenKummerLevel[B,C,1];
    ";
    // The scalar Frobenius is -5, so the 3-torsion is all fixed.
    let [n1, n2, n3] = [6u128.pow(4), 24u128.pow(4), 126u128.pow(4)];
    let h3 = (n1.pow(3) + 3 * n1 * n2 + 2 * n3) / 6 + 25 * n1 * n1 + 25 * 25 * n1;
    let expected = h3 / n1 + 80 * 25 * 25;
    assert_eq!(
        stdout(interpret(source)),
        format!("rvGenKummerLevel[B, C, 1] = {expected}, smooth generalized Kummer fiber, length=3, rational torsion order=81\n")
    );
}

#[test]
fn coprime_length_five_has_equal_fibers() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            rvJacHilbPoints[C,5,1];
            rvGenKummerPoints[C,5];
            rvGenKummerFiber[C,5](rvJacPoint[C](rvUnramElem[F]{{2}},o));
        "
        ))),
        concat!(
            "rvJacHilbPoints[C, 5, 1] = 215832\n",
            "rvGenKummerPoints[C, 5] = 17986, smooth generalized Kummer fiber\n",
            "rvGenKummerFiber[C, 5] = 17986, smooth generalized Kummer fiber\n",
        )
    );
}

#[test]
fn diagnoses_wrong_curves_zero_arguments_and_unrepresentable_counts() {
    for (statement, expected) in [
        (
            "rvJacHilbPoints[C,0,0];",
            "extension degree must be positive",
        ),
        (
            "rvGenKummerPoints[C,0];",
            "summation fiber length must be positive",
        ),
        (
            "rvGenKummerFiber[C,0](zero[RVJac(C)]);",
            "summation fiber length must be positive",
        ),
        (
            "rvJacHilbPoints[C,1000,1];",
            "checked 128-bit arithmetic range",
        ),
        (
            "rvhyperelliptic E[F]{z,o,z,o}; rvJacHilbPoints[E,1,1];",
            "require a genus-two curve",
        ),
        (
            "rvhyperelliptic D[F]{z,o,z,z,z,o}; rvGenKummerFiber[C,1](zero[RVJac(D)]);",
            "expected type RVJac(C)",
        ),
        (
            "rvGenKummerLevel[Unknown,C,1];",
            "unknown Jacobian torsion basis",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F3}{statement}")));
        assert!(error.contains(expected), "{error}");
    }
    assert_eq!(
        stdout(interpret(&format!(
            "{F3}rvJacHilbPoints[C,0,340282366920938463463374607431768211455];"
        ))),
        "rvJacHilbPoints[C, 0, 340282366920938463463374607431768211455] = 1\n",
    );
    let example = include_str!("../examples/hyperelliptic-generalized-kummer.asl");
    let error = diagnostic(interpret(&format!(
        "{example}\nrvhyperelliptic E[F]{{one,z,z,z,z,one}}; rvGenKummerLevel[B,E,1];"
    )));
    assert!(error.contains("not the declared base extension"), "{error}");
}
