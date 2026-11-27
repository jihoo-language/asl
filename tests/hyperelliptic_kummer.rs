mod common;

use std::collections::BTreeSet;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
";

#[test]
fn common_branch_factors_and_symmetric_difference() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer")),
        concat!(
            "rvJacKummer = B=[[1], [0], [0], [0], [1]], nonsquare orbits by degree=[(2, 2)]\n",
            "rvJacKummer = B=[[0], [2], [2], [1]], nonsquare orbits by degree=[(1, 1), (2, 1)]\n",
            "rvJacKummer = B=[[0], [2], [1], [1]], nonsquare orbits by degree=[(1, 1), (2, 1)]\n",
            "rvJacKummer = B=[[1]], nonsquare orbits by degree=[]\n",
            "rvJacKummer = B=[[1], [0], [0], [1], [2], [1]], nonsquare orbits by degree=[(2, 1), (3, 1)]\n",
        )
    );
}

#[test]
fn extension_checks_square_classes_in_each_new_branch_field() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}
            rvunram E[3,1]{{1,0,1}};
            rvunram O[3,1]{{1,2,0,1}};
            rvhyperelliptic C[F]{{z,o,o,o,t,o}};
            rvjacbase CE[C,E];
            rvjacbase CO[C,O];
            let ramified = rvJacPoint[C](t,z);
            assert rvJacDivisible2(ramified) == S.lo;
            assert rvJacDivisible2(rvJacExtend[CE](ramified)) == S.hi;
            assert rvJacDivisible2(rvJacExtend[CO](ramified)) == S.lo;

            rvhyperelliptic M[F]{{o,z,z,o,t,o}};
            rvjacbase ME[M,E];
            rvjacbase MO[M,O];
            let point = rvJacPoint[M](t,o);
            let even = rvJacExtend[ME](point);
            let odd = rvJacExtend[MO](point);
            assert rvJacDivisible2(even) == S.lo;
            assert rvJacDivisible2(odd) == S.lo;
            rvJacKummer(even);
            rvJacKummer(odd);
            assert rvJacKummerSupport[2](even) == rvUnramElem[E]{{1}};
            assert rvJacKummerSupport[3](even) == rvUnramElem[E]{{}};
            assert rvJacDivisible2(add(odd,odd)) == S.hi;
        "
        ))),
        concat!(
            "rvJacKummer = B=[[1, 0], [0, 0], [1, 0]], nonsquare orbits by degree=[(1, 2)]\n",
            "rvJacKummer = B=[[1, 0, 0], [0, 0, 0], [0, 0, 0], [1, 0, 0], [2, 0, 0], [1, 0, 0]], nonsquare orbits by degree=[(1, 3), (2, 1)]\n",
        )
    );
}

#[test]
fn support_coefficients_can_lie_outside_the_prime_field() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram K[3,1]{1,0,1};
            let z=rvUnramElem[K]{};
            let o=rvUnramElem[K]{1};
            let b=rvUnramElem[K]{1,1};
            rvhyperelliptic C[K]{z,b,neg(add(o,b)),o};
            let d=rvJacPoint[C](b,z);
            assert rvJacKummerSupport[1](d) == neg(b);
            assert rvJacKummerSupport[2](d) == o;
            assert rvJacDivisible2(d) == S.lo;
            rvJacKummer(d);
        "
        )),
        "rvJacKummer = B=[[0, 0], [2, 2], [1, 0]], nonsquare orbits by degree=[(1, 2)]\n"
    );
}

#[test]
fn finite_field_cardinalities_beyond_u64() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram E[65537,1]{4,1,0,1};
            let z=rvUnramElem[E]{};
            let o=rvUnramElem[E]{1};
            let five=rvUnramElem[E]{5};
            rvhyperelliptic C[E]{z,five,o,z,z,o};
            let d=rvJacPoint[C](z,z);
            // The quartic branch factor stays irreducible in the cubic extension.
            // Its residue cardinality exceeds u64, and its constant 5 is a nonsquare.
            assert rvJacDivisible2(d) == S.lo;
            assert rvJacKummerSupport[1](d) == five;
            assert rvJacKummerSupport[5](d) == o;

            rvunram F[65537,1]{3,1,0,0,0,1};
            let zero=rvUnramElem[F]{};
            let one=rvUnramElem[F]{1};
            let three=rvUnramElem[F]{3};
            let minus_four=rvUnramElem[F]{-4};
            rvhyperelliptic D[F]{zero,three,minus_four,one};
            let e=rvJacPoint[D](zero,zero);
            // The base field itself exceeds u64. The support is X(X-3).
            assert rvJacDivisible2(e) == S.lo;
            assert rvJacKummerSupport[0](e) == zero;
            assert rvJacKummerSupport[1](e) == neg(three);
            assert rvJacKummerSupport[2](e) == one;
            assert rvJacKummerSupport[3](e) == zero;
        "
        )),
        ""
    );
}

#[test]
fn exhaustive_comparison_with_doubling_on_a_genus_two_group() {
    let mut source = format!("{F3} rvhyperelliptic C[F]{{z,o,z,z,z,o}};");
    let mut count = 0;
    for degree in 0..=2 {
        for i in 0..3_i32.pow(degree as u32) {
            let mut u = ternary_digits(i, degree);
            u.push(1);
            for j in 0..3_i32.pow(degree as u32) {
                let v = ternary_digits(j, degree);
                if !valid_pair(&u, &v) {
                    continue;
                }
                writeln!(
                    source,
                    "let d{count}=rvjac[C]{{{}}}{{{}}}; show d{count}; show add(d{count},d{count}); show rvJacDivisible2(d{count});",
                    payload(&u),
                    payload(&v),
                )
                .unwrap();
                count += 1;
            }
        }
    }
    assert_eq!(count, 12);
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 3 * count);
    let mut actual_image = BTreeSet::new();
    let mut detected = BTreeSet::new();
    for row in lines.chunks_exact(3) {
        actual_image.insert(row[1]);
        match row[2] {
            "S.hi" => {
                detected.insert(row[0]);
            }
            "S.lo" => {}
            flag => panic!("unexpected divisibility result {flag}"),
        }
    }
    assert_eq!(actual_image.len(), 3);
    assert_eq!(detected, actual_image);
}

fn ternary_digits(mut n: i32, degree: usize) -> Vec<i32> {
    (0..degree)
        .map(|_| {
            let digit = n % 3;
            n /= 3;
            digit
        })
        .collect()
}

fn valid_pair(u: &[i32], v: &[i32]) -> bool {
    let mut remainder = [0, 1, 0, 0, 0, 1];
    for (i, a) in v.iter().enumerate() {
        for (j, b) in v.iter().enumerate() {
            remainder[i + j] -= a * b;
        }
    }
    let degree = u.len() - 1;
    for i in (degree..remainder.len()).rev() {
        let leading = remainder[i];
        for (j, c) in u.iter().enumerate() {
            remainder[i - degree + j] -= leading * c;
        }
    }
    remainder[..degree].iter().all(|&c| c % 3 == 0)
}

fn payload(coefficients: &[i32]) -> String {
    coefficients
        .iter()
        .map(|c| format!("[{c}]"))
        .collect::<Vec<_>>()
        .join(",")
}

#[test]
fn rejects_non_jacobian_inputs_and_mismatched_coefficient_fields() {
    for misuse in [
        "show rvJacDivisible2(o);",
        "show rvJacKummerSupport[0](S.hi);",
        "rvJacKummer(o);",
    ] {
        let error = diagnostic(interpret(&format!("{F3} {misuse}")));
        assert!(error.contains("RVJac divisor class"), "{error}");
    }
    let error = diagnostic(interpret(&format!(
        r"{F3}
        rvunram G[3,1]{{0,1}};
        rvhyperelliptic C[F]{{z,o,z,z,z,o}};
        assert rvJacKummerSupport[0](zero[RVJac(C)]) == rvUnramElem[G]{{1}};
        "
    )));
    assert!(error.contains("expected type RVUnram(F)"), "{error}");
}
