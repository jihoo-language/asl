mod common;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,t,z,t,o,o};
";

const F25: &str = r"
    rvunram W[5,1]{2,0,1};
    let z=rvUnramElem[W]{};
    let o=rvUnramElem[W]{1};
    rvhyperelliptic S[W]{z,o,z,z,z,o};
";

fn coefficients(line: &str) -> Vec<i128> {
    let output = line.split_once(" = ").unwrap().1.trim();
    output[1..output.len() - 1]
        .split(", ")
        .map(|term| term.parse().unwrap())
        .collect()
}

fn product(left: &[i128], right: &[i128]) -> Vec<i128> {
    let mut result = vec![0; left.len() + right.len() - 1];
    for (i, &a) in left.iter().enumerate() {
        for (j, &b) in right.iter().enumerate() {
            result[i + j] += a * b;
        }
    }
    result
}

// The known K3 denominator in H² is the cross quartic, two abelian factors
// 1-3T, four fixed nodes, and four node orbits of length three.
fn surface_h2() -> Vec<i128> {
    let mut p = vec![1, 8, 33, 72, 81];
    for _ in 0..6 {
        p = product(&p, &[1, -3]);
    }
    for _ in 0..4 {
        p = product(&p, &[1, 0, 0, -27]);
    }
    p
}

#[test]
fn runnable_example_and_nonsplit_surface_polynomial() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-hilbert-characteristic")),
        concat!(
            "rvKummerHilbChar[C, 1, 1, 2, 2, 12] = [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1]\n",
            "rvKummerHilbChar[C, 2, 1, 4, 2, 8] = [1, 1, 1, 0, 0, 0, 1, 1, 1]\n",
            "rvKummerHilbChar[C, 2, 1, 4, 5, 8] = [1, 0, 4, 0, 4, 4, 1, 1, 0]\n",
            "rvKummerHilbChar[C, 1, 1, 2, 3, 6] = [1, 2, 0, 0, 0, 0, 0]\n",
        )
    );
    let mut expected = surface_h2();
    expected.resize(27, 0);
    for ell in [2, 3, 5] {
        let output = stdout(interpret(&format!(
            "{F3} rvKummerHilbChar[C,1,1,2,{ell},26];"
        )));
        assert_eq!(
            coefficients(&output),
            expected
                .iter()
                .map(|value| value.rem_euclid(ell))
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn hilbert_square_matches_integral_newton_coefficients_from_symmetric_square() {
    let h2 = surface_h2();
    let order = 8;
    let mut sums = vec![0];
    for j in 1..=2 * order {
        let trace = -(1..=j.min(h2.len() - 1))
            .map(|i| h2[i] * if i == j { j as i128 } else { sums[j - i] })
            .sum::<i128>();
        sums.push(trace);
    }
    let traces: Vec<_> = (0..=order)
        .map(|j| {
            let q = 3i128.pow(j as u32);
            // Sym² H², the surface top class, and the length-two oscillator of H².
            (sums[j] * sums[j] + sums[2 * j]) / 2 + q * q + q * sums[j]
        })
        .collect();
    let mut exact = vec![1];
    for j in 1..=order {
        let weighted = -(1..=j).map(|i| exact[j - i] * traces[i]).sum::<i128>();
        assert_eq!(weighted % j as i128, 0);
        exact.push(weighted / j as i128);
    }
    for ell in [2, 3, 5, 11] {
        let output = stdout(interpret(&format!(
            "{F3} rvKummerHilbChar[C,2,1,4,{ell},{order}];"
        )));
        assert_eq!(
            coefficients(&output),
            exact
                .iter()
                .map(|value| value.rem_euclid(ell))
                .collect::<Vec<_>>()
        );
    }
    let mut expected = product(&h2, &[1, -3]);
    expected.resize(27, 0);
    let output = stdout(interpret(&format!("{F3} rvKummerHilbChar[C,2,1,2,2,26];")));
    assert_eq!(
        coefficients(&output),
        expected
            .iter()
            .map(|value| value.rem_euclid(2))
            .collect::<Vec<_>>()
    );
    let output = stdout(interpret(&format!("{F3} rvKummerHilbChar[C,2,1,6,7,26];")));
    assert_eq!(
        coefficients(&output),
        expected
            .iter()
            .enumerate()
            .map(|(j, value)| (value * 9i128.pow(j as u32)).rem_euclid(7))
            .collect::<Vec<_>>()
    );
}

#[test]
fn scalar_frobenius_has_binomial_coefficients_through_and_beyond_the_top() {
    for (n, d, dimension, order) in [(2, 1, 23, 27), (2, 2, 276, 280), (3, 3, 2554, 9)] {
        for ell in [2, 3, 7] {
            let i = 2 * d;
            let output = stdout(interpret(&format!(
                "{F25} rvKummerHilbChar[S,{n},1,{i},{ell},{order}];"
            )));
            let factor = 25i128.pow(d);
            let mut expected = vec![0; order + 1];
            expected[0] = 1;
            for _ in 0..dimension {
                for j in (1..=order).rev() {
                    expected[j] = (expected[j] - factor * expected[j - 1]).rem_euclid(ell);
                }
            }
            assert_eq!(coefficients(&output), expected);
        }
    }
}

#[test]
fn huge_extension_powers_do_not_use_exact_traces_or_multiply_exponents() {
    let huge = u128::MAX;
    let output = stdout(interpret(&format!(
        "{F3} rvKummerHilbChar[C,2,{huge},4,2,12]; rvKummerHilbChar[C,2,3,4,2,12];"
    )));
    let values: Vec<_> = output.lines().map(coefficients).collect();
    // The H¹ Weil polynomial mod 2 divides T^6-1, and the node cycles have
    // order dividing 3. The huge exponent is 3 modulo 6.
    assert_eq!(values[0], values[1]);
    assert_eq!(values[0], [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]);
    assert!(
        diagnostic(interpret(&format!("{F3} rvKummerHilbTrace[C,2,{huge},4];")))
            .contains("exact intermediate exceeds 512 bits")
    );
}

#[test]
fn grading_trivial_cases_and_nominal_prime_scalar() {
    let wide = u128::MAX;
    let output = stdout(interpret(&format!(
        "{F3}
        rvKummerHilbChar[C,0,{wide},0,3,4];
        rvKummerHilbChar[C,2,{wide},3,2,4];
        rvKummerHilbChar[C,2,{wide},10,2,4];
        rvKummerHilbChar[C,2,1,8,2,4];
        rvKummerHilbChar[C,{wide},1,2,2,8];
        rvKummerHilbChar[C,2,1,2,2,8];
        rvKummerHilbChar[C,2,{wide},4,2,0];
        rvKummerHilbChar[C,2,{wide},8,5,4];
        assert rvKummerHilbCharCoeff[C,0,1,0,2,0]==rvPrimeElt[2,1];
        assert rvKummerHilbCharCoeff[C,2,1,4,2,8]==rvPrimeElt[2,1];"
    )));
    let values: Vec<_> = output.lines().map(coefficients).collect();
    assert_eq!(values[0], [1, 2, 0, 0, 0]);
    assert_eq!(values[1], [1, 0, 0, 0, 0]);
    assert_eq!(values[2], values[1]);
    assert_eq!(values[3], [1, 1, 0, 0, 0]);
    assert_eq!(values[4], values[5]);
    assert_eq!(values[6], [1]);
    assert_eq!(values[7], [1, 4, 0, 0, 0]);
}

#[test]
fn invalid_moduli_extension_degrees_and_resource_requests_are_diagnosed() {
    for (source, expected) in [
        (
            "rvKummerHilbChar[C,1,0,2,2,3];",
            "extension degree must be positive",
        ),
        (
            "show rvKummerHilbCharCoeff[C,1,0,2,2,3];",
            "extension degree must be positive",
        ),
        ("rvKummerHilbChar[C,1,1,2,4,0];", "modulus must be prime"),
        ("rvKummerHilbChar[C,1,1,2,9,0];", "9 is not prime"),
        (
            "rvKummerHilbChar[C,1,1,2,2,22000];",
            "exceeds the output limit",
        ),
        (
            "rvKummerHilbChar[C,1000000000,1,2000000000,2,1];",
            "coefficient storage limit",
        ),
        (
            "rvKummerHilbChar[Missing,0,1,0,2,0];",
            "unknown residue hyperelliptic curve",
        ),
        (
            "show add(rvKummerHilbCharCoeff[C,2,1,4,3,0],rvPrimeElt[2,1]);",
            "expected type PrimeScalar(3)",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(problem.contains(expected), "{source} {problem}");
    }
}
