mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

type Index = (u32, u32, i32);
type Fourier = BTreeMap<Index, i128>;

const DMVV: &str = "rvK3EllipticGenus";
const IGUSA: &str = "rvIgusaPhi10";

fn coefficients(name: &str, indices: impl IntoIterator<Item = Index>) -> Fourier {
    let indices = indices.into_iter().collect::<Vec<_>>();
    let mut source = String::new();
    for &(p, q, r) in &indices {
        writeln!(source, "{name}[{p},{q},{r}];").unwrap();
    }
    let output = stdout(interpret(&source));
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), indices.len());
    indices
        .into_iter()
        .zip(lines)
        .map(|((p, q, r), line)| {
            let prefix = format!("{name}[{p}, {q}, {r}] = ");
            let coefficient = line.strip_prefix(&prefix).unwrap().parse().unwrap();
            ((p, q, r), coefficient)
        })
        .collect()
}

fn slice(name: &str, p: u32, q: u32, radius: i32) -> Vec<i128> {
    coefficients(name, (-radius..=radius).map(|r| (p, q, r)))
        .into_values()
        .collect()
}

#[test]
fn theta_seed_matches_the_primitive_k3_expansion_and_specialization() {
    assert_eq!(slice(DMVV, 1, 0, 2), [0, 2, 20, 2, 0]);
    assert_eq!(slice(DMVV, 1, 1, 3), [0, 20, -128, 216, -128, 20, 0]);
    assert_eq!(slice(DMVV, 1, 2, 3), [2, 216, -1026, 1616, -1026, 216, 2]);
    for q in 1..=5_u32 {
        let radius = (4 * q + 1).isqrt() as i32;
        assert_eq!(slice(DMVV, 1, q, radius).into_iter().sum::<i128>(), 0);
    }
    assert_eq!(slice(DMVV, 0, 0, 1), [0, 1, 0]);
    assert_eq!(slice(DMVV, 0, 3, 1), [0, 0, 0]);
}

#[test]
fn the_length_two_and_three_coefficients_match_independent_expansions() {
    // At q^0 length two has Sym²(c(0,y)) plus the m=2 factor.
    let second = slice(DMVV, 2, 0, 2);
    assert_eq!(second, [3, 42, 234, 42, 3]);
    assert_eq!(second.into_iter().sum::<i128>(), 324);
    let third = slice(DMVV, 3, 0, 3);
    assert_eq!(third.into_iter().sum::<i128>(), 3200);

    // At q^1 length two has c(0,y)c(1,y)+c(2,y).
    assert_eq!(slice(DMVV, 2, 1, 3), [42, 360, -3114, 5424, -3114, 360, 42]);
}

#[test]
fn the_hilbert_genus_obeys_jacobi_elliptic_translation() {
    let shifts = [
        (1, 0, 0, 1),
        (1, 1, 1, 1),
        (2, 0, -2, 1),
        (2, 0, 0, 1),
        (2, 0, 2, 1),
        (2, 1, -1, -1),
        (2, 2, 4, -1),
        (3, 0, 0, -1),
        (3, 1, -2, 1),
    ];
    let pairs = shifts.map(|(n, d, r, shift)| {
        let translated = d + shift * r + n * shift * shift;
        (
            (n as u32, d as u32, r),
            (n as u32, translated as u32, r + 2 * n * shift),
        )
    });
    let values = coefficients(DMVV, pairs.into_iter().flat_map(|(a, b)| [a, b]));
    for (original, translated) in pairs {
        assert_eq!(values[&original], values[&translated], "{original:?}");
    }
}

#[test]
fn the_igusa_weyl_normalization_and_first_product_slices_are_exact() {
    assert_eq!(slice(IGUSA, 1, 1, 2), [0, 1, -2, 1, 0]);
    // The next coefficient is the Weyl Laurent factor times -c(0,y).
    assert_eq!(slice(IGUSA, 1, 2, 2), [-2, -16, 36, -16, -2]);
    // At P²q² the remaining product contributes c(0,y)²-c(1,y).
    assert_eq!(
        slice(IGUSA, 2, 2, 4),
        [0, -16, 240, -240, 32, -240, 240, -16, 0]
    );
    assert_eq!(slice(IGUSA, 0, 2, 1), [0, 0, 0]);
    assert_eq!(slice(IGUSA, 3, 0, 1), [0, 0, 0]);
}

#[test]
fn igusa_is_siegel_symmetric_and_cancels_the_positive_length_dmv_product() {
    let rows = |p: u32, q: u32| (-7..=7).map(move |r| (p, q, r));
    let phi = coefficients(
        IGUSA,
        (1..=3).flat_map(|p| (1..=3).flat_map(move |q| rows(p, q))),
    );
    let genus = coefficients(
        DMVV,
        (0..=2).flat_map(|p| (0..=2).flat_map(move |q| rows(p, q))),
    );
    for (&(p, q, r), &value) in &phi {
        assert_eq!(value, phi[&(q, p, r)]);
        assert_eq!(value, phi[&(p, q, -r)]);
        if r * r >= 4 * (p * q) as i32 {
            assert_eq!(value, 0, "nonpositive discriminant at {:?}", (p, q, r));
        }
    }
    // Radius seven covers every coefficient in these retained slices
    // before cancellation. Radius fourteen covers the signed convolution.
    for p in 2..=3 {
        for q in 1..=3 {
            for r in -14..=14 {
                let coefficient: i128 = phi
                    .iter()
                    .filter(|&(&(m, n, _), _)| m <= p && n <= q)
                    .map(|(&(m, n, s), &value)| {
                        value * genus.get(&(p - m, q - n, r - s)).unwrap_or(&0)
                    })
                    .sum();
                assert_eq!(coefficient, 0, "{:?}", (p, q, r));
            }
        }
    }
}

#[test]
fn the_example_runs_and_indices_are_checked() {
    assert_eq!(
        stdout(example("k3-elliptic-genus")),
        concat!(
            "rvK3EllipticGenus[1, 1, -1] = -128\n",
            "rvK3EllipticGenus[2, 0, 0] = 234\n",
            "rvK3EllipticGenus[3, 0, 0] = 2048\n",
            "rvIgusaPhi10[1, 1, 0] = -2\n",
            "rvIgusaPhi10[1, 2, -1] = -16\n",
            "rvIgusaPhi10[2, 2, 0] = 32\n",
        )
    );
    for (source, expected) in [
        (
            "rvK3EllipticGenus[-1,0,0];",
            "a nonnegative decimal integer",
        ),
        ("rvIgusaPhi10[1,-2,0];", "a nonnegative decimal integer"),
        (
            "rvIgusaPhi10[1,1,-170141183460469231731687303715884105729];",
            "y exponent exceeds the signed i128 range",
        ),
        (
            "rvK3EllipticGenus[1000000,1000000,0];",
            "coefficient storage limit",
        ),
        (
            "rvIgusaPhi10[1000000,1000000,0];",
            "coefficient storage limit",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{error}");
    }
    assert_eq!(
        stdout(interpret(
            "rvIgusaPhi10[1,1,-170141183460469231731687303715884105728];"
        )),
        "rvIgusaPhi10[1, 1, -170141183460469231731687303715884105728] = 0\n"
    );
}
