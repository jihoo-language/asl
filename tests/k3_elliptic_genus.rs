mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

type Index = (u32, u32, i32);
type Fourier = BTreeMap<Index, i128>;

const DMVV: &str = "rvK3EllipticGenus";
const IGUSA: &str = "rvIgusaPhi10";
const SYMPLECTIC: &str = "rvK3SymplecticElliptic";
const TRANSLATE: &str = "rvKummerTranslateElliptic";
const ENRIQUES: &str = "rvKummerEnriquesElliptic";

const MARKING: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[L](zw,zw);
    let f1=rvJacPoint[L](ow,zw);
    let e2=add(add(e1,f1),rvJacPoint[L](neg(ow),zw));
    let f2=add(add(e1,f1),rvJacPoint[L](i,zw));
    rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};
";

fn coefficients(name: &str, indices: impl IntoIterator<Item = Index>) -> Fourier {
    report_coefficients("", name, "", indices)
}

fn report_coefficients(
    setup: &str,
    name: &str,
    marking: &str,
    indices: impl IntoIterator<Item = Index>,
) -> Fourier {
    let indices = indices.into_iter().collect::<Vec<_>>();
    let mut source = setup.to_owned();
    for &(p, q, r) in &indices {
        writeln!(source, "{name}[{marking}{p},{q},{r}];").unwrap();
    }
    let output = stdout(interpret(&source));
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), indices.len());
    indices
        .into_iter()
        .zip(lines)
        .map(|((p, q, r), line)| {
            let prefix = format!("{name}[{marking}{p}, {q}, {r}] = ");
            let coefficient = line.strip_prefix(&prefix).unwrap().parse().unwrap();
            ((p, q, r), coefficient)
        })
        .collect()
}

fn twined_slice(name: &str, marking: &str, n: u32, d: u32, radius: i32) -> Vec<i128> {
    report_coefficients(
        MARKING,
        name,
        marking,
        (-radius..=radius).map(|r| (n, d, r)),
    )
    .into_values()
    .collect()
}

fn slice(name: &str, p: u32, q: u32, radius: i32) -> Vec<i128> {
    coefficients(name, (-radius..=radius).map(|r| (p, q, r)))
        .into_values()
        .collect()
}

fn symplectic(order: u32, indices: impl IntoIterator<Item = Index>) -> Fourier {
    report_coefficients("", SYMPLECTIC, &format!("{order}, "), indices)
}

#[test]
fn cyclic_primitive_seeds_match_independent_jacobi_forms() {
    let fixed = [24, 8, 6, 4, 4, 2, 3, 2];
    for (h, fixed_count) in (1..=8).zip(fixed) {
        let indices = (0..=1).flat_map(|d| (-2..=2).map(move |r| (1, d, r)));
        let coefficients = symplectic(h, indices);
        assert_eq!(
            (-2..=2)
                .map(|r| coefficients[&(1, 0, r)])
                .collect::<Vec<_>>(),
            [0, 2, fixed_count - 4, 2, 0],
            "order {h}"
        );
        let q = (-2..=2)
            .map(|r| coefficients[&(1, 1, r)])
            .collect::<Vec<_>>();
        let k3 = [20, -128, 216, -128, 20];
        let expected = match h {
            1 => k3,
            2 | 3 | 5 | 7 => {
                // Prime order p: (2φ₀,₁ + 2p E₂⁽ᵖ⁾ φ₋₂,₁)/(p+1).
                // Here 2φ₀,₁ is the ordinary K3 genus, and E₂⁽ᵖ⁾ has
                // constant 1 and q coefficient 24/(p-1).
                let q0 = [0, 1, -2, 1, 0];
                let q1 = [-2, 8, -12, 8, -2];
                let p = i128::from(h);
                std::array::from_fn(|i| (k3[i] + 2 * p * (q1[i] + 24 / (p - 1) * q0[i])) / (p + 1))
            }
            // The first oscillator coefficient is x(2-y-y⁻¹).
            // For h=8, Tr(x/(2-x))=2 and Tr(x²/(2-x))=4 over Q(√2).
            4 => [0; 5],
            6 => [-2, 6, -8, 6, -2],
            8 => [-2, 8, -12, 8, -2],
            _ => unreachable!(),
        };
        assert_eq!(q, expected, "order {h}");
    }
}

#[test]
fn cyclic_hilbert_euler_traces_match_frame_shape_partitions() {
    let frames: [&[(usize, usize)]; 8] = [
        &[(1, 24)],
        &[(1, 8), (2, 8)],
        &[(1, 6), (3, 6)],
        &[(1, 4), (2, 2), (4, 4)],
        &[(1, 4), (5, 4)],
        &[(1, 2), (2, 2), (3, 2), (6, 2)],
        &[(1, 3), (7, 3)],
        &[(1, 2), (2, 1), (4, 1), (8, 2)],
    ];
    for (h, frame) in (1..=8).zip(frames) {
        let euler = symplectic(
            h,
            (0..=8).flat_map(|n| (-(n as i32)..=n as i32).map(move |r| (n, 0, r))),
        );
        let mut expected = [0; 9];
        expected[0] = 1;
        // Independently treat each colored cycle of length k as a freely
        // repeatable part of size km. This Euler calculation has no q or y.
        for &(k, colors) in frame {
            for part in (k..=8).step_by(k) {
                for _ in 0..colors {
                    for n in part..=8 {
                        expected[n] += expected[n - part];
                    }
                }
            }
        }
        let mut traces = [0; 9];
        for ((n, _, _), value) in euler {
            traces[n as usize] += value;
        }
        assert_eq!(traces, expected, "order {h}");
        match h {
            3 => assert_eq!(traces[2..=3], [27, 104]),
            5 => assert_eq!(traces[2], 14),
            7 => assert_eq!(traces[2], 9),
            _ => (),
        }
    }
}

#[test]
fn cyclic_orders_one_and_two_match_the_existing_genera() {
    let indices = [
        (0, 1, 0),
        (1, 2, -1),
        (1, 4, 3),
        (2, 2, -4),
        (3, 2, 0),
        (4, 1, 2),
    ];
    assert_eq!(symplectic(1, indices), coefficients(DMVV, indices));
    assert_eq!(
        symplectic(2, indices),
        report_coefficients(MARKING, TRANSLATE, "B, 9, ", indices)
    );
}

#[test]
fn order_three_hilbert_square_matches_the_geometric_lefschetz_trace() {
    let primitive = symplectic(3, (0..=4).flat_map(|d| (-4..=4).map(move |r| (1, d, r))));
    let square = symplectic(3, (0..=4).flat_map(|d| (-6..=6).map(move |r| (2, d, r))));
    for ((_, d, r), trace) in square {
        // There are C(6,2) reduced fixed pairs and two invariant double
        // points at each of the six fixed points. For a double point the
        // symplectic weight pairs are ζ^±1 and ζ^±2, which agree at order 3.
        // Thus all 27 isolated Hilbert-square points have local factor
        // (Z_g/6)^2, giving 4 Z_(g^[2]) = 3 Z_g^2.
        let product: i128 = primitive
            .iter()
            .filter(|&(&(_, q, _), _)| q <= d)
            .map(|(&(_, q, s), value)| value * primitive.get(&(1, d - q, r - s)).unwrap_or(&0))
            .sum();
        assert_eq!(4 * trace, 3 * product, "q^{d} y^{r}");
    }
}

#[test]
fn cyclic_hilbert_coefficients_are_jacobi_and_use_all_power_traces() {
    for h in 2..=8 {
        let values = symplectic(
            h,
            [
                (1, 0, 0),
                (1, 1, 0),
                (1, 2, 0),
                (2, 0, 0),
                (2, 1, 0),
                (2, 3, -4),
                (2, 2, 4),
                (2, 0, 4),
            ],
        );
        let primitive = &values[&(1, 0, 0)];
        let g_squared = if h.is_multiple_of(2) { h / 2 } else { h };
        let trace_squared = symplectic(g_squared, [(1, 0, 0)])[&(1, 0, 0)];
        // At q⁰ the constant term of the square of the Laurent trace is
        // the square of its constant plus 2·2². The symmetric square uses
        // the trace of g² at y², and the length-two oscillator adds c_g.
        assert_eq!(
            values[&(2, 0, 0)],
            (primitive * primitive + 8 + trace_squared) / 2 + primitive,
            "order {h}"
        );
        // At q¹ there is a product of two m=1 traces and one m=2 trace.
        let side = symplectic(h, [(1, 1, 1)])[&(1, 1, 1)];
        assert_eq!(
            values[&(2, 1, 0)],
            primitive * values[&(1, 1, 0)] + 4 * side + values[&(1, 2, 0)],
            "order {h}"
        );
        assert_eq!(values[&(2, 1, 0)], values[&(2, 3, -4)], "order {h}");
        assert_eq!(values[&(2, 0, 0)], values[&(2, 2, 4)], "order {h}");
        assert_eq!(values[&(2, 0, 4)], 0, "order {h}");
    }
}

#[test]
fn cyclic_syntax_examples_and_cutoff_diagnostics() {
    assert_eq!(
        stdout(example("k3-symplectic-elliptic-genus")),
        concat!(
            "rvK3SymplecticElliptic[1, 1, 0, 0] = 20\n",
            "rvK3SymplecticElliptic[2, 2, 0, 0] = 26\n",
            "rvK3SymplecticElliptic[3, 1, 1, -1] = -2\n",
            "rvK3SymplecticElliptic[3, 2, 0, 0] = 9\n",
            "rvK3SymplecticElliptic[8, 1, 1, 0] = -12\n",
        )
    );
    for (source, expected) in [
        (
            "rvK3SymplecticElliptic[0,1,0,0];",
            "symplectic K3 order must be between 1 and 8",
        ),
        (
            "rvK3SymplecticElliptic[9,1,0,0];",
            "symplectic K3 order must be between 1 and 8",
        ),
        (
            "rvK3SymplecticElliptic[3,1000000,1000000,0];",
            "coefficient storage limit",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{error}");
    }
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

#[test]
fn nikulin_primitive_and_hilbert_coefficients_match_lefschetz_and_symmetric_squares() {
    for (degree, radius, expected) in [
        (0, 1, vec![2, 4, 2]),
        (1, 2, vec![4, 0, -8, 0, 4]),
        (2, 3, vec![2, -8, -2, 16, -2, -8, 2]),
    ] {
        assert_eq!(
            twined_slice(TRANSLATE, "B, 1, ", 1, degree, radius),
            expected
        );
    }
    let second = twined_slice(TRANSLATE, "B, 1, ", 2, 0, 2);
    // Sym² uses half the squared trace plus half the trace of g² at y².
    // The length-two factor contributes the primitive g trace once more.
    assert_eq!(second, [3, 10, 26, 10, 3]);
    assert_eq!(second.into_iter().sum::<i128>(), 52);
    assert_eq!(
        twined_slice(TRANSLATE, "B, 1, ", 2, 1, 3),
        [10, 8, -10, -16, -10, 8, 10]
    );
    assert_eq!(
        twined_slice(TRANSLATE, "B, 1, ", 3, 0, 3)
            .into_iter()
            .sum::<i128>(),
        256
    );
}

#[test]
fn identity_and_nonzero_translation_depend_only_on_the_tame_action_type() {
    let indices = [(0, 0, 0), (1, 0, 0), (1, 1, -1), (2, 1, 0), (3, 0, 2)];
    assert_eq!(
        coefficients(DMVV, indices),
        report_coefficients(MARKING, TRANSLATE, "B, 0, ", indices)
    );
    let mut source = MARKING.to_owned();
    let mut expected = String::new();
    for mask in 1..16 {
        writeln!(source, "rvKummerTranslateElliptic[B,{mask},1,1,0];").unwrap();
        writeln!(
            expected,
            "rvKummerTranslateElliptic[B, {mask}, 1, 1, 0] = -8"
        )
        .unwrap();
    }
    source.push_str(
        r"
        rvunram V[3,1]{2,0,1,0,1};
        rvjacbase M[L,V];
        let ov=rvUnramElem[V]{1};
        rvjactorsionbasis E[M,2](neg(ov)){
            rvJacExtend[M](e1),rvJacExtend[M](e2),rvJacExtend[M](f1),rvJacExtend[M](f2)
        };
        rvKummerTranslateElliptic[E,8,1,1,0];
        rvKummerTranslateElliptic[E,13,2,0,0];
        ",
    );
    expected.push_str(concat!(
        "rvKummerTranslateElliptic[E, 8, 1, 1, 0] = -8\n",
        "rvKummerTranslateElliptic[E, 13, 2, 0, 0] = 26\n",
    ));
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn free_switch_uses_half_the_seed_at_doubled_indices() {
    let primitive = coefficients(DMVV, (0..=2).flat_map(|d| (-3..=3).map(move |r| (1, d, r))));
    let values = report_coefficients(
        MARKING,
        ENRIQUES,
        "B, C, 2, ",
        (0..=4).flat_map(|d| (-7..=7).map(move |r| (2, d, r))),
    );
    for (&(_, d, r), &value) in &values {
        let expected = if d.is_multiple_of(2) && r % 2 == 0 {
            primitive[&(1, d / 2, r / 2)] / 2
        } else {
            0
        };
        assert_eq!(value, expected, "q^{d} y^{r}");
    }
    assert_eq!(
        values
            .iter()
            .filter(|&(&(_, d, _), _)| d == 0)
            .map(|(_, value)| value)
            .sum::<i128>(),
        12
    );
    let fourth = twined_slice(ENRIQUES, "B, C, 2, ", 4, 0, 4);
    assert_eq!(fourth, [1, 0, 11, 0, 66, 0, 11, 0, 1]);
    assert_eq!(fourth.into_iter().sum::<i128>(), 90);
    let odd = report_coefficients(
        MARKING,
        ENRIQUES,
        "B, C, 2, ",
        [1, 3]
            .into_iter()
            .flat_map(|n| (0..=3).flat_map(move |d| (-5..=5).map(move |r| (n, d, r)))),
    );
    assert!(odd.values().all(|&value| value == 0));
    assert_eq!(twined_slice(ENRIQUES, "B, C, 2, ", 0, 0, 1), [0, 1, 0]);
}

#[test]
fn twining_requires_the_marked_translation_or_a_certified_free_switch() {
    assert_eq!(
        stdout(example("k3-twined-elliptic-genus")),
        concat!(
            "rvKummerTranslateElliptic[B, 0, 1, 0, 0] = 20\n",
            "rvKummerTranslateElliptic[B, 9, 1, 0, 0] = 4\n",
            "rvKummerTranslateElliptic[B, 9, 1, 1, 0] = -8\n",
            "rvKummerTranslateElliptic[B, 9, 2, 0, 0] = 26\n",
            "rvKummerEnriquesElliptic[B, C, 2, 1, 0, 0] = 0\n",
            "rvKummerEnriquesElliptic[B, C, 2, 2, 2, -2] = -64\n",
            "rvKummerEnriquesElliptic[B, C, 2, 4, 0, 0] = 66\n",
        )
    );
    // Mask 8 becomes Frobenius-stable after extension to the marking field.
    assert_eq!(
        twined_slice(ENRIQUES, "B, L, 8, ", 2, 0, 2),
        [1, 0, 10, 0, 1]
    );
    for (input, expected) in [
        (
            "rvKummerEnriquesElliptic[B,C,13,1,1,1];",
            "requires an even theta characteristic",
        ),
        (
            "rvKummerEnriquesElliptic[B,C,8,1,1,1];",
            "must be fixed by base Frobenius",
        ),
        ("rvKummerTranslateElliptic[B,16,1,0,0];", "between 0 and 15"),
        ("rvKummerEnriquesElliptic[B,C,16,0,0,0];", "between 0 and 15"),
        ("rvKummerTranslateElliptic[Unknown,0,1,0,0];", "unknown Jacobian torsion basis"),
        (
            "rvhyperelliptic A[F]{z,neg(o),z,o}; rvjactorsionbasis J[A,2](neg(o)){rvJacPoint[A](z,z),rvJacPoint[A](o,z)}; rvKummerTranslateElliptic[J,1,1,0,0];",
            "certified genus-two level-2 basis",
        ),
        (
            "rvhyperelliptic D[F]{z,o,z,z,z,o}; rvKummerEnriquesElliptic[B,D,2,2,0,0];",
            "not the declared base extension",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{MARKING}{input}")));
        assert!(error.contains(expected), "{error}");
    }
}
