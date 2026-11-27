mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

fn values(source: &str) -> Vec<i128> {
    stdout(interpret(source))
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

fn eta_inverse() -> [i128; 10] {
    // Separate 24-coloured partition implementation for the raw polar oracle.
    let mut result = [0; 10];
    result[0] = 1;
    for k in 1..result.len() {
        for _ in 0..24 {
            for n in k..result.len() {
                result[n] += result[n - k];
            }
        }
    }
    result
}

fn direct_polar(m: i32, n: i32, r: i32, inside: bool, eta: &[i128; 10]) -> i128 {
    let mut value = 0;
    // The tested indices have |r| ≤ 7; |s| ≥ 20 cannot contribute.
    for s in -20..=20 {
        let t = r - 2 * m * s;
        let positive_expansion = s > 0 || (s == 0 && inside);
        if (positive_expansion && t > 0) || (!positive_expansion && t < 0) {
            let degree = n - m * s * s - s * t;
            if degree >= -1 {
                value += i128::from(t.abs()) * eta[(degree + 1) as usize];
            }
        }
    }
    eta[(m + 1) as usize] * value
}

#[test]
fn the_raw_appell_lerch_orbit_separates_both_chambers() {
    let eta = eta_inverse();
    let indices = [
        (1, -1, 0),
        (1, -1, 1),
        (1, -1, -3),
        (1, 0, 2),
        (1, 1, 3),
        (1, 1, -3),
        (1, 2, 4),
        (1, 5, 5),
        (1, 5, -5),
        (2, -1, 2),
        (2, 5, 7),
        (2, 5, -7),
    ];
    let mut source = "rvIgusaInverse[2,5,0,inside];".to_owned();
    for (m, n, r) in indices {
        writeln!(source, "rvIgusaFinite[{m},{n},{r}]; rvIgusaPolar[{m},{n},{r},inside]; rvIgusaPolar[{m},{n},{r},outside]; rvIgusaInverse[{m},{n},{r},inside]; rvIgusaInverse[{m},{n},{r},outside]; rvIgusaWallJump[{m},{n},{r}];").unwrap();
    }
    let actual = values(&source);
    let rows = actual[1..].chunks_exact(6);
    assert_eq!(rows.len(), indices.len());
    for ((m, n, r), row) in indices.into_iter().zip(rows) {
        let [finite, polar_inside, polar_outside, full_inside, full_outside, jump]: [i128; 6] =
            row.try_into().unwrap();
        assert_eq!(polar_inside, direct_polar(m, n, r, true, &eta));
        assert_eq!(polar_outside, direct_polar(m, n, r, false, &eta));
        assert_eq!(polar_inside - polar_outside, jump);
        assert_eq!(full_inside - polar_inside, finite, "inside ({m},{n},{r})");
        assert_eq!(
            full_outside - polar_outside,
            finite,
            "outside ({m},{n},{r})"
        );
    }
}

#[test]
fn normalized_small_indices_match_dabholkar_murthy_zagier() {
    // DMZ (128), (249) give Δψ₁ᶠ = 3 E₄ φ₋₂,₁ − 648 ℋ.
    // H(0)=-1/12, H(3)=1/3, H(4)=1/2, H(7)=H(8)=1.
    let mut source = String::new();
    for (m, max_n) in [(1, 2), (2, 1)] {
        for n in 0..=max_n {
            for r in 0..=m {
                writeln!(source, "rvIgusaMockWeight2[{m},{n},{r}];").unwrap();
            }
        }
    }
    // At index 2 these are obtained independently by combining (128)
    // with the three coefficient rows in DMZ §9.2.
    assert_eq!(
        values(&source),
        [48, 3, -1800, 528, -22416, 11709, 648, 72, 4, -28352, 6648, 1056]
    );

    // The normalized first two rows also determine the corresponding
    // unnormalized rows through Δ⁻¹ = q⁻¹ + 24 + ... .
    assert_eq!(
        values(
            r"
        rvIgusaFinite[1,-1,0]; rvIgusaFinite[1,-1,1];
        rvIgusaFinite[1,0,0]; rvIgusaFinite[1,0,1];
        rvIgusaFinite[2,-1,0]; rvIgusaFinite[2,-1,1]; rvIgusaFinite[2,-1,2];
        rvIgusaFinite[1,-1,2]; rvIgusaMockWeight2[1,-1,0];
    "
        ),
        [48, 3, -648, 600, 648, 72, 4, 0, 0]
    );
}

#[test]
fn spectral_flow_and_elliptic_evenness_use_the_same_reduced_coefficients() {
    let cases = [(1, -1, 1), (1, 1, 0), (2, 0, 2), (2, 2, 0), (3, 1, 2)];
    let shifts: [i128; 3] = [-4, 3, 1_000_000_000_000_000_000];
    for name in ["rvIgusaFinite", "rvIgusaMockWeight2"] {
        let mut source = String::new();
        for (m, n, r) in cases {
            writeln!(source, "{name}[{m},{n},{r}]; {name}[{m},{n},{}];", -r).unwrap();
            for shift in shifts {
                let exponent = n + shift * r + m * shift * shift;
                let power = r + 2 * m * shift;
                writeln!(source, "{name}[{m},{exponent},{power}];").unwrap();
            }
        }
        for group in values(&source).chunks_exact(2 + shifts.len()) {
            assert!(group.iter().all(|&value| value == group[0]), "{name}");
        }
    }
    // The lower shear with entry -2 transports the outside cusp period
    // and its Fourier index by exactly the λ=2 elliptic shift.
    assert_eq!(
        values(
            r"
            rvsiegelchamber Base[2,-1,3]; rvsiegelchamber Moved[18,-7,3];
            rvIgusaChamberCoeff[Base,1,-1,1]; rvIgusaChamberCoeff[Moved,1,5,5];
            rvIgusaFinite[1,5,5]; rvIgusaInverse[1,5,5,outside]; rvIgusaPolar[1,5,5,outside];
        "
        ),
        [3, 3, 3, 315255, 315252]
    );
}

#[test]
fn normalization_is_exact_convolution_by_delta() {
    let tau = [1, -24, 252, -1472];
    for (m, n, r) in [(1, 1, 0), (1, 2, 2), (2, 2, -3), (2, 3, 5), (3, 3, 6)] {
        let mut source = format!("rvIgusaMockWeight2[{m},{n},{r}];");
        for j in 0..=n {
            writeln!(source, "rvIgusaFinite[{m},{},{r}];", n - j - 1).unwrap();
        }
        let actual = values(&source);
        let expected: i128 = tau
            .into_iter()
            .zip(&actual[1..])
            .map(|(tau, coefficient)| tau * coefficient)
            .sum();
        assert_eq!(actual[0], expected, "({m},{n},{r})");
    }
}

#[test]
fn wide_laurent_indices_and_nonpositive_support_are_checked() {
    let hi = i128::MAX;
    let lo = i128::MIN;
    assert_eq!(
        stdout(interpret(&format!(
            "rvIgusaPolar[1,-1,{lo},outside]; rvIgusaPolar[1,-1,{lo},inside];\
         rvIgusaPolar[1,{hi},{lo},inside];\
         rvIgusaFinite[1,{hi},{lo}]; rvIgusaMockWeight2[1,{hi},{hi}];"
        ))),
        format!(
            "rvIgusaPolar[1, -1, {lo}, outside] = 55125743441192031081066686403946450255872\n\
         rvIgusaPolar[1, -1, {lo}, inside] = 0\n\
         rvIgusaPolar[1, {hi}, {lo}, inside] = 1323017842588608745945600473694714806125376\n\
         rvIgusaFinite[1, {hi}, {lo}] = 0\n\
         rvIgusaMockWeight2[1, {hi}, {hi}] = 0\n"
        )
    );
    assert_eq!(
        values("rvIgusaPolar[1,-2,7,inside]; rvIgusaFinite[1,-2,-7]; rvIgusaMockWeight2[1,-1,1];"),
        [0, 0, 0]
    );
    for (source, expected) in [
        ("rvIgusaPolar[0,0,0,inside];", "index m must be positive"),
        ("rvIgusaFinite[-1,0,0];", "index m must be positive"),
        ("rvIgusaMockWeight2[0,0,1];", "index m must be positive"),
        ("rvIgusaPolar[1,0,0,other];", "must be inside or outside"),
        (
            "rvIgusaFinite[100000,100000,0];",
            "coefficient storage limit",
        ),
        (
            "rvIgusaMockWeight2[100000,100000,0];",
            "coefficient storage limit",
        ),
        (
            "rvIgusaPolar[100000,0,1,inside];",
            "coefficient storage limit",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}

#[test]
fn the_runnable_example_subtracts_the_correct_polar_coefficients() {
    assert_eq!(
        stdout(example("igusa-mock-jacobi")),
        concat!(
            "rvIgusaInverse[1, 1, 3, inside] = 315255\n",
            "rvIgusaPolar[1, 1, 3, inside] = 315252\n",
            "rvIgusaInverse[1, 1, 3, outside] = 327\n",
            "rvIgusaPolar[1, 1, 3, outside] = 324\n",
            "rvIgusaFinite[1, 1, 3] = 3\n",
            "rvIgusaFinite[1, -1, 1] = 3\n",
            "rvIgusaMockWeight2[1, 1, 0] = -1800\n",
        )
    );
}
