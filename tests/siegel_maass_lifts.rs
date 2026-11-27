mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

const MIXTURE: &str = "rvsiegelmaass Mix[16]{zero(1,0,1,3),minus(0,1,1,-2)};";

fn values(source: &str) -> Vec<i128> {
    stdout(interpret(source))
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn the_weight_ten_lift_agrees_with_the_borcherds_product() {
    let mut source = "rvsiegelmaass A[10]{minus(0,0,1,1)};".to_owned();
    for m in 1..=3 {
        for n in 1..=3 {
            let radius = (4_u32 * m * n).isqrt() as i32 + 1;
            for r in -radius..=radius {
                writeln!(
                    source,
                    "rvSiegelMaassCoeff[A,{m},{n},{r}]; rvIgusaPhi10[{m},{n},{r}];"
                )
                .unwrap();
            }
        }
    }
    source.push_str("rvSiegelMaassCoeff[A,4,4,0]; rvIgusaPhi10[4,4,0];");
    let values = values(&source);
    assert!(values.chunks_exact(2).all(|pair| pair[0] == pair[1]));
    assert_eq!(&values[values.len() - 2..], [-279040, -279040]);
}

#[test]
fn chi12_is_normalized_and_has_the_witt_restriction() {
    let mut source = "rvsiegelmaass C[12]{zero(0,0,1,1)};".to_owned();
    for r in -2..=2 {
        writeln!(
            source,
            "rvIgusaChi12[1,1,{r}]; rvSiegelMaassCoeff[C,1,1,{r}];"
        )
        .unwrap();
    }
    assert_eq!(values(&source), [0, 0, 1, 1, 10, 10, 1, 1, 0, 0]);

    // These Ramanujan tau values come from the ordinary Δ series.
    let tau = [0, 1, -24, 252, -1472, 4830];
    for (m, n) in [(1, 2), (1, 5), (2, 2), (2, 3), (3, 3), (3, 4), (4, 5)] {
        let radius = (4_u32 * m * n - 1).isqrt() as i32;
        let mut source = String::new();
        for r in -radius..=radius {
            writeln!(source, "rvIgusaChi12[{m},{n},{r}];").unwrap();
        }
        assert_eq!(
            values(&source).into_iter().sum::<i128>(),
            12 * tau[m as usize] * tau[n as usize],
            "Witt coefficient ({m}, {n})"
        );
    }
    assert_eq!(
        values("rvIgusaChi12[0,1,0]; rvIgusaChi12[1,1,2]; rvIgusaChi12[1,1,-3];"),
        [0, 0, 0]
    );
}

#[test]
fn mixed_weight_seeds_evaluate_both_jacobi_families_and_signed_powers() {
    let mut source = MIXTURE.to_owned();
    for (n, radius) in [(1, 2), (2, 2)] {
        for r in -radius..=radius {
            writeln!(source, "rvSiegelMaassJacobi[Mix,{n},{r}];").unwrap();
        }
    }
    // Expand E4=1+240q+…, E6=1−504q+…, Δ=q−24q²+…
    // against the first two Laurent rows of the two weak Jacobi forms.
    assert_eq!(values(&source), [0, 1, 34, 1, 0, 34, 1496, 4716, 1496, 34]);

    // E4³-E6² = 1728 Δ. This also tests a higher power of Δ and
    // cancellation between distinct monomials before a first nonzero row.
    let mut source = concat!(
        "rvsiegelmaass E[24]{zero(3,0,1,1),zero(0,2,1,-1)};",
        "rvsiegelmaass D[24]{zero(0,0,2,1728)};",
    )
    .to_owned();
    for (n, r) in [(1, 0), (2, 0), (2, -1), (3, 2), (4, 0), (8, -5)] {
        writeln!(
            source,
            "rvSiegelMaassJacobi[E,{n},{r}]; rvSiegelMaassJacobi[D,{n},{r}];"
        )
        .unwrap();
    }
    let values = values(&source);
    assert!(values.chunks_exact(2).all(|pair| pair[0] == pair[1]));
    assert_eq!(&values[..6], [0, 0, 17280, 17280, 1728, 1728]);
}

#[test]
fn maass_relations_and_large_unimodular_basis_changes_preserve_coefficients() {
    let mut source = MIXTURE.to_owned();
    source.push_str("rvSiegelMaassCoeff[Mix,4,4,4]; rvSiegelMaassJacobi[Mix,16,4]; rvSiegelMaassJacobi[Mix,4,2]; rvSiegelMaassJacobi[Mix,1,1];");
    source.push_str("rvSiegelMaassCoeff[Mix,6,4,-2]; rvSiegelMaassJacobi[Mix,24,-2]; rvSiegelMaassJacobi[Mix,6,-1];");
    let output = values(&source);
    assert_eq!(
        output[0],
        output[1] + 2_i128.pow(15) * output[2] + 4_i128.pow(15) * output[3]
    );
    assert_eq!(output[4], output[5] + 2_i128.pow(15) * output[6]);

    // Apply the determinant-one change of variables with rows (k,k−1)
    // and (k+1,k) to x²+y², then also test content two and a rank-one form.
    let k = 1_000_000_000_000_000_000_i128;
    let (m, n, r) = (
        k * k + (k + 1) * (k + 1),
        k * k + (k - 1) * (k - 1),
        4 * k * k,
    );
    let (a, b, s) = (k * k, (k - 1) * (k - 1), 2 * k * (k - 1));
    let output = values(&format!(
        "rvIgusaChi12[{m},{n},{r}]; rvIgusaChi12[{n},{m},-{r}]; rvIgusaChi12[{},{},{}]; rvIgusaChi12[2,2,0];\
         rvIgusaChi12[{a},{b},{s}]; {MIXTURE}\
         rvSiegelMaassJacobi[Mix,{a},{}]; rvSiegelMaassJacobi[Mix,{},{}]; rvSiegelMaassJacobi[Mix,1,0];",
        2*m, 2*n, 2*r, 2*k, a+1, 2*k
    ));
    assert_eq!(&output[..2], [10, 10]);
    assert_ne!(output[2], 0);
    assert_eq!(output[2], output[3]);
    assert_eq!(&output[4..], [0, 0, 34, 34]);
}

#[test]
fn declarations_validate_weights_and_resolve_exact_cancellation() {
    for (source, expected) in [
        ("rvsiegelmaass X[11]{};", "even weight at least 10"),
        ("rvsiegelmaass X[8]{};", "even weight at least 10"),
        ("rvsiegelmaass X[10]{minus(0,0,0,1)};", "positive Delta exponent"),
        ("rvsiegelmaass X[12]{minus(0,0,1,1)};", "does not have the declared weight"),
        ("rvsiegelmaass X[10]{cusp(0,0,1,1)};", "expected a minus or zero Jacobi monomial"),
        ("rvsiegelmaass X[10]{minus(-1,0,1,1)};", "nonnegative decimal integer"),
        ("rvsiegelmaass X[12]{}; rvsiegelmaass X[14]{};", "already defined"),
        ("rvSiegelMaassCoeff[Missing,0,0,0];", "unknown Maaß lift"),
        ("rvIgusaChi12[340282366920938463463374607431768211455,340282366920938463463374607431768211455,0];", "coefficient storage limit"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
    let huge = i128::MAX;
    let max = u128::MAX;
    assert_eq!(values(&format!(
        "rvsiegelmaass Z[10]{{minus(0,0,1,{huge}),minus(0,0,1,{huge}),minus(0,0,1,-{huge}),minus(0,0,1,-{huge})}};\
         rvSiegelMaassCoeff[Z,{max},{max},0]; rvSiegelMaassJacobi[Z,{max},0];\
         rvsiegelmaass High[120000000]{{zero(0,0,10000000,1)}}; rvSiegelMaassCoeff[High,2,2,0];"
    )), [0, 0, 0]);
    assert_eq!(
        stdout(interpret(&format!(
            "rvsiegelmaass W[12]{{zero(0,0,1,{huge}),zero(0,0,1,{huge})}}; rvSiegelMaassCoeff[W,1,1,1];"
        ))),
        format!("rvSiegelMaassCoeff[W, 1, 1, 1] = {}\n", max - 1)
    );
    let overflow = diagnostic(interpret(
        "rvsiegelmaass W[4012]{zero(1000,0,1,1)}; rvSiegelMaassCoeff[W,2,2,2];",
    ));
    assert!(
        overflow.contains("exact intermediate exceeds 512 bits"),
        "{overflow}"
    );
}

#[test]
fn the_example_runs_with_the_published_normalizations() {
    assert_eq!(
        stdout(example("siegel-maass-lifts")),
        concat!(
            "rvSiegelMaassJacobi[Lift10, 1, 0] = -2\n",
            "rvSiegelMaassCoeff[Lift10, 4, 4, 0] = -279040\n",
            "rvIgusaPhi10[4, 4, 0] = -279040\n",
            "rvIgusaChi12[1, 1, 0] = 10\n",
            "rvSiegelMaassJacobi[Lift16, 2, 0] = 4716\n",
            "rvSiegelMaassCoeff[Lift16, 1, 2, -1] = 1496\n",
        )
    );
}
