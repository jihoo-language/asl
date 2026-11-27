mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

type Weight = (i32, i32, i32);
type Fourier = BTreeMap<Weight, i128>;

fn values(source: &str) -> Vec<i128> {
    stdout(interpret(source))
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn leading_and_axis_rows_use_the_requested_laurent_chamber() {
    let mut source = String::new();
    for r in -4..=4 {
        writeln!(
            source,
            "rvIgusaInverse[-1,-1,{r},inside]; rvIgusaInverse[-1,-1,{r},outside];"
        )
        .unwrap();
    }
    let leading: Vec<_> = (-4i128..=4).flat_map(|r| [r.max(0), (-r).max(0)]).collect();
    assert_eq!(values(&source), leading);

    // The first q-axis row is 2/y+20+2y. The second is
    // 3(y²+y⁻²)+42(y+y⁻¹)+234 by symmetric square and length two.
    // In bidegree (1,1), the cross-axis product plus the (1,1) oscillator
    // is 24(y²+y⁻²)−48(y+y⁻¹)+624.
    assert_eq!(
        values(
            r"
        rvIgusaInverse[-1,0,0,inside]; rvIgusaInverse[-1,0,0,outside];
        rvIgusaInverse[-1,0,1,inside]; rvIgusaInverse[-1,0,1,outside];
        rvIgusaInverse[1,-1,0,inside]; rvIgusaInverse[1,-1,0,outside];
        rvIgusaInverse[1,-1,1,inside]; rvIgusaInverse[1,-1,1,outside];
        rvIgusaInverse[0,0,0,inside]; rvIgusaInverse[0,0,0,outside];
        rvIgusaInverse[0,0,1,inside]; rvIgusaInverse[0,0,1,outside];
        rvIgusaInverse[0,0,2,inside]; rvIgusaInverse[0,0,2,outside];
    "
        ),
        [2, 2, 24, 0, 48, 48, 327, 3, 0, 0, 600, 24, 1152, 0]
    );
}

#[test]
fn the_jump_matches_the_independent_eta_residue_and_symmetries() {
    // 24-colored partitions by successive unbounded factors, independent
    // of the eta logarithmic-derivative recurrence and the K3 product.
    let mut eta = [0i128; 10];
    eta[0] = 1;
    for part in 1..eta.len() {
        for _ in 0..24 {
            for n in part..eta.len() {
                eta[n] += eta[n - part];
            }
        }
    }
    let mut source = String::from("rvIgusaResidue[-1,8]; rvIgusaInverse[2,2,0,inside];");
    let pairs = [(-1, -1), (-1, 0), (0, 0), (1, 0), (-1, 2), (2, 1), (2, 2)];
    let powers = [-12, -3, -1, 0, 1, 3, 12];
    for (m, n) in pairs {
        writeln!(source, "rvIgusaResidue[{m},{n}];").unwrap();
        for r in powers {
            writeln!(source, "rvIgusaInverse[{m},{n},{r},inside]; rvIgusaInverse[{m},{n},{r},outside]; rvIgusaWallJump[{m},{n},{r}]; rvIgusaInverse[{m},{n},{},outside]; rvIgusaInverse[{n},{m},{r},inside];", -r).unwrap();
        }
    }
    let output = values(&source);
    assert_eq!(output[0], eta[9]);
    let mut output = output[2..].iter().copied();
    for (m, n) in pairs {
        let residue = eta[(m + 1) as usize] * eta[(n + 1) as usize];
        assert_eq!(output.next(), Some(residue));
        for r in powers {
            let [inside, outside, jump, reflected, exchanged] =
                std::array::from_fn(|_| output.next().unwrap());
            assert_eq!(inside - outside, jump, "{m}, {n}, {r}");
            assert_eq!(jump, i128::from(r) * residue, "{m}, {n}, {r}");
            assert_eq!(inside, reflected);
            assert_eq!(inside, exchanged);
            if r == 12 {
                assert_eq!((inside, outside), (jump, 0));
            }
            if r == 0 {
                assert_eq!(inside, outside);
            }
        }
    }
    assert_eq!(output.next(), None);
}

#[test]
fn convolution_with_phi10_is_one_in_each_chamber() {
    let phi_indices: Vec<_> = (1..=3)
        .flat_map(|p| (1..=3).flat_map(move |q| (1 - p - q..p + q).map(move |r| (p, q, r))))
        .collect();
    let mut source = String::new();
    for &(p, q, r) in &phi_indices {
        writeln!(source, "rvIgusaPhi10[{p},{q},{r}];").unwrap();
    }
    let phi: Fourier = phi_indices
        .into_iter()
        .zip(values(&source))
        .filter(|(_, c)| *c != 0)
        .collect();
    let inverse_indices: Vec<_> = (-1..=1)
        .flat_map(|p| (-1..=1).flat_map(move |q| (-18..=18).map(move |r| (p, q, r))))
        .collect();
    for chamber in ["inside", "outside"] {
        let mut source = format!("rvIgusaInverse[1,1,0,{chamber}];");
        for &(p, q, r) in &inverse_indices {
            writeln!(source, "rvIgusaInverse[{p},{q},{r},{chamber}];").unwrap();
        }
        let inverse: Fourier = inverse_indices
            .iter()
            .copied()
            .zip(values(&source).into_iter().skip(1))
            .collect();
        for p in 0..=2 {
            for q in 0..=2 {
                for r in -12..=12 {
                    let product: i128 = phi
                        .iter()
                        .filter(|&(&(a, b, _), _)| a <= p + 1 && b <= q + 1)
                        .map(|(&(a, b, s), c)| c * inverse[&(p - a, q - b, r - s)])
                        .sum();
                    assert_eq!(
                        product,
                        i128::from((p, q, r) == (0, 0, 0)),
                        "{chamber} at ({p},{q},{r})"
                    );
                }
            }
        }
    }
}

#[test]
fn arbitrarily_remote_laurent_tails_use_wide_arithmetic() {
    let hi = i128::MAX;
    let lo = i128::MIN;
    assert_eq!(
        stdout(interpret(&format!(
            "rvIgusaInverse[0,0,{lo},outside]; rvIgusaInverse[0,0,{hi},inside];\
         rvIgusaInverse[0,0,{hi},outside]; rvIgusaWallJump[0,0,{lo}];"
        ))),
        format!(
            "rvIgusaInverse[0, 0, {lo}, outside] = 98001321673230277477451886940349244899328\n\
         rvIgusaInverse[0, 0, {hi}, inside] = 98001321673230277477451886940349244898752\n\
         rvIgusaInverse[0, 0, {hi}, outside] = 0\n\
         rvIgusaWallJump[0, 0, {lo}] = -98001321673230277477451886940349244899328\n"
        )
    );
}

#[test]
fn diagnostics_and_the_example_observe_signed_support() {
    assert_eq!(
        stdout(example("igusa-meromorphic-chambers")),
        concat!(
            "rvIgusaInverse[-1, -1, 3, inside] = 3\n",
            "rvIgusaInverse[-1, -1, 3, outside] = 0\n",
            "rvIgusaInverse[0, 0, 1, inside] = 600\n",
            "rvIgusaInverse[0, 0, 1, outside] = 24\n",
            "rvIgusaResidue[0, 0] = 576\n",
            "rvIgusaWallJump[0, 0, 1] = 576\n",
        )
    );
    let hi = i128::MAX;
    assert_eq!(
        values(&format!(
            "rvIgusaInverse[-2,{hi},1,inside]; rvIgusaInverse[{hi},-2,-1,outside];\
         rvIgusaWallJump[-2,{hi},1]; rvIgusaWallJump[{hi},{hi},0]; rvIgusaResidue[-2,{hi}];"
        )),
        [0; 5]
    );
    for (source, expected) in [
        (
            "rvIgusaInverse[-1,-1,1,other];",
            "must be inside or outside",
        ),
        (
            "rvIgusaInverse[100000,100000,0,inside];",
            "coefficient storage limit",
        ),
        (
            "rvIgusaResidue[170141183460469231731687303715884105727,0];",
            "coefficient storage limit",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}
