mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_fermat_and_non_fermat_series_use_the_supplied_integral_lift() {
    let output = stdout(example("dwork-k3-formal-brauer"));
    assert!(output.contains("rvDworkK3Log[5, 0, 9] = [1: 1, 5: 24/5, 9: 280]"));
    assert!(output.contains("(2, 3): -48"));
    assert!(output.contains("rvDworkK3Inverse[5, 0, 9] = [1: -1]"));
    assert!(
        output.contains("rvDworkK3PSeries[5, 0, 5] = [1: 5, 5: -14976], mod 5=[5: 4], height=1")
    );
    assert!(output.contains("rvDworkK3Group[13, 2, 3] = [(1, 0): 1, (0, 1): 1, (1, 1): 8]"));
    assert!(output.contains("rvDworkK3Inverse[13, 2, 3] = [1: -1, 2: 8, 3: -64]"));
    let other = stdout(interpret(
        "rvDworkK3Log[13,-11,3]; rvDworkK3Group[13,-11,2]; rvDworkK3PSeries[13,-11,13];",
    ));
    assert!(other.contains("rvDworkK3Log[13, -11, 3] = [1: 1, 2: 22, 3: 1936/3]"));
    assert!(other.contains("(1, 1): -44"));
    assert!(other.contains("mod 13=[13: 7], height=1"));
}

#[test]
fn higher_hasse_witt_ratios_give_exact_frobenius_unit_residues() {
    let output = stdout(example("dwork-k3-formal-brauer"));
    let cases = [(5u64, 0, [4, 19, 19]), (13, 2, [7, 163, 1853])];
    for (p, psi, expected) in cases {
        for (index, value) in expected.into_iter().enumerate() {
            let precision = index + 1;
            let modulus = p.pow(precision as u32);
            let prefix = format!("rvDworkK3Unit[{p}, {psi}, {precision}] = mod {modulus},");
            let row = output
                .lines()
                .find(|line| line.starts_with(&prefix))
                .unwrap();
            assert!(row.ends_with(&format!("unit={value}")), "{row}");
            if p == 5 {
                assert_eq!((value * value + 6 * value + 25) % modulus, 0);
            }
        }
    }

    // The Witt Frobenius fixes all Z_p coefficients for a prime field.
    // The alternate lift gives the same certified root at the full precision.
    let other = stdout(interpret("rvDworkK3Unit[13,-11,3]; rvDworkK3Unit[5,5,3];"));
    assert!(other.lines().next().unwrap().ends_with("unit=1853"));
    assert!(other.lines().nth(1).unwrap().ends_with("unit=19"));

    // Independent projective point counts check the first Hasse–Witt digit.
    for (p, psi, count) in [(5, 0, 0i64), (7, 3, 56), (13, 2, 320)] {
        assert_eq!(points(p, psi), count);
        let output = stdout(interpret(&format!("rvDworkK3Unit[{p},{psi},1];")));
        assert!(output
            .trim_end()
            .ends_with(&format!("unit={}", (count - 1).rem_euclid(p))));
    }
}

fn points(p: i64, psi: i64) -> i64 {
    let mut count = 0;
    for pivot in 0..4 {
        for mut tail in 0..p.pow(3 - pivot) {
            let mut x = [0; 4];
            x[pivot as usize] = 1;
            for entry in &mut x[pivot as usize + 1..] {
                *entry = tail % p;
                tail /= p;
            }
            let value = x.iter().map(|&value| value.pow(4)).sum::<i64>()
                - 4 * psi * x.iter().product::<i64>();
            count += i64::from(value % p == 0);
        }
    }
    count
}

#[test]
fn height_distinguishes_certificates_from_truncated_zeroes() {
    let paths = stdout(example("dwork-k3-formal-brauer"));
    assert!(paths.contains("mod 3=[], height=infinite (supersingular Fermat quartic)"));
    assert!(paths.contains("mod 13=[], height lower bound=2"));
    assert!(paths.contains("Hasse–Witt=0, nonordinary, no Frobenius unit root"));

    let short = stdout(interpret(
        "rvDworkK3PSeries[5,0,4]; rvDworkK3PSeries[3,3,9];",
    ));
    assert!(short.contains("mod 5=[], height lower bound=1"));
    assert!(short.contains("mod 3=[], height=infinite (supersingular Fermat quartic)"));
}

#[test]
fn validates_smooth_reduction_and_checked_resources() {
    for (source, message) in [
        ("rvDworkK3Log[2,0,5];", "require an odd prime"),
        ("rvDworkK3Log[15,0,5];", "15 is not prime"),
        ("rvDworkK3Group[13,5,3];", "singular Dwork quartic"),
        ("rvDworkK3PSeries[13,12,3];", "singular Dwork quartic"),
        ("rvDworkK3Inverse[5,0,33];", "degree must be in 1..=32"),
        ("rvDworkK3Unit[5,0,0];", "precision must be positive"),
        ("rvDworkK3Unit[5,0,40];", "p^s must fit in u64"),
        (
            "rvDworkK3Log[13,170141183460469231731687303715884105727,32];",
            "exceeds the checked exact coefficient range",
        ),
    ] {
        let failure = diagnostic(interpret(source));
        assert!(failure.contains(message), "{failure}");
    }
}

#[test]
fn unramified_norm_uses_the_teichmuller_lift_and_the_twisted_denominator() {
    let output = stdout(example("dwork-k3-formal-brauer"));
    assert!(output.contains(
        "rvDworkK3ClosedUnit mod 81 = q=9, Teichmüller psi=[70, 70], Hasse–Witt=[0, 2], Cartier factor=[0, 59], unit=79"
    ));
    let output = stdout(interpret(
        "rvunram K[3,3]{1,0,1};
         rvDworkK3ClosedUnit(rvUnramElem[K]{1,-1});
         rvDworkK3ClosedUnit(rvUnramElem[K]{4,4});
         rvunram L[5,3]{2,0,1};
         rvDworkK3ClosedUnit(rvUnramElem[L]{});",
    ));
    let lines = output.lines().collect::<Vec<_>>();
    assert!(
        lines[0].ends_with("Cartier factor=[0, 22], unit=25"),
        "{output}"
    );
    assert!(
        lines[1].ends_with("Cartier factor=[0, 5], unit=25"),
        "{output}"
    );
    assert!(
        lines[2].ends_with("Cartier factor=[19, 0], unit=111"),
        "{output}"
    );

    // Independent projective counts of the same surface over F_9 and F_81.
    assert_eq!(79 % 9, (152 - 1) % 9);
    assert_eq!(79 * 79 % 81, (8024 - 1) % 81);

    let failure = diagnostic(interpret(
        "rvunram K[3,3]{1,0,1}; rvDworkK3ClosedUnit(rvUnramElem[K]{0,1});",
    ));
    assert!(failure.contains("singular Dwork quartic"));
}
