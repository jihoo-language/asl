mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn marking() -> &'static str {
    include_str!("../examples/hyperelliptic-kummer-weil-theta.asl")
        .split_once("// Weil representation and theta coefficients below.")
        .unwrap()
        .0
}

#[test]
fn low_norm_example_and_fourier_fraction() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-weil-theta")),
        concat!(
            "rvKummerWeilCoeff = 1/8\n",
            "rvKummerWeilCoeff = -1/8\n",
            "rvKummerWeilCoeff = 1\n",
            "rvKummerWeilCoeff = 32\n",
            "rvKummerWeilCoeff = 64\n",
            "rvKummerWeilCoeff = 1024\n",
            "rvKummerWeilCoeff = 8160\n",
        )
    );
}

#[test]
fn generators_satisfy_all_relations_and_character_traces_on_the_entire_basis() {
    let mut script = marking().to_owned();
    for i in 0..64 {
        // Each distinct label is a basis vector and its three transformed images.
        writeln!(script, "let d{i}=rvKummerDisc[B,{i}]; let a{i}=rvKummerWeilBasis(d{i}); let s{i}=rvKummerWeilS(a{i}); let t{i}=rvKummerWeilT(a{i}); let st{i}=rvKummerWeilS(t{i});").unwrap();
        writeln!(
            script,
            "assert rvKummerWeilS(s{i})==a{i}; assert rvKummerWeilT(t{i})==a{i};"
        )
        .unwrap();
        let mut cube = format!("st{i}");
        for _ in 0..2 {
            cube = format!("rvKummerWeilS(rvKummerWeilT({cube}))");
        }
        writeln!(script, "assert {cube}==a{i};").unwrap();
        writeln!(
            script,
            "assert rvKummerWeilFrob[C,1](a{i})==rvKummerWeilBasis(rvKummerDiscFrob[C,1](d{i}));"
        )
        .unwrap();
        for (op, prefix) in [("S", "s"), ("T", "t")] {
            writeln!(script, "assert rvKummerWeilFrob[C,1]({prefix}{i})==rvKummerWeil{op}(rvKummerWeilFrob[C,1](a{i}));").unwrap();
        }
        for prefix in ["s", "t", "st"] {
            writeln!(script, "rvKummerWeilCoeff({prefix}{i},d{i});").unwrap();
        }
    }
    let output = stdout(interpret(&script));
    let coefficients: Vec<_> = output
        .lines()
        .map(|line| {
            let value = line.strip_prefix("rvKummerWeilCoeff = ").unwrap();
            if let Some((n, d)) = value.split_once('/') {
                let divisor: i32 = d.parse().unwrap();
                assert_eq!(8 % divisor, 0);
                n.parse::<i32>().unwrap() * (8 / divisor)
            } else {
                value.parse::<i32>().unwrap() * 8
            }
        })
        .collect();
    assert_eq!(coefficients.len(), 3 * 64);
    let trace = (0..3)
        .map(|i| coefficients.iter().skip(i).step_by(3).sum::<i32>() / 8)
        .collect::<Vec<_>>();
    assert_eq!(trace, [8, 8, 1]);
    assert_eq!(
        coefficients
            .iter()
            .skip(1)
            .step_by(3)
            .filter(|&&a| a == 8)
            .count(),
        36
    );
    assert_eq!(
        coefficients
            .iter()
            .skip(1)
            .step_by(3)
            .filter(|&&a| a == -8)
            .count(),
        28
    );
    let standard = (64 - trace[2]) / 3;
    let trivial = (64 - 2 * standard + trace[0]) / 2;
    assert_eq!([trivial, trivial - trace[0], standard], [15, 7, 21]);
}

fn coset_weights(mask: u8) -> [u32; 32] {
    let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    std::array::from_fn(|a| {
        (0u8..16)
            .filter(|&v| {
                let quad = pairs
                    .iter()
                    .enumerate()
                    .filter(|&(i, &(l, r))| {
                        mask & (1 << i) != 0 && v & (1 << l) != 0 && v & (1 << r) != 0
                    })
                    .count()
                    % 2;
                let linear = ((a >> 1) & usize::from(v)).count_ones() as usize % 2;
                (quad + linear + (a & 1)) % 2 == 1
            })
            .count() as u32
    })
}

#[test]
fn all_cosets_match_direct_low_norm_patterns_and_theta_equivariance() {
    let mut script = marking().to_owned();
    let mut expected = String::new();
    for norm in [0, 2, 3, 4] {
        writeln!(script, "let m{norm}=rvKummerThetaShell[B,{norm}];").unwrap();
        writeln!(script, "assert rvKummerWeilFrob[C,1](m{norm})==m{norm};").unwrap();
        let sign = if norm % 2 == 0 { "" } else { "neg" };
        let target = if sign.is_empty() {
            format!("m{norm}")
        } else {
            format!("{sign}(m{norm})")
        };
        writeln!(script, "assert rvKummerWeilT(m{norm})=={target};").unwrap();
        for mask in 0..64 {
            // Patterns of squares 0, 1 and 4 suffice here. Each nonzero
            // coordinate independently chooses a sign.
            let count: u64 = coset_weights(mask)
                .iter()
                .map(|&w| match (norm, w) {
                    (0, 0) => 1,
                    (2, 0) => 32,
                    (2, 4) => 16,
                    (3, 6) => 64,
                    (4, 0) => 480,
                    (4, 4) => 384,
                    (4, 8) => 256,
                    _ => 0,
                })
                .sum();
            writeln!(script, "rvKummerWeilCoeff(m{norm},rvKummerDisc[B,{mask}]);").unwrap();
            writeln!(expected, "rvKummerWeilCoeff = {count}").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&script)), expected);
}

#[test]
fn wide_rational_normalization_preserves_fourier_cancellation() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
        let e=rvKummerWeilBasis(rvKummerDisc[B,0]);
        let f=rvKummerWeilS(e);
        let tiny=rvKummerWeilScale[8,340282366920938463463374607431768211455](f);
        assert rvKummerWeilS(tiny)==rvKummerWeilScale[8,340282366920938463463374607431768211455](e);
        let large=rvKummerWeilScale[85070591730234615865843651857942052864,1](f);
        assert rvKummerWeilS(large)==rvKummerWeilScale[85070591730234615865843651857942052864,1](e);
        assert add(f,neg(f))==zero[RVKummerWeil(B)];
        assert add(rvKummerWeilScale[1,3](e),rvKummerWeilScale[1,6](e))==rvKummerWeilScale[1,2](e);
        ",
            marking()
        ))),
        ""
    );
}

#[test]
fn rejects_unrelated_markings_bad_scalars_and_excessive_shells() {
    let setup = format!("{} let e=rvKummerWeilBasis(rvKummerDisc[B,0]);", marking());
    for (source, expected) in [
        ("show rvKummerWeilS(rvKummerDisc[B,0]);", "requires an RVKummerWeil vector"),
        ("show rvKummerWeilScale[1,0](e);", "denominator must be positive"),
        ("show rvKummerThetaShell[B,1000000];", "shell exceeds the polynomial work or storage limit"),
        ("show rvKummerThetaShell[Missing,2];", "unknown Jacobian torsion basis"),
        ("rvjactorsionbasis Twin[L,2](neg(ow)){r0,e2,r1,f2}; rvKummerWeilCoeff(e,rvKummerDisc[Twin,0]);", "expected type RVKummerDisc(B)"),
        ("rvhyperelliptic D[F]{z,o,z,z,o,o}; show rvKummerWeilFrob[D,1](e);", "not the declared base extension"),
    ] {
        let problem = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
