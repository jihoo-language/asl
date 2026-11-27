mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

fn value(line: &str) -> &str {
    line.split_once(" = ").unwrap().1
}

fn ratio(text: &str) -> (i128, i128) {
    let (n, d) = text.split_once('/').unwrap_or((text, "1"));
    (n.parse().unwrap(), d.parse().unwrap())
}

fn entries(text: &str) -> Vec<(i128, i128)> {
    text.split([',', '[', ']'])
        .filter_map(|s| {
            let s = s.trim();
            (!s.is_empty()).then(|| ratio(s))
        })
        .collect()
}

fn assert_characteristic(matrix: &str, trace: i128, norm: i128) {
    let values = entries(matrix);
    let [(a, a_d), (b, b_d), (c, c_d), (d, d_d)] = values[..] else {
        panic!("expected a 2x2 matrix, got {matrix}");
    };
    assert_eq!(a * d_d + d * a_d, trace * a_d * d_d, "{matrix}");
    assert_eq!(
        a * d * b_d * c_d - b * c * a_d * d_d,
        norm * a_d * b_d * c_d * d_d,
        "{matrix}"
    );
}

#[test]
fn dimensions_and_classical_weight_twelve_periods() {
    let mut program = String::new();
    let mut expected = String::new();
    for k in (2..=48).step_by(2) {
        let dim = (k / 12_i32 - i32::from(k % 12 == 2)).max(0);
        for sign in ["plus", "minus"] {
            writeln!(program, "rvPeriodDim[{k},{sign}];").unwrap();
            writeln!(expected, "rvPeriodDim[{k}, {sign}] = {dim}").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&program)), expected);

    let report = stdout(example("level-one-period-polynomials"));
    let lines: Vec<_> = report.lines().collect();
    assert_eq!(
        &lines[..3],
        [
            "rvPeriodDim[12, plus] = 1",
            "rvPeriodDim[12, minus] = 1",
            "rvPeriodCoeff[2] = -24",
        ]
    );
    assert_characteristic(value(lines[3]), 1080, -20468736);
    assert_characteristic(value(lines[4]), 1080, -20468736);
    assert_eq!(lines[5], "rvHeckeChar = [1, -1080, -20468736]");
}

#[test]
fn checked_literals_normalize_the_eisenstein_coboundary() {
    let source = "
        let p = rvPeriodPoly[12,plus](0,0,1,0,-3,0,3,0,-1,0,0);
        let q = rvPeriodPoly[12,plus](-3,0,(1,2),0,(-3,2),0,(3,2),0,(-1,2),0,3);
        let e = rvPeriodPoly[12,plus](-1,0,0,0,0,0,0,0,0,0,1);
        assert p == rvPeriodScale[2,1](q);
        assert add(q,neg(q)) == zero[RVPeriod(12,plus)];
        assert e == zero[RVPeriod(12,plus)];
        assert rvPeriodHecke(rvHeckeT[12,3],q) == rvPeriodScale[252,1](q);
        assert rvPeriodHecke(rvHeckeT[12,2],e) == e;
        assert rvPeriodPoly[4,plus](-1,0,1) == zero[RVPeriod(4,plus)];
        assert rvPeriodPoly[2,plus](0) == zero[RVPeriod(2,plus)];
        assert rvPeriodHecke(rvHeckeT[2,3],zero[RVPeriod(2,minus)]) == zero[RVPeriod(2,minus)];
        rvPeriodHeckeMatrix[2,minus](rvHeckeUnit[2]);
        rvPeriodCoeff[2](q);
        rvPeriodCoeff[10](q);
        rvPeriodCoeff[11](q);";
    assert_eq!(
        stdout(interpret(source)),
        "rvPeriodHeckeMatrix = []\nrvPeriodCoeff[2] = 1/2\nrvPeriodCoeff[10] = 0\nrvPeriodCoeff[11] = 0\n"
    );
}

#[test]
fn both_signs_match_the_cusp_eigenpacket_and_formal_composition() {
    let mut program = String::from("rvheckepacket E(rvHeckeT[24,2]);\n");
    for n in 1..=6 {
        writeln!(
            program,
            "rvPeriodHeckeMatrix[24,plus](rvHeckeT[24,{n}]);
             rvPeriodHeckeMatrix[24,minus](rvHeckeT[24,{n}]);
             rvHeckeEigenTrace(rvHeckeEigenCoeff[E,{n}]);
             rvHeckeEigenNorm(rvHeckeEigenCoeff[E,{n}]);",
        )
        .unwrap();
    }
    for sign in ["plus", "minus"] {
        writeln!(program,
            "let p{sign} = add(rvPeriodBasis[24,{sign},0],rvPeriodScale[3,2](rvPeriodBasis[24,{sign},1]));
             let a{sign} = rvHeckeT[24,2]; let b{sign} = rvHeckeT[24,3];
             assert rvPeriodHecke(rvHeckeMul(a{sign},b{sign}),p{sign}) ==
                rvPeriodHecke(a{sign},rvPeriodHecke(b{sign},p{sign}));
             assert rvPeriodHecke(rvHeckeMul(a{sign},a{sign}),p{sign}) ==
                rvPeriodHecke(a{sign},rvPeriodHecke(a{sign},p{sign}));
             assert rvPeriodHecke(rvHeckeUnit[24],p{sign}) == p{sign};",
        ).unwrap();
    }
    let report = stdout(interpret(&program));
    let lines: Vec<_> = report.lines().collect();
    assert_eq!(lines.len(), 24);
    for row in lines.chunks_exact(4) {
        let trace = value(row[2]).parse().unwrap();
        let norm = value(row[3]).parse().unwrap();
        assert_characteristic(value(row[0]), trace, norm);
        assert_characteristic(value(row[1]), trace, norm);
    }
}

fn choose(n: usize, k: usize) -> i128 {
    (1..=k).fold(1, |a, j| a * (n - j + 1) as i128 / j as i128)
}

fn direct_substitution(p: &[i128], [a, b, c, d]: [i128; 4]) -> Vec<i128> {
    let w = p.len() - 1;
    let mut output = vec![0; w + 1];
    for (j, &p_j) in p.iter().enumerate().filter(|&(_, &p_j)| p_j != 0) {
        for i in 0..=j {
            for r in 0..=w - j {
                // Expand independent binomial powers of both forms directly.
                output[i + r] += p_j
                    * choose(j, i)
                    * a.pow(i as u32)
                    * b.pow((j - i) as u32)
                    * choose(w - j, r)
                    * c.pow(r as u32)
                    * d.pow((w - j - r) as u32);
            }
        }
    }
    output
}

fn direct_hecke(p: &[i128], n: i128) -> Vec<i128> {
    let mut output = vec![0; p.len()];
    // Independent four-loop determinant sieve, restricted to tiny indices.
    for a in 1..=n {
        for d in 1..=n {
            for c in 0..a {
                for b in 0..d {
                    if a * d - b * c == n {
                        for (r, term) in output.iter_mut().zip(direct_substitution(p, [a, b, c, d]))
                        {
                            *r += term;
                        }
                    }
                }
            }
        }
    }
    output
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn read_polynomial(line: &str) -> Vec<(i128, i128)> {
    entries(line.split_once('{').unwrap().1.strip_suffix('}').unwrap())
}

#[test]
fn mixed_operator_and_rational_vectors_match_direct_manin_and_merel_checks() {
    for sign in ["plus", "minus"] {
        let bindings = format!(
            "let p = add(rvPeriodBasis[24,{sign},0],rvPeriodScale[3,2](rvPeriodBasis[24,{sign},1]));
             let a = add(rvHeckeT[24,2],add(rvHeckeScale[2](rvHeckeT[24,3]),rvHeckeScale[-3](rvHeckeT[24,4])));",
        );
        let report = stdout(interpret(&format!(
            "{bindings} show p; show rvPeriodHecke(a,p);"
        )));
        let mut lines = report.lines().map(read_polynomial);
        let source = lines.next().unwrap();
        let actual = lines.next().unwrap();
        let d = source.iter().fold(1, |d, &(_, b)| d / gcd(d, b) * b);
        let p: Vec<_> = source.iter().map(|&(a, b)| a * (d / b)).collect();
        assert!(p.iter().filter(|&&x| x != 0).count() > 2);
        for (&a, b) in p.iter().zip(direct_substitution(&p, [0, -1, 1, 0])) {
            assert_eq!(a + b, 0);
        }
        let u = direct_substitution(&p, [0, -1, 1, 1]);
        let u2 = direct_substitution(&p, [-1, -1, 1, 0]);
        assert!(p.iter().zip(u).zip(u2).all(|((&x, y), z)| x + y + z == 0));
        let mut expected = direct_hecke(&p, 2);
        for ((c, b), a) in expected
            .iter_mut()
            .zip(direct_hecke(&p, 3))
            .zip(direct_hecke(&p, 4))
        {
            *c += 2 * b - 3 * a;
        }
        if sign == "plus" {
            let w = expected.len() - 1;
            expected[0] += expected[w];
            expected[w] = 0;
        }
        for (j, ((a, b), c)) in actual.into_iter().zip(expected).enumerate() {
            assert_eq!(a * d, b * c, "{sign} X^{j}");
        }
    }
}

#[test]
fn malformed_manin_polynomials_and_weight_mismatches_are_rejected() {
    for (source, expected) in [
        ("rvPeriodDim[3,plus];", "even level-one weight"),
        ("rvPeriodDim[12,other];", "sign must be plus or minus"),
        ("show rvPeriodBasis[2,plus,0];", "basis index"),
        ("show rvPeriodPoly[12,plus](0);", "exactly k−1"),
        ("show rvPeriodPoly[4,plus](0,1,0);", "wrong X-parity"),
        ("show rvPeriodPoly[4,plus](1,0,0);", "Manin relation 1+S"),
        (
            "show rvPeriodPoly[12,plus](0,0,1,0,0,0,0,0,-1,0,0);",
            "Manin relation 1+U+U²",
        ),
        (
            "show rvPeriodScale[1,0](zero[RVPeriod(2,plus)]);",
            "denominator must be positive",
        ),
        (
            "show rvPeriodHecke(rvHeckeT[24,2],rvPeriodBasis[12,plus,0]);",
            "RVHecke(12)",
        ),
        (
            "show add(zero[RVPeriod(12,plus)],zero[RVPeriod(12,minus)]);",
            "RVPeriod(12,plus)",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}
