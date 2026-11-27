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

#[test]
fn haberland_normalization_is_alternating_rational_and_parabolic() {
    let source = "
        let p = rvPeriodPoly[12,plus](0,0,1,0,-3,0,3,0,-1,0,0);
        let q = rvPeriodPoly[12,minus](0,4,0,-25,0,42,0,-25,0,4,0);
        let e = rvPeriodPoly[12,plus](-1,0,0,0,0,0,0,0,0,0,1);
        rvPeriodPair(p,q); rvPeriodPair(q,p);
        rvPeriodPair(rvPeriodScale[1,7](p),rvPeriodScale[1,5](q));
        rvPeriodPair(p,p); rvPeriodPair(q,q); rvPeriodPair(e,q); rvPeriodPair(q,e);
        rvPeriodPairing[2]; rvPeriodPairing[12];";
    let report = stdout(interpret(source));
    let lines: Vec<_> = report.lines().collect();
    assert_eq!(
        lines[..8],
        [
            "rvPeriodPair = 6",
            "rvPeriodPair = -6",
            "rvPeriodPair = 6/35",
            "rvPeriodPair = 0",
            "rvPeriodPair = 0",
            "rvPeriodPair = 0",
            "rvPeriodPair = 0",
            "rvPeriodPairing[2] = denominator=1, matrix=[], smith=[], elementary=[]",
        ]
    );
    assert!(lines[8].starts_with("rvPeriodPairing[12] = denominator=1,"));
    assert!(lines[8].ends_with("smith=[6], elementary=[6]"));
    let demo = stdout(example("haberland-period-pairings"));
    let lines: Vec<_> = demo.lines().collect();
    assert_eq!(&lines[..2], &["rvPeriodPair = 6", "rvPeriodPair = -6"]);
    assert_eq!(lines[2], lines[3]);
    assert!(lines[4].starts_with("rvPeriodPairing[24] = denominator=3553,"));
}

fn determinant(matrix: &[Vec<i128>]) -> i128 {
    if matrix.is_empty() {
        return 1;
    }
    matrix[0]
        .iter()
        .enumerate()
        .map(|(j, coefficient)| {
            let minor: Vec<_> = matrix[1..]
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter(|&(i, _)| i != j)
                        .map(|(_, &x)| x)
                        .collect()
                })
                .collect();
            if j.is_multiple_of(2) {
                coefficient * determinant(&minor)
            } else {
                -coefficient * determinant(&minor)
            }
        })
        .sum()
}

fn saturation_index(basis: &[Vec<i128>]) -> i128 {
    let width = basis[0].len();
    let columns: Vec<_> = (1..=(width - 1) / 2)
        .filter(|&j| basis.iter().any(|row| row[j] != 0))
        .collect();
    let mut index = 0;
    for mask in 0..(1u32 << columns.len()) {
        if mask.count_ones() as usize == basis.len() {
            let minor: Vec<_> = basis
                .iter()
                .map(|row| {
                    columns
                        .iter()
                        .enumerate()
                        .filter(|&(i, _)| mask & (1u32 << i) != 0)
                        .map(|(_, &j)| row[j])
                        .collect()
                })
                .collect();
            index = gcd(index, determinant(&minor).abs());
        }
    }
    index
}

fn integral_rows(report: &str) -> Vec<Vec<i128>> {
    report
        .lines()
        .map(|line| {
            read_polynomial(line)
                .into_iter()
                .map(|(n, d)| {
                    assert_eq!(d, 1);
                    n
                })
                .collect()
        })
        .collect()
}

fn independently_pair(p: &[i128], q: &[i128]) -> (i128, i128) {
    let w = p.len() - 1;
    let left = direct_substitution(p, [1, 1, 0, 1]);
    let right = direct_substitution(p, [1, -1, 0, 1]);
    let scale = (0..=w)
        .map(|i| choose(w, i))
        .fold(1, |d, e| d / gcd(d, e) * e);
    let numerator = (0..=w)
        .map(|i| {
            let sign = if i.is_multiple_of(2) { 1 } else { -1 };
            sign * (left[i] - right[i]) * q[w - i] * (scale / choose(w, i))
        })
        .sum::<i128>();
    let common = gcd(numerator.abs(), scale);
    (numerator / common, scale / common)
}

#[test]
fn integral_period_basis_is_saturated_and_reports_the_actual_haberland_matrix() {
    for (weight, dimension, denominator, smith, elementary) in [
        (24, 2, 3553, "[6, 639540]", "[6/3553, 180]"),
        (
            36,
            3,
            7733198,
            "[3, 53130, 4871914740]",
            "[3/7733198, 105/15283, 630]",
        ),
    ] {
        let mut program = String::new();
        for sign in ["plus", "minus"] {
            for i in 0..dimension {
                writeln!(program, "show rvPeriodIntegralBasis[{weight},{sign},{i}];").unwrap();
            }
        }
        let basis = stdout(interpret(&program));
        let vectors = integral_rows(&basis);
        let (plus, minus) = vectors.split_at(dimension);
        assert_eq!(saturation_index(plus), 1);
        assert_eq!(saturation_index(minus), 1);
        for p in &vectors {
            assert_eq!([p[0], *p.last().unwrap()], [0, 0]);
            let s = direct_substitution(p, [0, -1, 1, 0]);
            let u = direct_substitution(p, [0, -1, 1, 1]);
            let u2 = direct_substitution(p, [-1, -1, 1, 0]);
            assert!(p.iter().zip(s).all(|(&x, y)| x + y == 0));
            assert!(p.iter().zip(u).zip(u2).all(|((&x, y), z)| x + y + z == 0));
        }
        let expected: Vec<_> = plus
            .iter()
            .flat_map(|p| minus.iter().map(|q| independently_pair(p, q)))
            .collect();
        let actual_denominator = expected.iter().fold(1, |d, &(_, e)| d / gcd(d, e) * e);
        assert_eq!(actual_denominator, denominator);
        let report = stdout(interpret(&format!("rvPeriodPairing[{weight}];")));
        let (_, matrix_and_tail) = report.split_once(", matrix=").unwrap();
        let (matrix, tail) = matrix_and_tail.split_once(", smith=").unwrap();
        assert_eq!(entries(matrix), expected);
        assert_eq!(tail.trim(), format!("{smith}, elementary={elementary}"));
    }

    // Clearing each rational basis vector separately really gives a
    // proper sublattice in weight 24, even when each vector is primitive.
    for sign in ["plus", "minus"] {
        let report = stdout(interpret(&format!(
            "show rvPeriodBasis[24,{sign},0]; show rvPeriodBasis[24,{sign},1];"
        )));
        let separately_cleared: Vec<_> = report
            .lines()
            .map(|line| {
                let terms = read_polynomial(line);
                let scale = terms.iter().fold(1, |d, &(_, e)| d / gcd(d, e) * e);
                terms.into_iter().map(|(n, d)| n * (scale / d)).collect()
            })
            .collect();
        assert_eq!(saturation_index(&separately_cleared), 323);
    }
}

#[test]
fn haberland_pairing_makes_mixed_hecke_operators_self_adjoint() {
    for (weight, dimension) in [(12, 1), (24, 2), (36, 3)] {
        let mut program = format!(
            "let a = add(rvHeckeT[{weight},2],add(rvHeckeScale[2](rvHeckeT[{weight},3]),rvHeckeScale[-3](rvHeckeT[{weight},4])));"
        );
        for i in 0..dimension {
            writeln!(
                program,
                "let p{i} = rvPeriodIntegralBasis[{weight},plus,{i}];
                 let q{i} = rvPeriodIntegralBasis[{weight},minus,{i}];"
            )
            .unwrap();
        }
        for i in 0..dimension {
            for j in 0..dimension {
                writeln!(
                    program,
                    "rvPeriodPair(rvPeriodHecke(a,p{i}),q{j});
                     rvPeriodPair(p{i},rvPeriodHecke(a,q{j}));"
                )
                .unwrap();
            }
        }
        let report = stdout(interpret(&program));
        let lines: Vec<_> = report.lines().collect();
        assert_eq!(lines.len(), 2 * dimension * dimension);
        for row in lines.chunks_exact(2) {
            assert_eq!(row[0], row[1], "weight {weight}");
        }
        assert!(lines.iter().any(|&line| line != "rvPeriodPair = 0"));
    }
}

#[test]
fn haberland_types_and_integral_basis_indices_are_checked() {
    for (source, expected) in [
        (
            "rvPeriodPair(zero[RVPeriod(12,plus)],zero[RVPeriod(24,minus)]);",
            "requires the same weight",
        ),
        (
            "rvPeriodPair(One.pt,One.pt);",
            "an RVPeriod(k,sign) polynomial",
        ),
        (
            "show rvPeriodIntegralBasis[2,minus,0];",
            "integral period-polynomial basis index",
        ),
        ("rvPeriodPairing[13];", "even level-one weight"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}
