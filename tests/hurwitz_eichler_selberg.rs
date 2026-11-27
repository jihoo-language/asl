mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

fn integer(line: &str) -> i128 {
    line.rsplit_once(" = ").unwrap().1.parse().unwrap()
}

fn scaled_class_number(line: &str) -> i128 {
    let value = line.rsplit_once(" = ").unwrap().1;
    let (num, den) = value.split_once('/').unwrap_or((value, "1"));
    let numerator: i128 = num.parse().unwrap();
    let denominator: i128 = den.parse().unwrap();
    assert_eq!(12 % denominator, 0);
    numerator * (12 / denominator)
}

fn brute_reduced_forms(n: i128) -> i128 {
    if n < 0 {
        return 0;
    }
    if n == 0 {
        return -1;
    }
    // Count all forms directly, imposing a positive middle coefficient
    // on the boundary. This includes imprimitive forms without any
    // square-divisor decomposition.
    let mut count = 0;
    for a in 1..=n {
        for b in -a..=a {
            if (b * b + n) % (4 * a) != 0 {
                continue;
            }
            let c = (b * b + n) / (4 * a);
            if a > c || (b < 0 && (b == -a || a == c)) {
                continue;
            }
            count += match (a == c, b) {
                (true, 0) => 6,
                (true, b) if b == a => 4,
                _ => 12,
            };
        }
    }
    count
}

#[test]
fn hurwitz_numbers_include_imprimitive_forms_and_both_stabilizers() {
    let mut source = String::new();
    for n in -4..=150 {
        writeln!(source, "rvHurwitzClass[{n}];").unwrap();
    }
    let report = stdout(interpret(&source));
    for (n, line) in (-4..=150).zip(report.lines()) {
        assert_eq!(scaled_class_number(line), brute_reduced_forms(n), "H({n})");
    }
    assert_eq!(report.lines().count(), 155);
    assert_eq!(
        stdout(interpret(
            "rvHurwitzClass[3]; rvHurwitzClass[4]; rvHurwitzClass[36];"
        )),
        "rvHurwitzClass[3] = 1/3\nrvHurwitzClass[4] = 1/2\nrvHurwitzClass[36] = 5/2\n"
    );
}

#[test]
fn jacobi_coefficients_depend_on_the_exact_discriminant() {
    // k² is in the signed literal range, while 4k² does not fit in u128.
    let k = 13_000_000_000_000_000_000_i128;
    let even = k * k + 1;
    let odd = k * k - k + 1;
    let r = 2 * k;
    let source = format!(
        "rvHurwitzJacobi[1,0]; rvHurwitzJacobi[{even},{r}]; rvHurwitzJacobi[{even},-{r}];\
         rvHurwitzJacobi[1,1]; rvHurwitzJacobi[{odd},{}]; rvHurwitzJacobi[{},{r}];\
         rvHurwitzJacobi[0,{r}]; rvHurwitzJacobi[-1,0];",
        r - 1,
        k * k
    );
    let report = stdout(interpret(&source));
    assert_eq!(
        report.lines().map(scaled_class_number).collect::<Vec<_>>(),
        [6, 6, 6, 4, 4, -1, 0, 0]
    );
}

fn delta(bound: usize) -> Vec<i128> {
    let mut q = vec![0; bound + 1];
    q[1] = 1;
    // Independent multiplication of Δ=q∏(1-qᵐ)²⁴, one linear
    // factor at a time. The implementation uses class numbers instead.
    for m in 1..bound {
        for _ in 0..24 {
            for j in (m..=bound).rev() {
                q[j] -= q[j - m];
            }
        }
    }
    q
}

fn e4(bound: usize) -> Vec<i128> {
    let mut q = vec![1];
    for n in 1..=bound {
        let sigma: i128 = (1..=n)
            .filter(|d| n.is_multiple_of(*d))
            .map(|d| (d as i128).pow(3))
            .sum();
        q.push(240 * sigma);
    }
    q
}

fn convolve(left: &[i128], right: &[i128], bound: usize) -> Vec<i128> {
    (0..=bound)
        .map(|n| (0..=n).map(|j| left[j] * right[n - j]).sum())
        .collect()
}

#[test]
fn elliptic_and_parabolic_hecke_terms_give_tau_and_weight_two_zero() {
    let max = 16;
    let tau = delta(max);
    let mut source = String::new();
    for n in 1..=max {
        writeln!(source, "rvHeckeTrace[2,{n}]; rvHeckeTrace[12,{n}];").unwrap();
    }
    let report = stdout(interpret(&source));
    for ((n, weight), line) in (1..=max)
        .flat_map(|n| [(n, 2), (n, 12)])
        .zip(report.lines())
    {
        assert_eq!(
            integer(line),
            if weight == 2 { 0 } else { tau[n] },
            "T_{n} on S_{weight}"
        );
    }
    assert_eq!(report.lines().count(), 2 * max);
}

#[test]
fn weight_twenty_four_traces_match_a_q_expansion_matrix() {
    let bound = 24;
    let delta = delta(bound);
    let e4 = e4(bound);
    let e4_cubed = convolve(&convolve(&e4, &e4, bound), &e4, bound);
    let f = convolve(&delta, &e4_cubed, bound);
    let g = convolve(&delta, &delta, bound);
    assert_eq!([f[1], f[2], g[1], g[2]], [1, 696, 0, 1]);
    let mut source = String::new();
    for n in 1..=bound / 2 {
        writeln!(source, "rvHeckeTrace[24,{n}];").unwrap();
    }
    let report = stdout(interpret(&source));
    for (n, line) in (1..=bound / 2).zip(report.lines()) {
        // In the basis f=ΔE₄³, g=Δ², the first diagonal entry is f_n.
        // The second is (T_n g)_2 − f_2 (T_n g)_1.
        let twice = if n.is_multiple_of(2) {
            (1_i128 << 23) * g[n / 2]
        } else {
            0
        };
        let trace = f[n] + g[2 * n] + twice - f[2] * g[n];
        assert_eq!(integer(line), trace, "T_{n} on S_24");
    }
    assert_eq!(report.lines().count(), bound / 2);
}

type Fourier = BTreeMap<(usize, i32), i128>;

fn jacobi_e4_times_phi_minus_two(bound: usize) -> Fourier {
    let mut series = Fourier::from([((0, -1), 1), ((0, 0), -2), ((0, 1), 1)]);
    for j in 1..=bound {
        for r in [-1, 1] {
            let factors = [(0, 0, 1), (j, r, -2), (2 * j, 2 * r, 1)];
            series = multiply_jacobi(&series, &factors, bound);
        }
        let factors: Vec<_> = (0..=bound / j)
            .map(|k| {
                let k_i = k as i128;
                (j * k, 0, (k_i + 1) * (k_i + 2) * (k_i + 3) / 6)
            })
            .collect();
        series = multiply_jacobi(&series, &factors, bound);
    }
    let eisenstein: Vec<_> = e4(bound)
        .into_iter()
        .enumerate()
        .map(|(n, c)| (n, 0, c))
        .collect();
    multiply_jacobi(&series, &eisenstein, bound)
}

fn multiply_jacobi(series: &Fourier, factor: &[(usize, i32, i128)], bound: usize) -> Fourier {
    let mut result = Fourier::new();
    for (&(n, r), &a) in series {
        for &(d, s, b) in factor {
            if n + d <= bound {
                *result.entry((n + d, r + s)).or_default() += a * b;
            }
        }
    }
    result
}

#[test]
fn zagiers_numbers_reconstruct_the_normalized_igusa_mock_jacobi_form() {
    let bound = 5;
    let modular = jacobi_e4_times_phi_minus_two(bound);
    let indices: Vec<_> = (0..=bound)
        .flat_map(|n| (-5..=5).map(move |r| (n, r)))
        .collect();
    let mut source = String::new();
    for &(n, r) in &indices {
        writeln!(
            source,
            "rvHurwitzJacobi[{n},{r}]; rvIgusaMockWeight2[1,{n},{r}];"
        )
        .unwrap();
    }
    let report = stdout(interpret(&source));
    let rows: Vec<_> = report.lines().collect();
    assert_eq!(rows.len(), 2 * indices.len());
    for ((n, r), row) in indices.into_iter().zip(rows.chunks_exact(2)) {
        let h = scaled_class_number(row[0]);
        let c = modular.get(&(n, r)).copied().unwrap_or_default();
        // DMZ equations (128) and (249): Δψ₁ᶠ=3E₄φ₋₂,₁−648ℋ.
        assert_eq!(integer(row[1]), 3 * c - 54 * h, "q^{n} y^{r}");
    }
}

#[test]
fn huge_identity_weights_and_range_diagnostics_stay_exact() {
    let k = u128::MAX - 1;
    let expected = k / 12 - u128::from(k % 12 == 2);
    let n = i128::MAX - 2;
    let lo = i128::MIN;
    assert_eq!(stdout(interpret(&format!("rvHeckeTrace[{k},1]; rvHurwitzClass[{n}]; rvHurwitzClass[{lo}];"))),
        format!("rvHeckeTrace[{k}, 1] = {expected}\nrvHurwitzClass[{n}] = 0\nrvHurwitzClass[{lo}] = 0\n"));
    for (source, expected) in [
        ("rvHeckeTrace[13,1];", "even weight at least 2"),
        ("rvHeckeTrace[0,1];", "even weight at least 2"),
        ("rvHeckeTrace[12,0];", "Hecke index n must be positive"),
        (
            "rvHeckeTrace[12,340282366920938463463374607431768211455];",
            "Hurwitz class number enumeration exceeds",
        ),
        (
            "rvHurwitzClass[170141183460469231731687303715884105727];",
            "Hurwitz class number enumeration exceeds",
        ),
        (
            "rvHeckeTrace[4096,2];",
            "exact intermediate exceeds 512 bits",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}

#[test]
fn the_example_shows_normalized_rationals_and_classical_traces() {
    assert_eq!(
        stdout(example("hurwitz-eichler-selberg")),
        concat!(
            "rvHurwitzClass[0] = -1/12\n",
            "rvHurwitzClass[12] = 4/3\n",
            "rvHurwitzClass[16] = 3/2\n",
            "rvHurwitzJacobi[1, 1] = 1/3\n",
            "rvHurwitzJacobi[3, 3] = 1/3\n",
            "rvHeckeTrace[2, 9] = 0\n",
            "rvHeckeTrace[12, 5] = 4830\n",
            "rvHeckeTrace[24, 1] = 2\n",
            "rvHeckeTrace[24, 2] = 1080\n",
        )
    );
}

#[test]
fn formal_hecke_operators_use_graded_integer_arithmetic() {
    let report = stdout(interpret(
        "let a = rvHeckeT[12,2];
         let b = rvHeckeT[12,3];
         let one = rvHeckeUnit[12];
         let z = zero[RVHecke(12)];
         assert one == rvHeckeT[12,1];
         assert rvHeckeMul(a,one) == a;
         assert rvHeckeMul(a,z) == z;
         assert rvHeckeMul(a,b) == rvHeckeT[12,6];
         assert rvHeckeMul(a,b) == rvHeckeMul(b,a);
         assert rvHeckeMul(a,a) == add(rvHeckeT[12,4],rvHeckeScale[2048](one));
         assert rvHeckeMul(rvHeckeT[12,4],rvHeckeT[12,6]) ==
            add(rvHeckeT[12,24],rvHeckeScale[2048](rvHeckeT[12,6]));
         assert rvHeckeMul(a,add(a,neg(b))) ==
            add(rvHeckeMul(a,a),neg(rvHeckeMul(a,b)));
         assert add(a,neg(a)) == z;
         assert rvHeckeScale[0](b) == z;
         assert rvHeckeSameAction(a,rvHeckeScale[-24](one)) == S.hi;
         assert rvHeckeSameAction(a,z) == S.lo;
         rvHeckeActionTrace(a);
         rvHeckeChar(a);
         rvHeckeCharCoeff[1](a);
         rvHeckeCharCoeff[2](a);",
    ));
    assert_eq!(
        report,
        "rvHeckeActionTrace = -24\nrvHeckeChar = [1, 24]\nrvHeckeCharCoeff[1] = 24\nrvHeckeCharCoeff[2] = 0\n"
    );
}

fn hecke24_matrix(f: &[i128], g: &[i128], index: usize) -> [[i128; 2]; 2] {
    let column = |form: &[i128]| {
        let first = hecke_q_coefficient(form, 24, index, 1);
        [
            first,
            hecke_q_coefficient(form, 24, index, 2) - f[2] * first,
        ]
    };
    [column(f), column(g)]
}

fn hecke_q_coefficient(form: &[i128], weight: u32, index: usize, j: usize) -> i128 {
    (1..=j)
        .filter(|&d| j.is_multiple_of(d) && index.is_multiple_of(d))
        .map(|d| (d as i128).pow(weight - 1) * form[(j * index) / (d * d)])
        .sum()
}

fn characteristic_2x2([[a, c], [b, d]]: [[i128; 2]; 2]) -> [i128; 3] {
    [1, -a - d, a * d - b * c]
}

#[test]
fn formal_hecke_characters_match_independent_q_expansion_action() {
    let bound = 12;
    let delta = delta(bound);
    let e4 = e4(bound);
    let e4_cubed = convolve(&convolve(&e4, &e4, bound), &e4, bound);
    let f = convolve(&delta, &e4_cubed, bound);
    let g = convolve(&delta, &delta, bound);
    let t2 = hecke24_matrix(&f, &g, 2);
    let t6 = hecke24_matrix(&f, &g, 6);
    assert_eq!(characteristic_2x2(t2), [1, -1080, -20468736]);
    let mixed = std::array::from_fn(|i| {
        std::array::from_fn(|j| 3 * i128::from(i == j) - 2 * t2[i][j] + t6[i][j])
    });
    let formal_relation = "let t = rvHeckeT[24,2];
        let one = rvHeckeUnit[24];
        let relation = add(add(rvHeckeMul(t,t),rvHeckeScale[-1080](t)),rvHeckeScale[-20468736](one));
        assert relation == zero[RVHecke(24)];";
    assert!(diagnostic(interpret(formal_relation)).contains("assertion failed"));

    let [one, trace, determinant] = characteristic_2x2(mixed);
    assert_eq!(
        stdout(example("graded-hecke-operators")),
        format!(
            "rvHeckeActionTrace = 1080\nrvHeckeChar = [1, -1080, -20468736]\n\
             RVHecke(24){{T_1: -12080128, T_2: -1080, T_4: 1}}\n\
             rvHeckeChar = [{one}, {trace}, {determinant}]\n"
        )
    );
}

#[test]
fn weight_thirty_six_characteristic_and_eigenpacket_match_cusp_form_action() {
    let bound = 9;
    let delta = delta(bound);
    let e4 = e4(bound);
    let e4_cubed = convolve(&convolve(&e4, &e4, bound), &e4, bound);
    let delta_squared = convolve(&delta, &delta, bound);
    let f = convolve(&delta, &convolve(&e4_cubed, &e4_cubed, bound), bound);
    let g = convolve(&delta_squared, &e4_cubed, bound);
    let h = convolve(&delta_squared, &delta, bound);
    assert_eq!([f[1], g[1], g[2], h[1], h[2], h[3]], [1, 0, 1, 0, 0, 1]);

    let hecke_matrix = |n| {
        [&f, &g, &h].map(|form| {
            let [u, v, w] = [1, 2, 3].map(|j| hecke_q_coefficient(form, 36, n, j));
            let v = v - u * f[2];
            [u, v, w - u * f[3] - v * g[3]]
        })
    };
    let t2 = hecke_matrix(2);
    let [a, b, c] = t2;
    let trace = a[0] + b[1] + c[2];
    let second = a[0] * b[1] - a[1] * b[0] + a[0] * c[2] - a[2] * c[0] + b[1] * c[2] - b[2] * c[1];
    let determinant = a[0] * (b[1] * c[2] - b[2] * c[1]) - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1]);
    let t3 = hecke_matrix(3);
    let t2_square: [[i128; 3]; 3] = t2
        .map(|column| std::array::from_fn(|i| (0..3).map(|j| t2[j][i] * column[j]).sum::<i128>()));
    for (j, column) in t3.iter().enumerate() {
        for (i, &value) in column.iter().enumerate() {
            // An independent Fourier check of 72 T_3 in the basis I,T_2,T_2².
            assert_eq!(
                72 * value,
                72 * (478011548 * i128::from(i == j) + 2697 * t2[j][i]) - t2_square[j][i],
            );
        }
    }
    let report = stdout(interpret(
        "let a = rvHeckeT[36,2]; rvHeckeChar(a); rvHeckeCharCoeff[3](a);
         rvheckepacket E(a);
         let a3 = rvHeckeEigenCoeff[E,3];
         rvHeckeEigenCoord[0](a3); rvHeckeEigenCoord[1](a3); rvHeckeEigenCoord[2](a3);
         rvHeckeEigenNorm(rvHeckeEigenGenerator[E]);",
    ));
    assert_eq!(
        report,
        format!(
            "rvHeckeChar = [1, {}, {second}, {}]\nrvHeckeCharCoeff[3] = {}\n\
         rvHeckeEigenCoord[0] = 478011548\nrvHeckeEigenCoord[1] = 2697\n\
         rvHeckeEigenCoord[2] = -1/72\nrvHeckeEigenNorm = {determinant}\n",
            -trace, -determinant, -determinant,
        )
    );
}

#[test]
fn hecke_characteristic_scalar_shortcuts_preserve_large_weight_and_empty_space() {
    let weight = u128::MAX - 1;
    let degree = weight / 12 - u128::from(weight % 12 == 2);
    let index = u128::MAX;
    let source = format!(
        "let t = rvHeckeT[2,{index}];
         let z = zero[RVHecke(2)];
         assert rvHeckeMul(t,rvHeckeUnit[2]) == t;
         assert rvHeckeSameAction(t,z) == S.hi;
         rvHeckeActionTrace(t);
         rvHeckeChar(t);
         let a = rvHeckeUnit[{weight}];
         rvHeckeActionTrace(a);
         rvHeckeCharCoeff[1](a);
         rvHeckeCharCoeff[{degree}](neg(a));
         rvHeckeCharCoeff[2](zero[RVHecke({weight})]);
         rvHeckeChar(rvHeckeScale[5](rvHeckeUnit[36]));
         rvHeckeChar(zero[RVHecke(36)]);"
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "rvHeckeActionTrace = 0\nrvHeckeChar = [1]\nrvHeckeActionTrace = {degree}\n\
             rvHeckeCharCoeff[1] = -{degree}\nrvHeckeCharCoeff[{degree}] = 1\n\
             rvHeckeCharCoeff[2] = 0\nrvHeckeChar = [1, -15, 75, -125]\nrvHeckeChar = [1, 0, 0, 0]\n"
        )
    );
}

#[test]
fn hecke_operators_diagnose_invalid_grades_and_arithmetic_limits() {
    let index = u128::MAX;
    for (source, expected) in [
        (
            "show rvHeckeT[3,1];".to_owned(),
            "even level-one weight at least 2",
        ),
        (
            "show zero[RVHecke(0)];".to_owned(),
            "even level-one weight at least 2",
        ),
        (
            "show rvHeckeT[12,0];".to_owned(),
            "Hecke index n must be positive",
        ),
        (
            "show rvHeckeMul(rvHeckeT[12,5],rvHeckeT[24,5]);".to_owned(),
            "RVHecke(12)",
        ),
        (
            "show rvHeckeSameAction(rvHeckeT[12,5],rvHeckeT[24,5]);".to_owned(),
            "RVHecke(12)",
        ),
        (
            "show add(rvHeckeT[12,5],rvHeckeT[24,5]);".to_owned(),
            "RVHecke(12)",
        ),
        (
            format!("show rvHeckeMul(rvHeckeT[12,{index}],rvHeckeT[12,2]);"),
            "index product exceeds unsigned 128-bit range",
        ),
        (
            "show rvHeckeMul(rvHeckeT[514,2],rvHeckeT[514,2]);".to_owned(),
            "exact intermediate exceeds 512 bits",
        ),
        (
            format!("rvHeckeChar(zero[RVHecke({})]);", index - 1),
            "exceeds the coefficient storage limit",
        ),
    ] {
        let error = diagnostic(interpret(&source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}

#[test]
fn eigenpacket_coordinates_and_conjugate_invariants_match_q_expansion() {
    let bound = 16;
    let delta = delta(bound);
    let e4 = e4(bound);
    let e4_cubed = convolve(&convolve(&e4, &e4, bound), &e4, bound);
    let f = convolve(&delta, &e4_cubed, bound);
    let g = convolve(&delta, &delta, bound);
    let a = hecke24_matrix(&f, &g, 2);
    assert_eq!(a[1][0], 1);
    let mut program = String::from("rvheckepacket E(rvHeckeT[24,2]);\n");
    let mut expected = String::new();
    for n in 1..=8 {
        let b = hecke24_matrix(&f, &g, n);
        let v = b[1][0];
        let u = b[0][0] - v * a[0][0];
        // These coordinates reproduce every entry in the independently
        // constructed T_n matrix, so they give both conjugate eigenvalues.
        let reconstructed =
            std::array::from_fn(|i| std::array::from_fn(|j| u * i128::from(i == j) + v * a[i][j]));
        assert_eq!(b, reconstructed);
        let [_, negative_trace, norm] = characteristic_2x2(b);
        writeln!(
            program,
            "let a{n} = rvHeckeEigenCoeff[E,{n}];
             rvHeckeEigenCoord[0](a{n}); rvHeckeEigenCoord[1](a{n});
             rvHeckeEigenTrace(a{n}); rvHeckeEigenNorm(a{n});",
        )
        .unwrap();
        writeln!(
            expected,
            "rvHeckeEigenCoord[0] = {u}\nrvHeckeEigenCoord[1] = {v}\nrvHeckeEigenTrace = {}\nrvHeckeEigenNorm = {norm}",
            -negative_trace
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&program)), expected);
}

#[test]
fn eigenpacket_arithmetic_preserves_hecke_products_and_rational_coordinates() {
    assert_eq!(
        stdout(example("hecke-eigenpackets")),
        concat!(
            "rvHeckeEigenData[E] = weight 24, [1, -1080, -20468736]\n",
            "RVHeckeEigen(E){[195660, -48]/1}\n",
            "rvHeckeEigenTrace = 339480\n",
            "rvHeckeEigenNorm = -19020146544\n",
            "RVHeckeEigen(F){[978636, -48]/5}\n",
            "rvHeckeEigenCoord[1] = -48/5\n",
        )
    );
    let source = "
        let a = rvHeckeT[24,2];
        rvheckepacket E(add(rvHeckeScale[2](a),rvHeckeScale[7](rvHeckeUnit[24])));
        let x = rvHeckeEigenGenerator[E];
        let half = rvHeckeEigenScalar[E,1,2];
        let a2 = rvHeckeEigenCoeff[E,2];
        assert a2 == rvHeckeEigenMul(half,add(x,rvHeckeEigenScalar[E,-7,1]));
        assert rvHeckeEigenMul(a2,zero[RVHeckeEigen(E)]) == zero[RVHeckeEigen(E)];
        assert add(a2,neg(a2)) == zero[RVHeckeEigen(E)];
        assert add(half,rvHeckeEigenScalar[E,2,4]) == rvHeckeEigenCoeff[E,1];
        let b = add(rvHeckeT[24,3],rvHeckeT[24,4]);
        assert rvHeckeEigenImage[E](rvHeckeMul(a,b)) ==
            rvHeckeEigenMul(a2,add(rvHeckeEigenCoeff[E,3],rvHeckeEigenCoeff[E,4]));
        assert rvHeckeEigenImage[E](add(rvHeckeMul(a,a),
            add(rvHeckeScale[-1080](a),rvHeckeScale[-20468736](rvHeckeUnit[24]))))
            == zero[RVHeckeEigen(E)];
        rvHeckeEigenNorm(half);
        rvHeckeEigenTrace(half);
        rvHeckeEigenNorm(zero[RVHeckeEigen(E)]);
        rvHeckeEigenCoord[1](a2);
        rvHeckeEigenCoord[2](a2);
        rvHeckeEigenNorm(a2);";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvHeckeEigenNorm = 1/4\nrvHeckeEigenTrace = 1\n",
            "rvHeckeEigenNorm = 0\nrvHeckeEigenCoord[1] = 1/2\n",
            "rvHeckeEigenCoord[2] = 0\nrvHeckeEigenNorm = -20468736\n",
        )
    );
}

#[test]
fn degree_one_eigenpackets_return_tau_and_reduce_the_generator() {
    let source = "
        rvheckepacket Delta(rvHeckeT[12,2]);
        rvHeckeEigenData[Delta];
        assert rvHeckeEigenGenerator[Delta] == rvHeckeEigenScalar[Delta,-24,1];
        assert rvHeckeEigenCoeff[Delta,1] == rvHeckeEigenScalar[Delta,1,1];
        assert rvHeckeEigenCoeff[Delta,5] == rvHeckeEigenScalar[Delta,4830,1];
        assert rvHeckeEigenCoeff[Delta,6] == rvHeckeEigenScalar[Delta,-6048,1];
        assert rvHeckeEigenCoeff[Delta,9] == rvHeckeEigenScalar[Delta,-113643,1];
        rvHeckeEigenNorm(rvHeckeEigenCoeff[Delta,5]);
        rvHeckeEigenTrace(rvHeckeEigenCoeff[Delta,5]);";
    assert_eq!(
        stdout(interpret(source)),
        "rvHeckeEigenData[Delta] = weight 12, [1, 24]\nrvHeckeEigenNorm = 4830\nrvHeckeEigenTrace = 4830\n",
    );
}

#[test]
fn eigenpacket_certificate_detects_repeated_eigenvalues_and_nominal_mismatch() {
    for (source, expected) in [
        (
            "rvheckepacket E(rvHeckeUnit[2]);",
            "requires a nonzero cusp space",
        ),
        (
            "rvheckepacket E(rvHeckeUnit[24]);",
            "squarefree characteristic polynomial",
        ),
        (
            "let t = rvHeckeT[24,2]; rvheckepacket E(add(rvHeckeMul(t,t),rvHeckeScale[-1080](t)));",
            "squarefree characteristic polynomial",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); show rvHeckeEigenImage[E](rvHeckeT[12,2]);",
            "RVHecke(24)",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); rvheckepacket F(rvHeckeT[24,2]);
             show rvHeckeEigenMul(rvHeckeEigenCoeff[E,2],rvHeckeEigenCoeff[F,2]);",
            "RVHeckeEigen(E)",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); rvheckepacket E(rvHeckeT[24,2]);",
            "already defined",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); show rvHeckeEigenScalar[E,1,0];",
            "denominator must be positive",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); show rvHeckeEigenCoeff[E,0];",
            "Hecke index n must be positive",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}

#[test]
fn sturm_certifies_classical_congruences_and_distinguishes_residual_branches() {
    assert_eq!(
        stdout(example("hecke-eisenstein-congruences")),
        concat!(
            "rvHeckeEisensteinSturm[R] = congruent modulo 691 (Sturm bound 1)\n",
            "rvHeckeEisensteinSturm[S16] = congruent modulo 3617 (Sturm bound 1)\n",
            "rvHeckeEisensteinSturm[Good] = congruent modulo 103 (Sturm bound 2)\n",
            "rvHeckeEisensteinSturm[Other] = mismatch at q^2 modulo 103 (eigenform 70, Eisenstein 83; Sturm bound 2)\n",
            "rvHeckeEisensteinSturm[Large] = congruent modulo 2294797 (Sturm bound 2)\n",
            "PrimeScalar(103){96}\n",
        )
    );
    assert_eq!(
        stdout(interpret(
            "rvheckepacket W(rvHeckeT[24,2]);
             rvheckeresidual R[W,2294797,791659]; rvHeckeEisensteinSturm[R];
             rvheckepacket D(rvHeckeT[12,2]);
             rvheckeresidual X[D,13,-24]; rvHeckeEisensteinSturm[X];
             rvheckeresidual Y[D,17,-24]; rvHeckeEisensteinSturm[Y];"
        )),
        concat!(
            "rvHeckeEisensteinSturm[R] = mismatch at q^2 modulo 2294797 (eigenform 791659, Eisenstein 1504218; Sturm bound 2)\n",
            "rvHeckeEisensteinSturm[X] = unavailable modulo 13: Eisenstein coefficient q^0 is nonintegral\n",
            "rvHeckeEisensteinSturm[Y] = mismatch at q^0 modulo 17 (eigenform 0, Eisenstein 14; Sturm bound 1)\n",
        )
    );
}

fn bernoulli_by_recurrence(weight: usize) -> (i128, i128) {
    fn gcd(mut a: i128, mut b: i128) -> i128 {
        a = a.abs();
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    }
    fn reduce(a: i128, b: i128) -> (i128, i128) {
        let g = gcd(a, b);
        (a / g, b / g)
    }
    // Independent Bernoulli recurrence from t/(exp(t)-1).
    let mut numbers = vec![(1, 1)];
    for n in 1..=weight {
        let (mut a, mut b) = (0, 1);
        let mut binomial = 1;
        for (j, &(x, y)) in numbers.iter().enumerate() {
            (a, b) = reduce(a * y + b * binomial * x, b * y);
            binomial = binomial * (n + 1 - j) as i128 / (j + 1) as i128;
        }
        numbers.push(reduce(-a, b * (n + 1) as i128));
    }
    numbers[weight]
}

#[test]
fn residual_bernoulli_constants_match_independent_rational_recurrence() {
    let mut program = String::new();
    let mut expected = String::new();
    for weight in (4..=26).step_by(2) {
        let (mut a, mut b) = bernoulli_by_recurrence(weight);
        a = -a;
        b *= 2 * weight as i128;
        for p in [5, 7, 11, 13, 17, 19, 103, 691] {
            let (mut num, mut den) = (a, b);
            while den % p == 0 && num % p == 0 {
                (num, den) = (num / p, den / p);
            }
            let source = format!("show rvHeckeEisensteinCoeff[{weight},{p},0];");
            if den % p == 0 {
                assert!(
                    diagnostic(interpret(&source)).contains("is nonintegral"),
                    "{source}"
                );
            } else {
                let result = (0..p).find(|r| (r * den - num) % p == 0).unwrap();
                writeln!(program, "{source}").unwrap();
                writeln!(expected, "PrimeScalar({p}){{{result}}}").unwrap();
            }
        }
    }
    assert_eq!(stdout(interpret(&program)), expected);
}

#[test]
fn residual_fourier_values_match_independent_modular_form_expansions_beyond_sturm() {
    let bound = 14;
    let d = delta(bound);
    let e = e4(bound);
    let e_cubed = convolve(&convolve(&e, &e, bound), &e, bound);
    let w16 = convolve(&d, &e, bound);
    let f24 = convolve(&d, &e_cubed, bound);
    let g24 = convolve(&d, &d, bound);
    let mut program = String::from(
        "rvheckepacket D(rvHeckeT[12,2]);
         rvheckepacket M(rvHeckeT[16,2]);
         rvheckepacket W(rvHeckeT[24,2]);
         rvheckeresidual Alpha[D,691,-24];
         rvheckeresidual Beta[M,3617,216];
         rvheckeresidual Gamma[W,103,83];
         rvheckeresidual Other[W,103,70];
         rvheckeresidual Large[W,2294797,1504218];",
    );
    let mut expected = String::new();
    for n in 1..=bound {
        for (name, weight, p, a_n, congruent) in [
            ("Alpha", 12, 691, d[n], true),
            ("Beta", 16, 3617, w16[n], true),
            ("Gamma", 24, 103, f24[n] + (83 - 696) * g24[n], true),
            ("Other", 24, 103, f24[n] + (70 - 696) * g24[n], false),
            (
                "Large",
                24,
                2294797,
                f24[n] + (1504218 - 696) * g24[n],
                true,
            ),
        ] {
            writeln!(program, "show rvHeckeResidualCoeff[{name},{n}];").unwrap();
            let residue = a_n.rem_euclid(p);
            writeln!(expected, "PrimeScalar({p}){{{residue}}}").unwrap();
            if congruent {
                let sigma = (1..=n)
                    .filter(|d| n.is_multiple_of(*d))
                    .map(|d| (d as i128).pow(weight - 1))
                    .sum::<i128>();
                assert_eq!(residue, sigma.rem_euclid(p), "{name} coefficient {n}");
                writeln!(program, "assert rvHeckeResidualCoeff[{name},{n}] == rvHeckeEisensteinCoeff[{weight},{p},{n}];").unwrap();
            }
        }
    }
    assert_eq!(stdout(interpret(&program)), expected);
}

#[test]
fn hensel_evaluation_preserves_integral_hecke_values_with_nonunit_coordinates() {
    let source = "
        let a = rvHeckeT[36,2];
        rvheckepacket W(add(rvHeckeMul(a,a),rvHeckeScale[-4](a)));
        rvheckeresidual R[W,11,10];
        let value = rvHeckeEigenCoeff[W,2];
        show value;
        show rvHeckeResidualReduce[R](value);
        show rvHeckeResidualCoeff[R,2];
        rvheckepacket Q(rvHeckeT[24,2]);
        rvheckeresidual S[Q,103,83];
        let x = rvHeckeEigenGenerator[Q];
        let p2 = rvHeckeEigenMul(add(x,rvHeckeEigenScalar[Q,-83,1]),rvHeckeEigenScalar[Q,1,103]);
        let p3 = rvHeckeEigenMul(add(x,rvHeckeEigenScalar[Q,-4306,1]),rvHeckeEigenScalar[Q,1,10609]);
        show rvHeckeResidualReduce[S](p2);
        show rvHeckeResidualReduce[S](p3);";
    let report = stdout(interpret(source));
    let mut lines = report.lines();
    let power_basis = lines.next().unwrap();
    let denominator = power_basis
        .rsplit_once('/')
        .unwrap()
        .1
        .strip_suffix('}')
        .unwrap()
        .parse::<i128>()
        .unwrap();
    assert_eq!(denominator % 11, 0, "{power_basis}");
    assert_eq!(
        lines.collect::<Vec<_>>(),
        [
            "PrimeScalar(11){7}",
            "PrimeScalar(11){7}",
            "PrimeScalar(103){41}",
            "PrimeScalar(103){16}",
        ]
    );
}

#[test]
fn residual_maps_reject_invalid_roots_nominal_mismatches_and_precision_limits() {
    let delta_packet = "rvheckepacket D(rvHeckeT[12,2]);";
    let cubic =
        "let a = rvHeckeT[36,2]; rvheckepacket W(add(rvHeckeMul(a,a),rvHeckeScale[-4](a)));";
    for (source, expected) in [
        (format!("{delta_packet} rvheckeresidual R[D,689,-24];"), "not prime"),
        (format!("{delta_packet} rvheckeresidual R[D,3,-24];"), "prime at least 5"),
        (format!("{delta_packet} rvheckeresidual R[D,691,24];"), "is not a root"),
        (format!("{cubic} rvheckeresidual R[W,11,8];"), "root must be simple"),
        (
            format!("{delta_packet} rvheckeresidual R[D,691,-24]; rvheckepacket E(rvHeckeT[12,2]); show rvHeckeResidualReduce[R](rvHeckeEigenCoeff[E,2]);"),
            "RVHeckeEigen(D)",
        ),
        (
            format!("{delta_packet} rvheckeresidual R[D,691,-24]; show rvHeckeResidualReduce[R](rvHeckeEigenScalar[D,1,691]);"),
            "is nonintegral",
        ),
        (
            format!("{delta_packet} rvheckeresidual R[D,691,-24]; show rvHeckeResidualReduce[R](rvHeckeEigenScalar[D,1,{}]);", 691_i128.pow(7)),
            "prime-power precision",
        ),
    ] {
        let error = diagnostic(interpret(&source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}

#[test]
fn quadratic_order_example_gives_integral_coordinates_and_sharp_local_depths() {
    assert_eq!(
        stdout(example("quadratic-hecke-orders")),
        concat!(
            "rvHeckeQuadraticOrder[E] = D_K=144169, omega=(1+sqrt(144169))/2, theta=(528, 24), theta_index=24, theta_discriminant=83041344, conductor=24, order_discriminant=83041344, congruence_ideal=(24*sqrt(144169)), congruence_norm=83041344\n",
            "rvHeckeQuadraticCoeff[E, 2] = (528, 24)\n",
            "rvHeckeQuadraticCoeff[E, 3] = (170316, -1152)\n",
            "rvHeckeQuadraticLocal[E, 2] = split, e=1, f=1, primes=2, depth=3, first_sharp_index=2\n",
            "rvHeckeQuadraticLocal[E, 3] = split, e=1, f=1, primes=2, depth=1, first_sharp_index=2\n",
            "rvHeckeQuadraticLocal[E, 144169] = ramified, e=2, f=1, primes=1, depth=1, first_sharp_index=2\n",
            "rvHeckeQuadraticOrder[F] = D_K=144169, omega=(1+sqrt(144169))/2, theta=(2647, 120), theta_index=120, theta_discriminant=2076033600, conductor=24, order_discriminant=83041344, congruence_ideal=(24*sqrt(144169)), congruence_norm=83041344\n",
            "rvHeckeQuadraticCoeff[F, 2] = (528, 24)\n",
        )
    );
}

fn series_power(series: &[i128], exponent: u32, bound: usize) -> Vec<i128> {
    let mut result = vec![0; bound + 1];
    result[0] = 1;
    for _ in 0..exponent {
        result = convolve(&result, series, bound);
    }
    result
}

fn integral_cusp_basis(weight: u32, bound: usize) -> [Vec<i128>; 2] {
    let delta = delta(bound);
    let e4 = e4(bound);
    let e6: Vec<i128> = std::iter::once(1)
        .chain((1..=bound).map(|n| {
            -504 * (1..=n)
                .filter(|d| n.is_multiple_of(*d))
                .map(|d| (d as i128).pow(5))
                .sum::<i128>()
        }))
        .collect();
    [1, 2].map(|j| {
        let b = (weight - 12 * j) % 4 / 2;
        let a = (weight - 12 * j - 6 * b) / 4;
        convolve(
            &series_power(&delta, j, bound),
            &convolve(
                &series_power(&e4, a, bound),
                &series_power(&e6, b, bound),
                bound,
            ),
            bound,
        )
    })
}

#[test]
fn quadratic_orders_in_all_degree_two_weights_match_independent_integral_cusp_bases() {
    let fixtures = [
        (24, 144169, 24),
        (28, 18209, 216),
        (30, 51349, 192),
        (32, 18295489, 24),
        (34, 2356201, 144),
        (38, 63737521, 96),
    ];
    let bound = 10;
    let mut source = String::new();
    let mut expected = String::new();
    for (k, d, g) in fixtures {
        let [f, h] = integral_cusp_basis(k, bound);
        assert_eq!([f[1], h[1], h[2]], [1, 0, 1]);
        let image_col = |form: &[i128]| {
            let first = hecke_q_coefficient(form, k, 2, 1);
            [first, hecke_q_coefficient(form, k, 2, 2) - f[2] * first]
        };
        let [_, negative_trace, norm] = characteristic_2x2([image_col(&f), image_col(&h)]);
        let trace = -negative_trace;
        let discriminant = trace * trace - 4 * norm;
        assert_eq!(discriminant, g * g * d);
        let offset = (trace - g) / 2;
        writeln!(
            source,
            "rvheckepacket E{k}(rvHeckeT[{k},2]); rvHeckeQuadraticOrder[E{k}];"
        )
        .unwrap();
        writeln!(expected, "rvHeckeQuadraticOrder[E{k}] = D_K={d}, omega=(1+sqrt({d}))/2, theta=({offset}, {g}), theta_index={g}, theta_discriminant={discriminant}, conductor={g}, order_discriminant={discriminant}, congruence_ideal=({g}*sqrt({d})), congruence_norm={discriminant}").unwrap();
        for n in 0..=bound {
            writeln!(source, "rvHeckeQuadraticCoeff[E{k},{n}];").unwrap();
            let a = f[n] + (offset - f[2]) * h[n];
            let b = g * h[n];
            writeln!(expected, "rvHeckeQuadraticCoeff[E{k}, {n}] = ({a}, {b})").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn field_orientation_and_coefficient_order_survive_large_generator_changes() {
    let scale = 100_000_000_000_000_000_000_i128;
    let program = format!(
        "let t = rvHeckeT[24,2]; let one = rvHeckeUnit[24];
         rvheckepacket Huge(add(rvHeckeScale[{scale}](t),rvHeckeScale[7](one)));
         rvheckepacket Neg(add(rvHeckeScale[-5](t),rvHeckeScale[7](one)));
         rvHeckeQuadraticOrder[Huge]; rvHeckeQuadraticOrder[Neg];
         rvHeckeQuadraticCoeff[Huge,3]; rvHeckeQuadraticCoeff[Neg,3];
         rvHeckeQuadraticIntegral(rvHeckeEigenCoeff[Huge,2]);
         let x = rvHeckeEigenGenerator[Huge];
         let a = add(x,rvHeckeEigenScalar[Huge,{},1]);
         rvHeckeQuadraticIntegral(rvHeckeEigenMul(a,rvHeckeEigenScalar[Huge,1,{}]));",
        -(528 * scale + 7),
        24 * scale
    );
    let report = stdout(interpret(&program));
    let lines: Vec<_> = report.lines().collect();
    assert!(lines[0].contains(&format!("theta_index={},", 24 * scale)));
    assert!(lines[1].contains("theta=(-2753, 120), theta_index=120,"));
    assert!(lines[..2]
        .iter()
        .all(|line| line.contains("conductor=24, order_discriminant=83041344")));
    assert_eq!(
        lines[2..],
        [
            "rvHeckeQuadraticCoeff[Huge, 3] = (170316, -1152)",
            "rvHeckeQuadraticCoeff[Neg, 3] = (169164, 1152)",
            "rvHeckeQuadraticIntegral = (528, 24)",
            "rvHeckeQuadraticIntegral = (0, 1)",
        ]
    );
}

#[test]
fn quadratic_local_splitting_and_congruence_witnesses_handle_inert_and_ramified_primes() {
    let report = stdout(interpret(
        "rvheckepacket E(rvHeckeT[24,2]); rvHeckeQuadraticLocal[E,23];
         rvHeckeQuadraticLocal[E,4294967311];
         rvheckepacket F(rvHeckeT[28,2]); rvHeckeQuadraticLocal[F,3];
         rvHeckeQuadraticLocal[F,131]; rvHeckeQuadraticLocal[F,139];
         rvheckepacket G(rvHeckeT[30,2]); rvHeckeQuadraticLocal[G,2];",
    ));
    assert_eq!(report, concat!(
        "rvHeckeQuadraticLocal[E, 23] = inert, e=1, f=2, primes=1, depth=0, first_sharp_index=2\n",
        "rvHeckeQuadraticLocal[E, 4294967311] = inert, e=1, f=2, primes=1, depth=0, first_sharp_index=2\n",
        "rvHeckeQuadraticLocal[F, 3] = inert, e=1, f=2, primes=1, depth=3, first_sharp_index=2\n",
        "rvHeckeQuadraticLocal[F, 131] = ramified, e=2, f=1, primes=1, depth=1, first_sharp_index=2\n",
        "rvHeckeQuadraticLocal[F, 139] = ramified, e=2, f=1, primes=1, depth=1, first_sharp_index=2\n",
        "rvHeckeQuadraticLocal[G, 2] = inert, e=1, f=2, primes=1, depth=6, first_sharp_index=2\n",
    ));
}

#[test]
fn quadratic_orders_validate_degree_local_primes_and_maximal_integrality() {
    for (source, expected) in [
        (
            "rvheckepacket E(rvHeckeT[12,2]); rvHeckeQuadraticOrder[E];",
            "requires a packet of degree two",
        ),
        (
            "rvheckepacket E(rvHeckeT[36,2]); rvHeckeQuadraticOrder[E];",
            "requires a packet of degree two",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); rvHeckeQuadraticLocal[E,9];",
            "not prime",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); rvHeckeQuadraticLocal[E,4];",
            "not prime",
        ),
        (
            "rvheckepacket E(rvHeckeT[24,2]); rvHeckeQuadraticIntegral(rvHeckeEigenScalar[E,1,2]);",
            "is nonintegral in the quadratic maximal order",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}: {error}");
    }
}
