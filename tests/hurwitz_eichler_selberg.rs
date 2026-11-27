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
fn weight_thirty_six_characteristic_uses_three_exact_newton_steps() {
    let bound = 6;
    let delta = delta(bound);
    let e4 = e4(bound);
    let e4_cubed = convolve(&convolve(&e4, &e4, bound), &e4, bound);
    let delta_squared = convolve(&delta, &delta, bound);
    let f = convolve(&delta, &convolve(&e4_cubed, &e4_cubed, bound), bound);
    let g = convolve(&delta_squared, &e4_cubed, bound);
    let h = convolve(&delta_squared, &delta, bound);
    assert_eq!([f[1], g[1], g[2], h[1], h[2], h[3]], [1, 0, 1, 0, 0, 1]);

    let [a, b, c] = [&f, &g, &h].map(|form| {
        let [u, v, w] = [1, 2, 3].map(|j| hecke_q_coefficient(form, 36, 2, j));
        let v = v - u * f[2];
        [u, v, w - u * f[3] - v * g[3]]
    });
    let trace = a[0] + b[1] + c[2];
    let second = a[0] * b[1] - a[1] * b[0] + a[0] * c[2] - a[2] * c[0] + b[1] * c[2] - b[2] * c[1];
    let determinant = a[0] * (b[1] * c[2] - b[2] * c[1]) - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1]);
    let report = stdout(interpret(
        "let a = rvHeckeT[36,2]; rvHeckeChar(a); rvHeckeCharCoeff[3](a);",
    ));
    assert_eq!(
        report,
        format!(
            "rvHeckeChar = [1, {}, {second}, {}]\nrvHeckeCharCoeff[3] = {}\n",
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
