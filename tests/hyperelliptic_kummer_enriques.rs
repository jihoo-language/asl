mod common;

use common::{diagnostic, example, interpret, stdout};

const SHEAR: &str = r"
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

const CUBIC: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,2,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=neg(o);
    rvhyperelliptic C[F]{z,t,t,t,t,o};
    rvjacbase L[C,W];
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let a=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[L](zw,zw);
    let f1=rvJacPoint[L](ow,zw);
    let e2=add(add(e1,f1),rvJacPoint[L](a,zw));
    let f2=add(add(e1,f1),rvJacPoint[L](rvUnramFrob[1](a),zw));
    rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};
";

#[test]
fn quotient_counts_include_nonsplit_fibers_and_normalized_zeta() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-enriques")),
        concat!(
            "rvKummerEnriquesPoints[B, C, 2, 1] = points=22, normalized trace=4, split=16, nonsplit=6, twisted cover=12\n",
            "rvKummerEnriquesPoints[B, C, 2, 2] = points=172, normalized trace=10, split=124, nonsplit=48, twisted cover=96\n",
            "rvKummerEnriquesPoints[B, L, 2, 1] = points=172, normalized trace=10, split=124, nonsplit=48, twisted cover=96\n",
            "rvKummerEnriquesZeta[B, C, 2] = 1/((1-T)R(3*T)(1-9*T)), R=[1, -4, 3, 8, -14, 0, 14, -8, -3, 4, -1]\n",
            "rvKummerEnriquesZeta[B, L, 2] = 1/((1-T)R(9*T)(1-81*T)), R=[1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1]\n",
        )
    );
}

#[test]
fn actual_even_switch_and_twisted_frobenius_give_the_reported_fibers() {
    let output = stdout(interpret(&format!(
        r"{SHEAR}
        rvKummerEnriquesPoints[B,C,2,1];
        let f=\x:RVKummerK3(L)=>rvKummerK3Frob[C,1](x);
        let s=\x:RVKummerK3(L)=>rvKummerK3Translate[B,2](rvKummerSwitch(x));
        let twisted=\x:RVKummerK3(L)=>s(f(x));
        assert twisted==(\x:RVKummerK3(L)=>f(s(x)));
        show f;
        show twisted;
        show s;
        "
    )));
    let mut lines = output.lines();
    assert_eq!(lines.next().unwrap(), "rvKummerEnriquesPoints[B, C, 2, 1] = points=22, normalized trace=4, split=16, nonsplit=6, twisted cover=12");
    for expected in [32, 12, 0] {
        let map = lines.next().unwrap();
        let points: Vec<_> = map
            .match_indices("RVKummerK3(L){")
            .map(|(start, _)| {
                let end = start + map[start..].find('}').unwrap() + 1;
                &map[start..end]
            })
            .collect();
        assert_eq!(points.len(), 2 * 248);
        assert_eq!(
            points
                .chunks_exact(2)
                .filter(|pair| pair[0] == pair[1])
                .count(),
            expected
        );
    }
    assert!(lines.next().is_none());
}

#[test]
fn point_counts_track_order_three_frobenius() {
    let output = stdout(interpret(&format!(
        r"{CUBIC}
        rvKummerEnriquesPoints[B,C,5,1];
        rvKummerEnriquesPoints[B,C,5,2];
        rvKummerEnriquesPoints[B,C,5,3];
        rvKummerEnriquesZeta[B,C,5];
        "
    )));
    let mut lines = output.lines();
    for (degree, points, trace) in [(1, 13, 1), (2, 91, 1), (3, 1000, 10)] {
        assert!(
            lines.next().unwrap().starts_with(&format!(
            "rvKummerEnriquesPoints[B, C, 5, {degree}] = points={points}, normalized trace={trace},"
        )),
            "{output}"
        );
    }
    assert_eq!(lines.next().unwrap(), "rvKummerEnriquesZeta[B, C, 5] = 1/((1-T)R(3*T)(1-9*T)), R=[1, -1, 0, -3, 3, 0, 3, -3, 0, -1, 1]");
    assert!(lines.next().is_none());
}

#[test]
fn rational_branch_points_give_full_trace_and_half_the_actual_k3_fiber_count() {
    let output = stdout(interpret(
        r"
        rvunram F[7,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let t=rvUnramElem[F]{2};
        let h=rvUnramElem[F]{4};
        rvhyperelliptic C[F]{z,h,z,t,z,o};
        let e1=rvJacPoint[C](z,z);
        let f1=rvJacPoint[C](o,z);
        let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
        let f2=add(add(e1,f1),rvJacPoint[C](t,z));
        rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};
        rvKummerEnriquesPoints[B,C,2,1];
        show (\x:RVKummerK3(C)=>rvKummerIsExceptional(x));
        ",
    ));
    let mut lines = output.lines();
    let report = lines.next().unwrap();
    let cover_points = lines.next().unwrap().matches("RVKummerK3(C)").count();
    assert!(cover_points.is_multiple_of(2));
    let split = cover_points / 2;
    let nonsplit = 120 - split;
    assert_eq!(report, format!("rvKummerEnriquesPoints[B, C, 2, 1] = points=120, normalized trace=10, split={split}, nonsplit={nonsplit}, twisted cover={}", 2 * nonsplit));
    assert!(nonsplit > 0);
}

#[test]
fn invariant_zeta_agrees_with_the_determinant_of_the_projected_divisor_action() {
    for (source, columns, mask) in [(SHEAR, [1, 2, 4, 10], 2), (CUBIC, [1, 8, 4, 10], 5)] {
        let output = stdout(interpret(&format!(
            "{source} rvKummerEnriquesZeta[B,C,{mask}];"
        )));
        let (_, coefficients) = output.trim().split_once("R=[").unwrap();
        let coefficients: Vec<i64> = coefficients
            .strip_suffix(']')
            .unwrap()
            .split(", ")
            .map(|value| value.parse().unwrap())
            .collect();
        assert_eq!(coefficients.len(), 11);
        let image = |node: usize| {
            columns
                .iter()
                .enumerate()
                .filter(|(i, _)| node & (1 << i) != 0)
                .fold(0, |sum, (_, value)| sum ^ value)
        };
        // Compute det(1−t F (1+σ)/2) directly on the whole 17-dimensional span
        // modulo 101. The complement of the projector contributes only ones.
        const INV4: i64 = 76;
        for t in 0..=10 {
            let mut matrix = [[0; 17]; 17];
            for (i, row) in matrix.iter_mut().enumerate() {
                row[i] = 1;
            }
            matrix[0][0] -= 8 * INV4 * t;
            for node in 0..16 {
                let row = image(node) + 1;
                matrix[0][node + 1] -= INV4 * t;
                matrix[row][0] += 2 * INV4 * t;
                matrix[row][node + 1] -= 2 * INV4 * t;
                for other in 0..16 {
                    if [0, 1, 4, 7, 13, 15].contains(&(other ^ node ^ mask)) {
                        matrix[image(other) + 1][node + 1] += INV4 * t;
                    }
                }
            }
            let expected = determinant_mod_101(matrix);
            let actual = coefficients.iter().rev().fold(0, |value, coefficient| {
                (t * value + coefficient).rem_euclid(101)
            });
            assert_eq!(actual, expected, "mask {mask}, parameter {t}");
        }
    }
}

fn determinant_mod_101(mut matrix: [[i64; 17]; 17]) -> i64 {
    let mut result = 1;
    for column in 0..17 {
        let Some(pivot_row) = (column..17).find(|&row| matrix[row][column] % 101 != 0) else {
            return 0;
        };
        if pivot_row != column {
            matrix.swap(pivot_row, column);
            result = -result;
        }
        let pivot = matrix[column][column].rem_euclid(101);
        result = (result * pivot).rem_euclid(101);
        let inverse = (1..101).find(|&i| (i * pivot) % 101 == 1).unwrap();
        for row in (column + 1)..17 {
            let factor = (inverse * matrix[row][column]).rem_euclid(101);
            for j in column..17 {
                matrix[row][j] = (matrix[row][j] - factor * matrix[column][j]).rem_euclid(101);
            }
        }
    }
    result
}

#[test]
fn rejects_odd_nonfixed_bad_base_and_unbounded_inputs() {
    for (source, expected) in [
        ("rvKummerEnriquesPoints[B,C,2,0];", "extension degree must be positive"),
        ("rvKummerEnriquesZeta[B,C,16];", "between 0 and 15"),
        ("rvKummerEnriquesZeta[B,C,13];", "requires an even theta characteristic"),
        ("rvKummerEnriquesPoints[B,C,8,2];", "must be fixed by base Frobenius"),
        ("rvKummerEnriquesPoints[B,C,2,340282366920938463463374607431768211455];", "point count exceeds unsigned 128-bit range"),
        ("rvhyperelliptic D[F]{z,o,z,z,z,o}; rvKummerEnriquesZeta[B,D,2];", "not the declared base extension"),
        ("rvhyperelliptic A[F]{z,neg(o),z,o}; rvjactorsionbasis J[A,2](neg(o)){rvJacPoint[A](z,z),rvJacPoint[A](o,z)}; rvKummerEnriquesZeta[J,C,2];", "certified genus-two level-2 basis"),
    ] {
        let error = diagnostic(interpret(&format!("{SHEAR}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
