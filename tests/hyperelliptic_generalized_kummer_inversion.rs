mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const ODD: &str = include_str!("../examples/hyperelliptic-generalized-kummer-inversion.asl");
const EVEN: &str = include_str!("../examples/hyperelliptic-generalized-kummer-translation.asl");
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
fn example_distinguishes_geometric_inversion_and_its_frobenius_powers() {
    assert_eq!(
        stdout(example("hyperelliptic-generalized-kummer-inversion")),
        concat!(
            "rvGenKummerInversionEuler[2, 1] = 24\n",
            "rvGenKummerInversionEuler[2, 2] = 8\n",
            "rvGenKummerInversionEuler[3, 3] = 60\n",
            "rvGenKummerInversionHodge[3, 3, 2, 1] = -4\n",
            "rvGenKummerInversionPoints[B, C, 1] = 1191, affine torsion fixed points=0\n",
            "rvGenKummerInversionPoints[B, C, 1] = 1416, affine torsion fixed points=9\n",
            "rvGenKummerInversionPoints[B, C, 2] = 543051, affine torsion fixed points=0\n",
            "rvGenKummerInversionPoints[B, C, 2] = 593676, affine torsion fixed points=81\n",
            "rvGenKummerInversionFrob[B, S, 1, 3] = 1000\n",
            "rvGenKummerInversionPoints[B, S, 1] = 491676, affine torsion fixed points=1\n",
        )
    );
}

#[test]
fn hodge_characters_give_integral_involution_eigenspaces_and_the_fixed_euler_number() {
    for n in 1..=6 {
        let dimension = 2 * (n - 1);
        let mut source = String::new();
        for p in 0..=dimension {
            for q in 0..=dimension {
                writeln!(source, "rvGenKummerHodge[{n},{p},{q}];").unwrap();
            }
        }
        for h in (1..=n).filter(|h| n % h == 0) {
            writeln!(source, "rvGenKummerInversionEuler[{n},{h}];").unwrap();
            for p in 0..=dimension {
                for q in 0..=dimension {
                    writeln!(source, "rvGenKummerInversionHodge[{n},{h},{p},{q}];").unwrap();
                    if n == 2 {
                        writeln!(source, "rvGenKummerTranslateHodge[{n},{h},{p},{q}];").unwrap();
                    }
                }
            }
        }
        let output = stdout(interpret(&source));
        let mut values = output.lines().map(number);
        let width = dimension + 1;
        let dimensions: Vec<_> = values.by_ref().take(width * width).collect();
        let coefficients = inverted_product(n);
        for h in (1..=n).filter(|h| n % h == 0) {
            let euler = values.next().unwrap();
            let trace: Vec<_> = (0..dimensions.len())
                .map(|_| {
                    let value = values.next().unwrap();
                    if n == 2 {
                        assert_eq!(value, values.next().unwrap());
                    }
                    value
                })
                .collect();
            let mut lefschetz = 0;
            for p in 0..=dimension {
                for q in 0..=dimension {
                    let index = p * width + q;
                    let value = trace[index];
                    let rank = dimensions[index];
                    assert!(value.abs() <= rank);
                    assert_eq!((rank - value) % 2, 0);
                    assert_eq!(value, trace[q * width + p]);
                    assert_eq!(value, trace[(dimension - p) * width + dimension - q]);
                    lefschetz += if (p + q) % 2 == 0 { value } else { -value };
                }
            }
            let correction = if n % 2 == 0 {
                if n / h % 2 == 0 {
                    15
                } else {
                    -1
                }
            } else {
                0
            };
            assert_eq!(
                (coefficients[n] + correction * coefficients[n / 2]) / 16,
                euler
            );
            assert_eq!(lefschetz, euler);
            match (n, h) {
                (1, 1) => assert_eq!(euler, 1),
                (2, 1) => assert_eq!(euler, 24),
                (2, 2) => assert_eq!(euler, 8),
                (3, _) => {
                    assert_eq!(euler, 60);
                    assert_eq!(trace[2 * width + 1], -4);
                    assert_eq!(trace[2 * width + 2], 16);
                }
                (4, 1 | 2) => assert_eq!(euler, 464),
                (4, 4) => assert_eq!(euler, 320),
                _ => (),
            }
        }
        assert!(values.next().is_none());
    }
}

#[test]
fn length_two_agrees_with_actual_resolved_kummer_translation_and_frobenius() {
    let output = stdout(interpret(&format!(
        r"{EVEN}
        let frob=\x:RVKummerK3(D)=>rvKummerK3Frob[C,1](x);
        let moving=\x:RVKummerK3(D)=>rvKummerK3Translate[B,8](frob(x));
        assert (\x:RVKummerK3(D)=>moving(moving(x)))==(\x:RVKummerK3(D)=>rvKummerK3Translate[B,2](x));
        show (\x:RVKummerK3(D)=>rvKummerK3Translate[B,2](frob(x)));
        show (\x:RVKummerK3(D)=>rvKummerK3Translate[B,1](frob(x)));
        rvGenKummerInversionPoints[B,C,1](e2);
        rvGenKummerInversionPoints[B,C,1](e1);
        rvGenKummerInversionPoints[B,C,2](f2);
        rvGenKummerTranslatePoints[B,C,2](e2);
        rvGenKummerInversionPoints[B,C,2](e2);
        rvGenKummerLevel[B,C,2];
        "
    )));
    let mut lines = output.lines().skip(7);
    let direct: Vec<_> = lines.by_ref().take(2).map(fixed_in_map).collect();
    let reported: Vec<_> = lines.map(number).collect();
    assert_eq!(direct, reported[..2]);
    assert_eq!(reported, [32, 8, 104, 104, 248, 248]);
}

#[test]
fn moving_shifts_can_be_solvable_and_insoluble_under_order_three_frobenius() {
    let mut source = format!(
        r"{CUBIC}
        assert rvJacFrob[C,1](e2) == f2;
        assert rvJacFrob[C,3](e2) == e2;
        let soluble = add(e2,f2);
        let insoluble = add(soluble,e1);
        assert rvJacFrob[C,1](soluble) == add(f2,rvJacFrob[C,2](e2));
        "
    );
    for r in 1..=3 {
        writeln!(source, "rvGenKummerLevel[B,C,{r}];").unwrap();
        writeln!(source, "rvGenKummerInversionPoints[B,C,{r}](soluble);").unwrap();
        writeln!(source, "rvGenKummerInversionPoints[B,C,{r}](insoluble);").unwrap();
    }
    let output = stdout(interpret(&source));
    let lines: Vec<_> = output.lines().collect();
    for (index, triplet) in lines.chunks_exact(3).enumerate() {
        let q = 3i128.pow(index as u32 + 1);
        let kernel = if index == 2 { 16 } else { 4 };
        // The invariant e1 component obstructs odd powers. It cancels
        // from the square, which retains only a shift in the image.
        let other = if index == 1 { kernel } else { 0 };
        assert_eq!(number(triplet[0]), number(triplet[1]));
        assert_eq!(
            number(triplet[1]) - number(triplet[2]),
            q * (kernel - other)
        );
        assert!(triplet[1].ends_with(&format!("affine torsion fixed points={kernel}")));
        assert!(triplet[2].ends_with(&format!("affine torsion fixed points={other}")));
    }
    assert_eq!(lines.len(), 9);
}

#[test]
fn marked_odd_level_has_the_inversion_sign_and_a_single_middle_correction() {
    let cases = [
        ("C", 1, "e1", 5i128, 9, 0),
        ("C", 1, "fixed", 5, 9, 9),
        ("C", 2, "e1", 25, 81, 0),
        ("C", 2, "fixed", 25, 81, 81),
        ("S", 1, "e1", 25, 81, 1),
        ("S", 2, "e1", 625, 81, 81),
    ];
    let mut source = ODD.to_owned();
    for (base, r, t, _, _, _) in cases {
        for i in 0..=8 {
            writeln!(source, "rvGenKummerTrace[B,{base},{r},{i}];").unwrap();
            writeln!(source, "rvGenKummerInversionFrob[B,{base},{r},{i}]({t});").unwrap();
        }
        writeln!(source, "rvGenKummerInversionPoints[B,{base},{r}]({t});").unwrap();
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines().skip(10);
    for (base, r, _, q, ordinary_torsion, affine_torsion) in cases {
        let mut lefschetz = 0;
        for i in 0..=8 {
            let ordinary = number(lines.next().unwrap());
            let inverted = number(lines.next().unwrap());
            let expected = if i == 4 {
                ordinary + q.pow(2) * (affine_torsion - ordinary_torsion)
            } else if r * i % 2 == 0 {
                ordinary
            } else {
                -ordinary
            };
            assert_eq!(inverted, expected);
            if base == "S" && r == 1 && i == 3 {
                assert_eq!(inverted, 1000);
            }
            lefschetz += if i % 2 == 0 { inverted } else { -inverted };
        }
        let report = lines.next().unwrap();
        assert_eq!(number(report), lefschetz);
        assert!(report.ends_with(&format!("affine torsion fixed points={affine_torsion}")));
    }
    assert!(lines.next().is_none());
}

#[test]
fn validates_markings_orders_and_wide_ranges() {
    let max = u128::MAX;
    let output = stdout(interpret(&format!(
        "{ODD}
         rvGenKummerInversionHodge[{max},{max},1,2];
         rvGenKummerInversionFrob[B,C,{max},0](e1);
         rvGenKummerInversionFrob[B,C,{max},9](e1);"
    )));
    assert_eq!(
        output.lines().skip(10).map(number).collect::<Vec<_>>(),
        [-4, 1, 0]
    );
    for (setup, statement, expected) in [
        (
            "",
            "rvGenKummerInversionEuler[0,1];",
            "length must be positive",
        ),
        ("", "rvGenKummerInversionEuler[4,0];", "positive divisor"),
        (
            "",
            "rvGenKummerInversionHodge[4,3,100,100];",
            "positive divisor",
        ),
        (
            "",
            "rvGenKummerInversionEuler[1000000,1];",
            "coefficient storage limit",
        ),
        (
            EVEN,
            "rvGenKummerInversionPoints[B,C,0](e2);",
            "extension degree must be positive",
        ),
        (
            EVEN,
            "rvGenKummerInversionPoints[B,C,1](rvJacPoint[C](z,z));",
            "expected type RVJac(D)",
        ),
        (
            ODD,
            "rvGenKummerInversionFrob[B,C,1,9](rvJacPoint[S](rvUnramElem[W]{},rvUnramElem[W]{}));",
            "requires a class killed by 3",
        ),
        (
            ODD,
            "rvGenKummerInversionPoints[B,C,1000](e1);",
            "exact intermediate exceeds 512 bits",
        ),
        (
            EVEN,
            "rvhyperelliptic X[F]{z,o,z,z,z,o}; rvGenKummerInversionFrob[B,X,1,5](e2);",
            "not the declared base extension",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn inverted_product(n: usize) -> Vec<i128> {
    let mut result = vec![0; n + 1];
    result[0] = 1;
    for j in 1..=n {
        for _ in 0..8 {
            for i in j..=n {
                result[i] += result[i - j];
            }
            for i in (j..=n).rev() {
                result[i] += result[i - j];
            }
        }
    }
    result
}

fn fixed_in_map(map: &str) -> i128 {
    let points: Vec<_> = map
        .match_indices("RVKummerK3(D){")
        .map(|(start, _)| {
            let end = start + map[start..].find('}').unwrap() + 1;
            &map[start..end]
        })
        .collect();
    assert_eq!(points.len(), 2 * 248);
    points
        .chunks_exact(2)
        .filter(|pair| pair[0] == pair[1])
        .count() as i128
}

fn number(report: &str) -> i128 {
    report
        .split_once(" = ")
        .unwrap()
        .1
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
