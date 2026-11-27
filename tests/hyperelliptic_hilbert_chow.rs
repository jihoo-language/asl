mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const FIXTURE: &str = include_str!("../examples/hyperelliptic-hilbert-chow-strata.asl");
const ODD_LEVEL: &str = include_str!("../examples/hyperelliptic-generalized-kummer-cohomology.asl");

#[test]
fn example_reports_signed_local_constants_and_compact_component_counts() {
    assert_eq!(
        stdout(example("hyperelliptic-hilbert-chow-strata")),
        concat!(
            "rvHilbStratumData[Double] = n=2, L=1, e=2, fiber dimension=1, codimension=2, intersection multiplier=-2\n",
            "rvHilbStratumData[TwoDoubles] = n=4, L=2, e=2, fiber dimension=2, codimension=4, intersection multiplier=4\n",
            "rvK3StratumHodge[Double, 1, 1] = 1\n",
            "rvJacStratumHodge[Double, 2, 1] = 2\n",
            "rvGenKummerStratumHodge[Double, 1, 1, 1] = 16\n",
            "rvGenKummerStratumHodge[Double, 2, 1, 1] = 0\n",
            "rvGenKummerStratumHodge[TwoDoubles, 2, 2, 2] = 16\n",
            "rvGenKummerStratumHodge[TwoDoubles, 4, 2, 2] = 0\n",
            "rvNakajimaPair = -2\n",
            "rvGenKummerStratumFrob[B, C, Double, 1, 2] = 24\n",
            "rvGenKummerStratumFrob[B, C, Double, 1, 2] = 0\n",
        ),
    );
}

#[test]
fn all_partitions_through_five_reconstruct_hilbert_and_kummer_hodge_characters() {
    for n in 0..=5 {
        let shapes = partitions(n, 1);
        let mut source = String::new();
        for (index, shape) in shapes.iter().enumerate() {
            let parts: Vec<_> = shape.iter().map(|size| format!("({size},1)")).collect();
            writeln!(source, "rvhilbstratum P{index}{{{}}};", parts.join(",")).unwrap();
        }
        let mut comparisons = 0;
        for p in 0..=2 * n {
            for q in 0..=2 * n {
                for (total, summand) in [
                    ("rvK3HilbHodge", "rvK3StratumHodge"),
                    ("rvJacHilbHodge", "rvJacStratumHodge"),
                ] {
                    writeln!(source, "{total}[{n},{p},{q}];").unwrap();
                    for index in 0..shapes.len() {
                        writeln!(source, "{summand}[P{index},{p},{q}];").unwrap();
                    }
                    comparisons += 1;
                }
                for h in (1..=n).filter(|&h| n % h == 0) {
                    writeln!(source, "rvGenKummerTranslateHodge[{n},{h},{p},{q}];").unwrap();
                    for index in 0..shapes.len() {
                        writeln!(source, "rvGenKummerStratumHodge[P{index},{h},{p},{q}];").unwrap();
                    }
                    comparisons += 1;
                }
            }
        }
        let output = stdout(interpret(&source));
        let values: Vec<_> = output.lines().map(number).collect();
        assert_eq!(values.len(), comparisons * (shapes.len() + 1));
        for group in values.chunks_exact(shapes.len() + 1) {
            assert_eq!(group[0], group[1..].iter().sum::<i128>(), "at length {n}");
        }
    }
}

#[test]
fn projectors_preserve_grades_and_are_orthogonal_for_the_wick_pairing() {
    let output = stdout(interpret(&format!(
        r"{FIXTURE}
        rvhilbstratum Open{{(1,2)}};
        rvhilbstratum Empty{{}};
        rvhilbstratum Huge{{(1,1000000000000000000000000000000)}};
        let mix = add(rvNakajimaScale[3](diagonal),rvNakajimaScale[-2](open));
        let diagonalPart = rvNakajimaStratum[Double](mix);
        let openPart = rvNakajimaStratum[Open](mix);
        assert diagonalPart == rvNakajimaScale[3](diagonal);
        assert openPart == rvNakajimaScale[-2](open);
        assert add(diagonalPart,openPart) == mix;
        assert rvNakajimaStratum[Double](diagonalPart) == diagonalPart;
        assert rvNakajimaStratum[Open](openPart) == openPart;
        assert rvNakajimaStratum[Double](openPart) == zero[RVNakajima(B,2,2)];
        assert rvNakajimaStratum[Empty](diagonalPart) == zero[RVNakajima(B,2,2)];
        assert rvNakajimaStratum[Empty](rvNakajimaVac[B]) == rvNakajimaVac[B];
        assert rvNakajimaStratum[Empty](zero[RVNakajima(B,-1,3)]) == zero[RVNakajima(B,-1,3)];
        let giant = rvNakajimaMonomial[B]{{(1,0,1000000000000000000000000000000)}};
        assert rvNakajimaStratum[Huge](giant) == giant;
        let dualDiagonal = rvNakajimaMonomial[B]{{(2,23,1)}};
        let dualOpen = rvNakajimaMonomial[B]{{(1,23,1),(1,1,1)}};
        let dual = add(dualDiagonal,dualOpen);
        rvNakajimaPair(diagonalPart,rvNakajimaStratum[Open](dual));
        rvNakajimaPair(openPart,rvNakajimaStratum[Double](dual));
        rvNakajimaPair(diagonalPart,dual);
        rvNakajimaPair(openPart,dual);
        rvNakajimaPair(mix,dual);
        "
    )));
    assert_eq!(
        output.lines().skip(11).map(number).collect::<Vec<_>>(),
        [0, 0, -6, 4, -2],
    );
}

#[test]
fn full_diagonal_and_divisible_shapes_use_the_actual_component_permutation() {
    let output = stdout(interpret(
        "rvhilbstratum Triple{(3,1)};
         rvhilbstratum TwoDoubles{(2,2)};
         rvGenKummerStratumHodge[Triple,1,2,2];
         rvGenKummerStratumHodge[Triple,3,2,2];
         rvGenKummerStratumHodge[Triple,1,2,3];
         rvGenKummerStratumHodge[TwoDoubles,1,2,2];
         rvGenKummerStratumHodge[TwoDoubles,2,2,2];
         rvGenKummerStratumHodge[TwoDoubles,4,2,2];",
    ));
    assert_eq!(
        output.lines().map(number).collect::<Vec<_>>(),
        [81, 0, 0, 16, 16, 0]
    );
}

#[test]
fn weighted_symmetric_power_frobenius_summands_reconstruct_the_total_trace() {
    for (setup, skipped, n, curve, moving) in [
        (FIXTURE, 11, 2, "L", "moving"),
        (ODD_LEVEL, 9, 3, "S", "e1"),
    ] {
        let shapes = partitions(n, 1);
        let mut source = setup.to_owned();
        for (index, shape) in shapes.iter().enumerate() {
            let parts: Vec<_> = shape.iter().map(|size| format!("({size},1)")).collect();
            writeln!(source, "rvhilbstratum P{index}{{{}}};", parts.join(",")).unwrap();
        }
        let zero = format!("zero[RVJac({curve})]");
        for r in 1..=2 {
            for t in [zero.as_str(), moving] {
                for i in 0..=4 * (n - 1) + 1 {
                    writeln!(source, "rvGenKummerTranslateFrob[B,C,{r},{i}]({t});").unwrap();
                    for index in 0..shapes.len() {
                        writeln!(source, "rvGenKummerStratumFrob[B,C,P{index},{r},{i}]({t});")
                            .unwrap();
                    }
                }
            }
        }
        let output = stdout(interpret(&source));
        let values: Vec<_> = output.lines().skip(skipped).map(number).collect();
        let width = 4 * (n - 1) + 2;
        let stride = shapes.len() + 1;
        assert_eq!(values.len(), 4 * width * stride);
        for (index, group) in values.chunks_exact(stride).enumerate() {
            let division = index / width;
            let degree = index % width;
            assert_eq!(group[0], group[1..].iter().sum::<i128>());
            // The last partition generated is the full diagonal. The others
            // contribute equally with and without the translation.
            if division % 2 == 1 {
                let untwisted = (index - width) * stride;
                assert_eq!(
                    &group[1..stride - 1],
                    &values[untwisted + 1..untwisted + stride - 1]
                );
            }
            let diagonal = *group.last().unwrap();
            if degree != 2 * (n - 1) {
                assert_eq!(diagonal, 0);
            } else {
                let expected = if n == 2 {
                    [24, 0, 144, 0]
                } else {
                    [225, 225, 50625, 0]
                };
                assert_eq!(diagonal, expected[division]);
            }
        }
    }
}

#[test]
fn huge_sparse_shapes_and_formal_ranges_are_checked() {
    let max = u128::MAX;
    let shift = max - 1;
    let output = stdout(interpret(&format!(
        "rvhilbstratum Empty{{(9,0)}};
         rvHilbStratumData[Empty];
         rvK3StratumHodge[Empty,0,0];
         rvJacStratumHodge[Empty,1,0];
         rvhilbstratum Max{{({max},1)}};
         rvHilbStratumData[Max];
         rvK3StratumHodge[Max,{shift},{shift}];
         rvGenKummerStratumHodge[Max,{max},{shift},{shift}];
         rvhilbstratum Doubles{{(2,1000000000000000000)}};
         rvHilbStratumData[Doubles];
         rvGenKummerStratumHodge[Doubles,2,1000000000000000000,1000000000000000000];
         rvhilbstratum Ones{{(1,{max})}};
         rvGenKummerStratumHodge[Ones,{max},1,1];"
    )));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines[0], "rvHilbStratumData[Empty] = n=0, L=0, e=0, fiber dimension=0, codimension=0, intersection multiplier=1");
    assert!(lines[3].contains("codimension=680564733841876926926749214863536422908"));
    assert!(lines[6].ends_with("intersection multiplier=2^1000000000000000000"));
    assert_eq!(
        lines
            .iter()
            .enumerate()
            .filter(|(i, _)| ![0, 3, 6].contains(i))
            .map(|(_, line)| number(line))
            .collect::<Vec<_>>(),
        [1, 0, 1, 0, 16, 4]
    );

    for (source, message) in [
        ("rvhilbstratum S{(0,1)};", "part size must be positive"),
        ("rvhilbstratum S{(340282366920938463463374607431768211455,2)};", "partition weight exceeds"),
        ("rvhilbstratum S{}; rvhilbstratum S{(1,1)};", "is already defined"),
        ("rvJacStratumHodge[S,0,0];", "unknown Hilbert–Chow stratum"),
        ("rvhilbstratum S{}; rvGenKummerStratumHodge[S,1,0,0];", "length must be positive"),
        ("rvhilbstratum S{(2,2)}; rvGenKummerStratumHodge[S,3,0,0];", "positive divisor"),
        ("rvhilbstratum S{(1,10)}; show rvNakajimaStratum[S](S.hi);", "a Nakajima state"),
        ("rvhilbstratum S{(1,1000000)}; rvK3StratumHodge[S,10000,10000];", "coefficient table exceeds the storage limit"),
        ("rvhilbstratum S{(340282366920938463463374607431768211455,1)}; rvGenKummerStratumHodge[S,1,340282366920938463463374607431768211454,340282366920938463463374607431768211454];", "exceeds unsigned 128-bit range"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(message), "{error}");
    }
    for (source, message) in [
        (
            "rvGenKummerStratumFrob[B,C,TwoDoubles,1,2](moving);",
            "weight must equal the prime level",
        ),
        (
            "rvGenKummerStratumFrob[B,C,Double,0,2](moving);",
            "extension degree must be positive",
        ),
        (
            "rvGenKummerStratumFrob[B,C,Double,1,2](rvJacPoint[C](z,z));",
            "expected type RVJac(L)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{FIXTURE}\n{source}")));
        assert!(error.contains(message), "{error}");
    }
}

fn partitions(n: usize, min: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![Vec::new()];
    }
    let mut result = Vec::new();
    for size in min..=n {
        for mut tail in partitions(n - size, size) {
            tail.insert(0, size);
            result.push(tail);
        }
    }
    result
}

fn number(line: &str) -> i128 {
    line.split_once(" = ").unwrap().1.parse().unwrap()
}
