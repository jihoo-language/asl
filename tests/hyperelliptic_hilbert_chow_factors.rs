mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const FIXTURE: &str = include_str!("../examples/hyperelliptic-hilbert-chow-factors.asl");
const INITIAL: usize = 10;

#[test]
fn example_reports_a_rank_one_factor_and_noncommuting_middle_traces() {
    assert_eq!(
        stdout(example("hyperelliptic-hilbert-chow-factors")),
        concat!(
            "rvK3StratumFrobTrace[Moving, C, Diagonal, 1, 4] = -6\n",
            "rvK3StratumFrobTrace[Moving, C, Diagonal, 2, 4] = 198\n",
            "rvK3StratumFrobTrace[Moving, C, Distinct, 1, 4] = 22\n",
            "rvK3FrobHilbTrace[Moving, C, 2, 1, 4] = 16\n",
            "rvK3StratumEulerFactor[Moving, C, Diagonal, 1, 2, 0] = 1\n",
            "rvK3StratumEulerFactor[Moving, C, Diagonal, 1, 2, 1] = -3\n",
            "rvK3StratumEulerFactor[Moving, C, Diagonal, 1, 2, 2] = 0\n",
            "rvK3StratumEulerFactor[Moving, C, Diagonal, 1, 4, 1] = 6\n",
            "rvK3StratumEulerFactor[Moving, C, Diagonal, 1, 4, 2] = -81\n",
            "rvK3StratumEulerFactor[Idle, C, Diagonal, 1, 4, 2] = 1431\n",
        ),
    );
}

#[test]
fn all_partition_traces_reconstruct_hilbert_cohomology() {
    for n in 0..=4 {
        let shapes = partitions(n, 1);
        let mut source =
            format!("{FIXTURE} rvk3autom Hyperbolic[B]{{project(0),project(1),project(2)}};");
        for (index, shape) in shapes.iter().enumerate() {
            let parts: Vec<_> = shape.iter().map(|size| format!("({size},1)")).collect();
            writeln!(source, "rvhilbstratum P{index}{{{}}};", parts.join(",")).unwrap();
        }
        for (word, base, r) in [
            ("Moving", "C", 1),
            ("Moving", "C", 2),
            ("Hyperbolic", "C", 1),
            ("Idle", "L", 1),
        ] {
            for i in 0..=4 * n + 1 {
                writeln!(source, "rvK3FrobHilbTrace[{word},{base},{n},{r},{i}];").unwrap();
                for index in 0..shapes.len() {
                    writeln!(
                        source,
                        "rvK3StratumFrobTrace[{word},{base},P{index},{r},{i}];"
                    )
                    .unwrap();
                }
            }
        }
        let output = stdout(interpret(&source));
        let values: Vec<_> = output.lines().skip(INITIAL).map(number).collect();
        assert_eq!(values.len(), 4 * (4 * n + 2) * (shapes.len() + 1));
        for group in values.chunks_exact(shapes.len() + 1) {
            assert_eq!(group[0], group[1..].iter().sum::<i128>(), "at weight {n}");
        }
    }
}

#[test]
fn stratum_determinants_multiply_to_the_hilbert_square_factor() {
    const STOP: usize = 6;
    let actions = [
        ("Moving", 1),
        ("Moving", 2),
        ("Idle", 1),
        ("Idle", 2),
        ("Hyperbolic", 1),
    ];
    let mut source =
        format!("{FIXTURE} rvk3autom Hyperbolic[B]{{project(0),project(1),project(2)}};");
    for (word, r) in actions {
        for i in (0..=8).step_by(2) {
            for k in 1..=STOP {
                writeln!(source, "rvK3FrobHilbTrace[{word},C,2,{},{i}];", r * k).unwrap();
            }
            for shape in ["Diagonal", "Distinct"] {
                for j in 0..=STOP {
                    writeln!(
                        source,
                        "rvK3StratumEulerFactor[{word},C,{shape},{r},{i},{j}];"
                    )
                    .unwrap();
                }
            }
            if word == "Idle" {
                writeln!(source, "rvKummerHilbChar[C,2,{r},{i},101,{STOP}];").unwrap();
            }
        }
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines().skip(INITIAL);
    for (word, _) in actions {
        for _ in 0..5 {
            let traces: Vec<_> = lines.by_ref().take(STOP).map(number).collect();
            let diagonal: Vec<_> = lines.by_ref().take(STOP + 1).map(number).collect();
            let distinct: Vec<_> = lines.by_ref().take(STOP + 1).map(number).collect();
            let product: Vec<_> = (0..=STOP)
                .map(|j| (0..=j).map(|k| diagonal[k] * distinct[j - k]).sum::<i128>())
                .collect();
            assert_eq!(product, newton(&traces));
            if word == "Idle" {
                let residues: Vec<i128> = rhs(lines.next().unwrap())
                    .strip_prefix('[')
                    .unwrap()
                    .strip_suffix(']')
                    .unwrap()
                    .split(", ")
                    .map(|value| value.parse().unwrap())
                    .collect();
                assert_eq!(
                    product
                        .iter()
                        .map(|value| value.rem_euclid(101))
                        .collect::<Vec<_>>(),
                    residues
                );
            }
        }
    }
    assert!(lines.next().is_none());
}

#[test]
fn diagonal_middle_factor_has_the_actual_surface_rank_and_tate_weight() {
    let output = stdout(interpret(&format!(
        "{FIXTURE}
         rvK3StratumEulerFactor[Moving,C,Diagonal,1,4,22];
         rvK3StratumEulerFactor[Moving,C,Diagonal,1,4,23];
         rvK3StratumFrobTrace[Idle,C,Diagonal,1,4];
         rvKummerHilbTrace[C,1,1,2];
         rvK3StratumFrobTrace[Idle,C,Diagonal,2,6];
         rvK3StratumEulerFactor[Idle,C,Diagonal,2,6,1];"
    )));
    let values: Vec<_> = output.lines().skip(INITIAL).map(number).collect();
    assert_eq!(values[0], 3i128.pow(44));
    assert_eq!(values[1], 0);
    assert_eq!(values[2], 3 * values[3]);
    assert_eq!(values[4..], [729, -729]);
}

#[test]
fn elementary_factors_and_unused_huge_inputs_avoid_power_work() {
    let max = u128::MAX;
    let output = stdout(interpret(&format!(
        "{FIXTURE}
         rvhilbstratum Empty{{}};
         rvhilbstratum Ones{{(1,{max})}};
         rvhilbstratum BigPart{{({max},1)}};
         rvK3StratumEulerFactor[Moving,C,Empty,{max},0,0];
         rvK3StratumEulerFactor[Moving,C,Empty,{max},0,1];
         rvK3StratumEulerFactor[Moving,C,Empty,{max},0,{max}];
         rvK3StratumFrobTrace[Moving,C,Empty,{max},0];
         rvK3StratumFrobTrace[Moving,C,Empty,{max},2];
         rvK3StratumFrobTrace[Moving,C,Ones,{max},0];
         rvK3StratumFrobTrace[Moving,C,Ones,1,2];
         rvK3StratumEulerFactor[Moving,C,Diagonal,{max},4,23];
         rvK3StratumEulerFactor[Moving,C,Diagonal,{max},2,2];
         rvK3StratumEulerFactor[Moving,C,Ones,{max},5,0];
         rvK3StratumEulerFactor[Moving,C,Ones,{max},5,1];
         rvK3StratumFrobTrace[Moving,C,BigPart,{max},{max}];"
    )));
    assert_eq!(
        output.lines().skip(INITIAL).map(number).collect::<Vec<_>>(),
        [1, -1, 0, 1, 0, 1, -2, 0, 0, 1, 0, 0]
    );
    for (statement, expected) in [
        ("rvK3StratumFrobTrace[Moving,C,Diagonal,0,2];", "iteration count must be positive"),
        ("rvK3StratumEulerFactor[Moving,C,Unlisted,1,0,0];", "unknown Hilbert–Chow stratum"),
        ("rvhyperelliptic Wrong[F]{o,neg(o),z,z,z,o}; rvK3StratumEulerFactor[Moving,Wrong,Diagonal,1,0,0];", "not the declared base extension"),
        ("rvK3StratumEulerFactor[Moving,C,Diagonal,1000,4,1];", "exact intermediate exceeds 512 bits"),
        ("rvK3StratumEulerFactor[Moving,C,Diagonal,40,6,1];", "exceeds signed 128-bit range"),
        ("rvhilbstratum Large{(1,1000000)}; rvK3StratumEulerFactor[Moving,C,Large,1,10000,2];", "coefficient storage limit"),
    ] {
        let error = diagnostic(interpret(&format!("{FIXTURE}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn partitions(n: usize, minimum: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![Vec::new()];
    }
    let mut result = Vec::new();
    for size in minimum..=n {
        for mut tail in partitions(n - size, size) {
            tail.insert(0, size);
            result.push(tail);
        }
    }
    result
}

fn newton(traces: &[i128]) -> Vec<i128> {
    let mut coefficients = vec![1];
    for j in 1..=traces.len() {
        let sum: i128 = (1..=j).map(|k| traces[k - 1] * coefficients[j - k]).sum();
        assert_eq!(sum % j as i128, 0);
        coefficients.push(-sum / j as i128);
    }
    coefficients
}

fn number(line: &str) -> i128 {
    rhs(line).parse().unwrap()
}

fn rhs(line: &str) -> &str {
    line.split_once(" = ").unwrap().1
}
