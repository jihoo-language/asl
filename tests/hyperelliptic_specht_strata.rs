mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const FIXTURE: &str = include_str!("../examples/hyperelliptic-specht-strata.asl");
const INITIAL: usize = 6;

#[test]
fn example_reports_intrinsic_specht_traces_and_an_euler_coefficient() {
    assert_eq!(
        stdout(example("hyperelliptic-specht-strata")),
        concat!(
            "rvK3SpechtFrobTrace[Moving, C, Symmetric, 1, 4] = 22\n",
            "rvK3SpechtFrobTrace[Moving, C, Alternating, 1, 4] = 0\n",
            "rvK3SpechtFrobTrace[Moving, C, Standard, 1, 4] = 13\n",
            "rvK3SpechtEulerFactor[Moving, C, Standard, 1, 2, 2] = -9\n",
            "rvK3SpechtFrobTrace[Moving, C, Row, 1, 2] = -2\n",
            "rvK3StratumFrobTrace[Moving, C, Weighted, 1, 6] = -18\n",
        ),
    );
}

#[test]
fn specht_multiplicity_spaces_recombine_to_the_ordered_powers() {
    let degrees: [&[(&[usize], i128)]; 3] = [
        &[(&[2], 1), (&[1, 1], 1)],
        &[(&[3], 1), (&[2, 1], 2), (&[1, 1, 1], 1)],
        &[
            (&[4], 1),
            (&[3, 1], 3),
            (&[2, 2], 2),
            (&[2, 1, 1], 3),
            (&[1, 1, 1, 1], 1),
        ],
    ];
    for (n, diagrams) in (2..=4).zip(degrees) {
        let mut source = format!("{FIXTURE} rvhilbstratum Ordered{{(1,{n})}};");
        for (index, (diagram, _)) in diagrams.iter().enumerate() {
            let rows: Vec<_> = diagram.iter().map(|r| r.to_string()).collect();
            writeln!(
                source,
                "rvspechtstratum D{index}[Ordered]{{(1,[{}])}};",
                rows.join(",")
            )
            .unwrap();
        }
        for r in 1..=2 {
            for exponent in [r, 3 * r] {
                for i in (0..=4).step_by(2) {
                    writeln!(source, "rvK3FrobHilbTrace[Moving,C,1,{exponent},{i}];").unwrap();
                }
            }
            for i in 0..=4 * n {
                for index in 0..diagrams.len() {
                    writeln!(source, "rvK3SpechtFrobTrace[Moving,C,D{index},{r},{i}];").unwrap();
                }
            }
        }
        let output = stdout(interpret(&source));
        let mut values = output.lines().skip(INITIAL).map(number);
        for _ in 1..=2 {
            let surface: Vec<_> = values.by_ref().take(3).collect();
            let third: Vec<_> = values.by_ref().take(3).collect();
            let ordered = (0..n).fold(vec![1], |product, _| product_of(&product, &surface));
            for i in 0..=4 * n {
                let pieces: Vec<_> = values.by_ref().take(diagrams.len()).collect();
                let recombined = pieces
                    .iter()
                    .zip(diagrams)
                    .map(|(&trace, &(_, dimension))| trace * dimension)
                    .sum::<i128>();
                let total = if i % 2 == 0 { ordered[i / 2] } else { 0 };
                assert_eq!(recombined, total, "ordered {n} in degree {i}");
                if n == 3 {
                    let adams = if i % 6 == 0 { third[i / 6] } else { 0 };
                    assert_eq!(3 * pieces[1], total - adams, "S_(2,1) in degree {i}");
                }
            }
        }
        assert!(values.next().is_none());
    }
}

#[test]
fn independent_size_groups_tensor_their_schur_functors() {
    let mut source = format!(
        "{FIXTURE} rvhilbstratum Mixed{{(1,2),(2,3)}};
         rvspechtstratum Both[Mixed]{{(2,[2,1]),(1,[1,1])}};"
    );
    for r in 1..=2 {
        for (name, dimension) in [("Alternating", 4), ("Standard", 6), ("Both", 10)] {
            for half in 0..=dimension {
                writeln!(
                    source,
                    "rvK3SpechtFrobTrace[Moving,C,{name},{r},{}];",
                    2 * half
                )
                .unwrap();
            }
        }
    }
    let output = stdout(interpret(&source));
    let mut values = output.lines().skip(INITIAL).map(number);
    for _ in 1..=2 {
        let exterior: Vec<_> = values.by_ref().take(5).collect();
        let standard: Vec<_> = values.by_ref().take(7).collect();
        let both: Vec<_> = values.by_ref().take(11).collect();
        assert_eq!(both, product_of(&exterior, &standard));
    }
    assert!(values.next().is_none());
}

#[test]
fn trivial_decorations_recover_hilbert_factors_after_the_tate_shift() {
    let mut source = FIXTURE.to_owned();
    for r in 1..=2 {
        for i in 0..=8 {
            writeln!(source, "rvK3SpechtFrobTrace[Moving,C,Row,{r},{i}];").unwrap();
            writeln!(
                source,
                "rvK3StratumFrobTrace[Moving,C,Weighted,{r},{}];",
                i + 4
            )
            .unwrap();
            for j in 0..=3 {
                writeln!(source, "rvK3SpechtEulerFactor[Moving,C,Row,{r},{i},{j}];").unwrap();
                writeln!(
                    source,
                    "rvK3StratumEulerFactor[Moving,C,Weighted,{r},{},{j}];",
                    i + 4
                )
                .unwrap();
            }
        }
    }
    let output = stdout(interpret(&source));
    let mut values = output.lines().skip(INITIAL).map(number);
    for r in 1..=2 {
        for _ in 0..=8 {
            let intrinsic = values.next().unwrap();
            assert_eq!(values.next().unwrap(), 3i128.pow(2 * r) * intrinsic);
            for j in 0..=3 {
                let intrinsic = values.next().unwrap();
                assert_eq!(values.next().unwrap(), 3i128.pow(2 * r * j) * intrinsic);
            }
        }
    }
    assert!(values.next().is_none());
}

#[test]
fn row_bounds_wide_shapes_and_formal_factors_are_exact() {
    let max = u128::MAX;
    let bottom = max - 1;
    let tall = ["1"; 25].join(",");
    let full = ["1"; 24].join(",");
    let output = stdout(interpret(&format!(
        "{FIXTURE}
         rvhilbstratum None{{}}; rvspechtstratum Empty[None]{{}};
         rvhilbstratum Wide{{(1,{max})}}; rvspechtstratum Hook[Wide]{{(1,[{bottom},1])}};
         rvhilbstratum Full{{(1,24)}}; rvspechtstratum Exterior[Full]{{(1,[{full}])}};
         rvhilbstratum Tall{{(1,25)}}; rvspechtstratum Vanishes[Tall]{{(1,[{tall}])}};
         rvK3SpechtFrobTrace[Moving,C,Empty,{max},0];
         rvK3SpechtEulerFactor[Moving,C,Empty,{max},0,1];
         rvK3SpechtEulerFactor[Moving,C,Empty,{max},0,{max}];
         rvK3SpechtFrobTrace[Moving,C,Hook,{max},0];
         rvK3SpechtFrobTrace[Moving,C,Hook,1,2];
         rvK3SpechtFrobTrace[Moving,C,Exterior,1,48];
         rvK3SpechtFrobTrace[Moving,C,Exterior,{max},46];
         rvK3SpechtEulerFactor[Moving,C,Exterior,{max},48,2];
         rvK3SpechtFrobTrace[Moving,C,Vanishes,{max},50];
         rvK3SpechtEulerFactor[Moving,C,Vanishes,{max},50,0];
         rvK3SpechtEulerFactor[Moving,C,Vanishes,{max},50,1];
         rvK3SpechtEulerFactor[Moving,C,Alternating,{max},2,23];"
    )));
    assert_eq!(
        output.lines().skip(INITIAL).map(number).collect::<Vec<_>>(),
        [1, -1, 0, 0, -2, 3i128.pow(24), 0, 0, 0, 1, 0, 0],
    );
}

#[test]
fn rejects_malformed_diagrams_mismatched_names_and_arithmetic_bounds() {
    for (statement, expected) in [
        ("rvspechtstratum Bad[Square]{(1,[1])};", "must equal the part multiplicity"),
        ("rvspechtstratum Bad[Square]{(1,[])};", "rows must be positive"),
        ("rvspechtstratum Bad[Square]{(1,[0,2])};", "rows must be positive"),
        ("rvspechtstratum Bad[Triple]{(2,[1,2])};", "weakly decreasing"),
        ("rvspechtstratum Bad[Square]{(2,[2])};", "size must occur"),
        ("rvspechtstratum Bad[Square]{(1,[2]),(1,[1,1])};", "already specified"),
        ("rvspechtstratum Row[Square]{};", "already defined"),
        ("rvspechtstratum Bad[Unlisted]{};", "unknown Hilbert–Chow stratum"),
        ("rvK3SpechtEulerFactor[Moving,C,Square,1,0,0];", "unknown Specht stratum"),
        ("rvK3SpechtFrobTrace[Moving,C,Standard,0,2];", "iteration count must be positive"),
        ("rvK3SpechtFrobTrace[Moving,C,Standard,1000,2];", "exact intermediate exceeds 512 bits"),
        ("rvhyperelliptic Wrong[F]{o,neg(o),z,z,z,o}; rvK3SpechtEulerFactor[Moving,Wrong,Row,1,0,0];", "not the declared base extension"),
    ] {
        let error = diagnostic(interpret(&format!("{FIXTURE}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn product_of(left: &[i128], right: &[i128]) -> Vec<i128> {
    let mut product = vec![0; left.len() + right.len() - 1];
    for (i, &a) in left.iter().enumerate() {
        for (j, &b) in right.iter().enumerate() {
            product[i + j] += a * b;
        }
    }
    product
}

fn number(line: &str) -> i128 {
    line.split_once(" = ").unwrap().1.parse().unwrap()
}
