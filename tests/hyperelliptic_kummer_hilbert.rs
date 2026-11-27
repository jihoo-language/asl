mod common;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    rvhyperelliptic C[F]{z,t,z,t,o,o};
";

const F25: &str = r"
    rvunram W[5,1]{2,0,1};
    let z=rvUnramElem[W]{};
    let o=rvUnramElem[W]{1};
    rvhyperelliptic S[W]{z,o,z,z,z,o};
";

fn integer_reports(out: &str) -> Vec<i128> {
    out.lines()
        .map(|line| line.split_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn runnable_example_uses_nonsplit_and_scalar_frobenius() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-hilbert")),
        concat!(
            "rvKummerHilbPoints[C, 2, 1] = 327\n",
            "rvKummerHilbPoints[C, 3, 1] = 4480\n",
            "rvKummerHilbTrace[C, 2, 1, 2] = 13\n",
            "rvKummerHilbTrace[C, 2, 1, 4] = 115\n",
            "rvKummerHilbTrace[C, 2, 1, 6] = 117\n",
            "rvKummerHilbBetti[2, 4] = 276\n",
            "rvKummerHilbBetti[3, 4] = 299\n",
            "rvKummerHilbTrace[S, 2, 1, 4] = 172500\n",
            "rvKummerHilbPoints[S, 2, 1] = 923076\n"
        )
    );
}

#[test]
fn point_counts_agree_with_cycle_partitions_and_the_trace_sum() {
    for r in [1, 2] {
        let mut source = format!(
            "{F3} rvKummerK3Points[C,{r}]; rvKummerK3Points[C,{}]; rvKummerK3Points[C,{}];",
            2 * r,
            3 * r
        );
        source.push_str(&format!(
            "rvKummerHilbPoints[C,2,{r}]; rvKummerHilbPoints[C,3,{r}];"
        ));
        for i in 0..=12 {
            source.push_str(&format!("rvKummerHilbTrace[C,3,{r},{i}];"));
        }
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        let [a, b, c]: [i128; 3] = std::array::from_fn(|_| {
            lines
                .next()
                .unwrap()
                .split_once("resolved=")
                .unwrap()
                .1
                .parse()
                .unwrap()
        });
        let q = 3i128.pow(r);
        let counts = integer_reports(&lines.collect::<Vec<_>>().join("\n"));
        assert_eq!(counts[0], (a * a + b) / 2 + q * a);
        assert_eq!(
            counts[1],
            (a.pow(3) + 3 * a * b + 2 * c) / 6 + q * a * a + q * q * a
        );
        let traces = &counts[2..];
        assert_eq!(traces.iter().sum::<i128>(), counts[1]);
        for i in 0..=6 {
            assert_eq!(traces[12 - i], q.pow(6 - i as u32) * traces[i]);
        }
        if r == 1 {
            assert_eq!(
                traces,
                [1, 0, 13, 0, 154, 0, 1144, 0, 1386, 0, 1053, 0, 729]
            );
        }
    }
}

#[test]
fn betti_numbers_and_euler_characteristics_match_the_k3_hilbert_scheme() {
    let mut source = String::new();
    for n in [2, 3] {
        for i in 0..=4 * n {
            source.push_str(&format!("rvKummerHilbBetti[{n},{i}];"));
        }
    }
    let betti = integer_reports(&stdout(interpret(&source)));
    let (square, cube) = betti.split_at(9);
    assert_eq!(square, [1, 0, 23, 0, 276, 0, 23, 0, 1]);
    assert_eq!(cube, [1, 0, 23, 0, 299, 0, 2554, 0, 299, 0, 23, 0, 1]);
    assert_eq!(square.iter().sum::<i128>(), 324);
    assert_eq!(cube.iter().sum::<i128>(), 3200);
}

#[test]
fn superspecial_frobenius_weights_every_piece_by_its_half_degree() {
    let mut source = F25.to_owned();
    for i in (0..=12).step_by(2) {
        source.push_str(&format!("rvKummerHilbTrace[S,3,1,{i}];"));
    }
    source.push_str("rvKummerHilbPoints[S,3,1];");
    let values = integer_reports(&stdout(interpret(&source)));
    let betti = [1, 23, 299, 2554, 299, 23, 1];
    for (d, (trace, dimension)) in values.iter().zip(betti).enumerate() {
        assert_eq!(*trace, dimension * 25i128.pow(d as u32));
    }
    assert_eq!(values[7], 625_640_576);
}

#[test]
fn surface_trace_can_be_negative_and_higher_degree_queries_are_independent() {
    let output = stdout(interpret(
        r"
        rvunram F[11,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let a=rvUnramElem[F]{4};
        rvhyperelliptic C[F]{o,z,o,a,o,o};
        rvKummerHilbTrace[C,1,1,2];
        rvKummerHilbTrace[C,2,1,2];
        rvKummerHilbTrace[C,2,1,6];
        rvKummerHilbPoints[C,1,1];
        ",
    ));
    assert_eq!(integer_reports(&output), [-4, 7, 847, 118]);

    let max = u128::MAX;
    let output = stdout(interpret(&format!(
        "{F3} rvKummerHilbPoints[C,0,{max}]; rvKummerHilbTrace[C,{max},{max},0];
         rvKummerHilbTrace[C,{max},{max},3]; rvKummerHilbTrace[C,{max},1,2];
         rvKummerHilbTrace[C,128,1,2]; rvKummerHilbTrace[C,2,{max},10];
         rvKummerHilbBetti[0,0]; rvKummerHilbBetti[0,1]; rvKummerHilbBetti[{max},2];"
    )));
    assert_eq!(integer_reports(&output), [1, 1, 0, 13, 13, 0, 1, 0, 23]);
    assert!(
        diagnostic(interpret(&format!("{F3} rvKummerHilbPoints[C,128,1];")))
            .contains("point count exceeds unsigned 128-bit range")
    );
}

#[test]
fn arithmetic_bounds_and_inputs_are_checked() {
    let output = stdout(interpret(
        r"
        rvunram F[19,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        rvhyperelliptic C[F]{z,o,z,z,z,o};
        rvKummerHilbPoints[C,1,15];
        ",
    ));
    assert_eq!(
        output,
        "rvKummerHilbPoints[C, 1, 15] = 230466617897195215130198816793320095000\n"
    );
    for (source, expected) in [
        (
            "rvKummerHilbPoints[C,0,0];",
            "extension degree must be positive",
        ),
        (
            "rvKummerHilbTrace[C,2,0,0];",
            "extension degree must be positive",
        ),
        (
            "rvKummerHilbPoints[C,-1,1];",
            "expected a nonnegative decimal integer",
        ),
        (
            "rvKummerHilbTrace[C,100,1,400];",
            "Frobenius trace exceeds signed 128-bit range",
        ),
        (
            "rvKummerHilbBetti[1000000000,2000000000];",
            "coefficient storage limit",
        ),
        (
            "rvhyperelliptic E[F]{o,o,z,o}; rvKummerHilbTrace[E,0,1,0];",
            "smooth monic quintic",
        ),
        (
            "rvKummerHilbPoints[Missing,0,1];",
            "unknown residue hyperelliptic curve",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
