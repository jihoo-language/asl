mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const FIXTURE: &str = include_str!("../examples/hyperelliptic-k3-semilinear-hilbert.asl");
const EXAMPLE_LINES: usize = 9;

#[test]
fn example_reports_the_enriques_square_and_two_other_actions() {
    assert_eq!(
        stdout(example("hyperelliptic-k3-semilinear-hilbert")),
        concat!(
            "rvK3FrobHilbFixed[Enriques, C, 1, 1] = 12\n",
            "rvK3FrobHilbFixed[Enriques, C, 2, 1] = 232\n",
            "rvK3FrobHilbTrace[Enriques, C, 2, 1, 0] = 1\n",
            "rvK3FrobHilbTrace[Enriques, C, 2, 1, 2] = 5\n",
            "rvK3FrobHilbTrace[Enriques, C, 2, 1, 4] = 100\n",
            "rvK3FrobHilbTrace[Enriques, C, 2, 1, 6] = 45\n",
            "rvK3FrobHilbTrace[Enriques, C, 2, 1, 8] = 81\n",
            "rvK3FrobHilbFixed[Shift, C, 2, 2] = 10452\n",
            "rvK3FrobHilbFixed[Hyperbolic, C, 2, 1] = 5344\n",
        ),
    );
}

#[test]
fn identity_word_recovers_the_ordinary_hilbert_counts_and_all_grades() {
    let mut source = format!("{FIXTURE} rvk3autom Idle[B]{{}};");
    for (base, r) in [("C", 1), ("C", 2), ("C", 3), ("L", 1)] {
        for n in 0..=4 {
            writeln!(
                source,
                "rvKummerHilbPoints[{base},{n},{r}]; rvK3FrobHilbFixed[Idle,{base},{n},{r}];"
            )
            .unwrap();
            for i in 0..=4 * n + 1 {
                writeln!(
                    source,
                    "rvKummerHilbTrace[{base},{n},{r},{i}]; rvK3FrobHilbTrace[Idle,{base},{n},{r},{i}];"
                )
                .unwrap();
            }
        }
    }
    let output = stdout(interpret(&source));
    let values: Vec<_> = output.lines().skip(EXAMPLE_LINES).map(value).collect();
    assert_eq!(values.len(), 440);
    for pair in values.chunks_exact(2) {
        assert_eq!(pair[0], pair[1]);
    }
}

#[test]
fn enriques_hilbert_square_agrees_with_actual_surface_cycles() {
    let output = stdout(interpret(&format!(
        r"{FIXTURE}
        let psi=\x:RVKummerK3(L)=>rvK3AutoPoint[Enriques](rvKummerK3Frob[C,1](x));
        assert (\x:RVKummerK3(L)=>psi(psi(x)))==(\x:RVKummerK3(L)=>x);
        show psi;
        rvK3FrobHilbFixed[Enriques,C,1,2];
        "
    )));
    let mut lines = output.lines();
    let length_one = value(lines.next().unwrap());
    let length_two = value(lines.next().unwrap());
    let map = lines.nth(EXAMPLE_LINES - 2).unwrap();
    let points: Vec<_> = map
        .match_indices("RVKummerK3(L){")
        .map(|(start, _)| {
            let end = start + map[start..].find('}').unwrap() + 1;
            &map[start..end]
        })
        .collect();
    let quadratic_count = value(lines.next().unwrap());
    assert_eq!(points.len() as i128, 2 * quadratic_count);
    let fixed = points
        .chunks_exact(2)
        .filter(|pair| pair[0] == pair[1])
        .count() as i128;
    assert_eq!(fixed, length_one);

    // Choose two fixed supports or one length-two orbit. A double support
    // has 3+1 fixed tangent directions under its semilinear P¹ action.
    let different_supports = fixed * (fixed - 1) / 2 + (quadratic_count - fixed) / 2;
    assert_eq!(length_two, different_supports + 4 * fixed);
    assert!(lines.next().is_none());
}

#[test]
fn squares_use_powers_of_the_composite_and_realize_the_symmetric_square_trace() {
    let mut source = FIXTURE.to_owned();
    source.push_str("rvK3Entropy[Hyperbolic]; rvKummerHilbBetti[2,2]; rvKummerHilbBetti[2,4];");
    for word in ["Enriques", "Shift", "Hyperbolic"] {
        for r in 1..=2 {
            writeln!(source, "rvK3FrobFixed[{word},C,{r}];").unwrap();
            writeln!(source, "rvK3FrobFixed[{word},C,{}];", 2 * r).unwrap();
            writeln!(source, "rvK3FrobHilbFixed[{word},C,2,{r}];").unwrap();
            for i in (0..=8).step_by(2) {
                writeln!(source, "rvK3FrobHilbTrace[{word},C,2,{r},{i}];").unwrap();
            }
            writeln!(source, "rvK3FrobHilbTrace[{word},C,2,{},2];", 2 * r).unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines().skip(EXAMPLE_LINES);
    assert!(lines.next().unwrap().contains("greater than 1"));
    let h2_rank = value(lines.next().unwrap());
    let h4_rank = value(lines.next().unwrap());
    assert_eq!(h4_rank, h2_rank * (h2_rank + 1) / 2);
    for word in ["Enriques", "Shift", "Hyperbolic"] {
        for r in 1..=2 {
            let q = 3i128.pow(r);
            let c1 = value(lines.next().unwrap());
            let c2 = value(lines.next().unwrap());
            let count = value(lines.next().unwrap());
            let traces: Vec<_> = lines.by_ref().take(5).map(value).collect();
            let next_h2 = value(lines.next().unwrap());
            assert_eq!(count, (c1 * c1 + c2) / 2 + q * c1);
            assert_eq!(count, traces.iter().sum::<i128>());
            assert_eq!(traces[2], (traces[1] * traces[1] + next_h2) / 2);
            assert_eq!(traces[3], q.pow(2) * traces[1]);
            assert_eq!(traces[4], q.pow(4));
            if word == "Shift" && r == 1 {
                // W^2 is the identity but (W F)^2 is a nonzero translation
                // followed by F^2. The ghost uses its 104 fixed points.
                assert_eq!((c1, c2, count), (8, 104, 108));
            }
        }
    }
    assert!(lines.next().is_none());
}

#[test]
fn the_hilbert_cube_uses_the_third_surface_power_and_length_three_mode() {
    let output = stdout(interpret(&format!(
        "{FIXTURE}
         rvK3FrobFixed[Shift,C,1];
         rvK3FrobFixed[Shift,C,2];
         rvK3FrobFixed[Shift,C,3];
         rvK3FrobHilbFixed[Shift,C,3,1];
         rvK3FrobHilbTrace[Shift,C,1,1,2];"
    )));
    let values: Vec<_> = output.lines().skip(EXAMPLE_LINES).map(value).collect();
    let [c1, c2, c3, count, h2] = values.as_slice() else {
        panic!("expected three surface counts, one Hilbert point count, and one trace");
    };
    let q = 3i128;
    assert_eq!(
        *count,
        (c1.pow(3) + 3 * c1 * c2 + 2 * c3) / 6 + q * c1.pow(2) + q.pow(2) * c1,
    );
    assert_eq!(*h2, c1 - 1 - q.pow(2));
    assert_eq!(*h2, -2);
}

#[test]
fn elementary_degrees_and_wide_arithmetic_are_checked_independently_of_counts() {
    let max = u128::MAX;
    let output = stdout(interpret(&format!(
        "{FIXTURE}
         rvK3FrobHilbFixed[Enriques,C,0,{max}];
         rvK3FrobHilbTrace[Enriques,C,0,{max},0];
         rvK3FrobHilbTrace[Enriques,C,2,{max},1];
         rvK3FrobHilbTrace[Enriques,C,2,{max},10];
         rvK3FrobHilbTrace[Enriques,C,{max},1,2];
         rvK3FrobHilbTrace[Enriques,C,{max},1,4];
         rvK3FrobHilbTrace[Enriques,C,4,1,4];
         rvk3autom Idle[B]{{}};
         rvK3FrobHilbTrace[Idle,C,1,41,2];
         rvKummerHilbTrace[C,1,41,2];"
    )));
    let results: Vec<_> = output.lines().skip(EXAMPLE_LINES).map(value).collect();
    assert_eq!(results[..5], [1, 1, 0, 0, 5]);
    assert_eq!(results[5], results[6]);
    assert_eq!(results[7], results[8]);
    assert!(results[7] != 0);

    for (statement, expected) in [
        (
            "rvK3FrobHilbFixed[Enriques,C,0,0];",
            "iteration count must be positive",
        ),
        (
            "rvK3FrobHilbFixed[Enriques,C,10000000000000,1];",
            "coefficient storage limit",
        ),
        (
            "rvK3FrobHilbFixed[Enriques,C,2,30];",
            "unsigned 128-bit range",
        ),
        (
            "rvK3FrobHilbTrace[Enriques,C,2,30,8];",
            "signed 128-bit range",
        ),
        (
            "rvK3FrobHilbTrace[Enriques,C,2,1000,8];",
            "exact intermediate exceeds 512 bits",
        ),
        (
            "rvhyperelliptic Other[F]{o,neg(o),z,z,z,o}; rvK3FrobHilbTrace[Enriques,Other,0,1,0];",
            "not the declared base extension",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{FIXTURE}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn value(line: &str) -> i128 {
    line.split_once(" = ").unwrap().1.parse().unwrap()
}
