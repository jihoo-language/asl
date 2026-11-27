mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EVEN: &str = include_str!("../examples/hyperelliptic-generalized-kummer-translation.asl");
const ODD: &str = include_str!("../examples/hyperelliptic-generalized-kummer-cohomology.asl");

#[test]
fn example_reports_a_character_and_affine_frobenius_counts() {
    assert_eq!(
        stdout(example("hyperelliptic-generalized-kummer-translation")),
        concat!(
            "rvGenKummerTranslateChar[3, 3, 4] = 27\n",
            "rvGenKummerTranslateFixed[3, 3] = 27 components of type K_0, Euler characteristic=27\n",
            "rvGenKummerTranslateFixed[4, 2] = 8 components of type K_1, Euler characteristic=192\n",
            "rvGenKummerTranslatePoints[B, C, 1] = 32, affine torsion fixed points=8\n",
            "rvGenKummerTranslatePoints[B, C, 1] = 8, affine torsion fixed points=0\n",
            "rvGenKummerTranslateFrob[B, C, 1, 2] = -2\n",
            "rvGenKummerTranslateFrob[B, C, 1, 2] = 22\n",
        )
    );
}

#[test]
fn geometric_characters_satisfy_duality_and_the_fixed_locus_euler_formula() {
    let mut source = String::new();
    for n in [1, 2, 3, 4, 6] {
        let last = 4 * (n - 1);
        for i in 0..=last {
            writeln!(source, "rvGenKummerBetti[{n},{i}];").unwrap();
        }
        for h in (1..=n).filter(|h| n % h == 0) {
            writeln!(source, "rvGenKummerTranslateFixed[{n},{h}];").unwrap();
            for i in 0..=last {
                writeln!(source, "rvGenKummerTranslateChar[{n},{h},{i}];").unwrap();
            }
        }
    }

    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    for n in [1i128, 2, 3, 4, 6] {
        let width = (4 * (n - 1) + 1) as usize;
        let betti: Vec<_> = lines.by_ref().take(width).map(number).collect();
        for h in (1..=n).filter(|h| n % h == 0) {
            let length = n / h;
            let euler = n.pow(3) * (1..=length).filter(|d| length % d == 0).sum::<i128>();
            let fixed_locus = lines.next().unwrap();
            assert!(fixed_locus.ends_with(&format!(
                "{} components of type K_{}, Euler characteristic={euler}",
                h.pow(3),
                length - 1,
            )));
            let trace: Vec<_> = lines.by_ref().take(width).map(number).collect();
            assert_eq!(alternating(&trace), euler);
            assert!(trace.iter().eq(trace.iter().rev()));
            for (&character, &dimension) in trace.iter().zip(&betti) {
                assert!(character.abs() <= dimension);
            }
            match (n, h) {
                (_, 1) => assert_eq!(trace, betti),
                (2, 2) => assert_eq!(trace, [1, 0, 6, 0, 1]),
                (3, 3) => assert_eq!(trace, [1, 0, 7, 8, 27, 8, 7, 0, 1]),
                (4, 2) => assert_eq!((trace[4], euler), (betti[4], 192)),
                (4, 4) => assert_eq!((trace[4], euler), (betti[4] - 16, 64)),
                _ => (),
            }
        }
    }
    assert!(lines.next().is_none());
}

#[test]
fn affine_level_two_counts_match_actual_resolved_kummer_maps() {
    let output = stdout(interpret(&format!(
        r"{EVEN}
        let frobenius = \x:RVKummerK3(D)=>rvKummerK3Frob[C,1](x);
        let soluble = \x:RVKummerK3(D)=>rvKummerK3Translate[B,2](frobenius(x));
        let insoluble = \x:RVKummerK3(D)=>rvKummerK3Translate[B,1](frobenius(x));
        show soluble;
        show insoluble;
        rvGenKummerTranslatePoints[B,C,1](f2);
        rvGenKummerTranslatePoints[B,C,2](e2);
        rvGenKummerTranslateFrob[B,C,2,2](e2);
        rvGenKummerTranslatePoints[B,C,2](zero[RVJac(D)]);
        rvGenKummerLevel[B,C,2];
        "
    )));
    let mut lines = output.lines().skip(7);
    for expected in [32, 8] {
        let map = lines.next().unwrap();
        let points: Vec<_> = map
            .match_indices("RVKummerK3(D){")
            .map(|(start, _)| {
                let end = start + map[start..].find('}').unwrap() + 1;
                &map[start..end]
            })
            .collect();
        // Both translations here are F3-rational and square to the identity,
        // so every twisted fixed point is already in the quadratic cover.
        assert_eq!(points.len(), 2 * 248);
        assert_eq!(
            points
                .chunks_exact(2)
                .filter(|pair| pair[0] == pair[1])
                .count(),
            expected,
        );
    }
    assert_eq!(
        lines.map(number).collect::<Vec<_>>(),
        [8, 104, 22, 248, 248]
    );
}

#[test]
fn prime_level_three_twist_changes_only_the_middle_trace() {
    let mut source = ODD.to_owned();
    for i in 0..=8 {
        writeln!(
            source,
            "rvGenKummerTrace[B,C,2,{i}]; rvGenKummerTranslateFrob[B,C,2,{i}](e1);"
        )
        .unwrap();
    }
    source.push_str("rvGenKummerTranslatePoints[B,C,2](e1);");
    let output = stdout(interpret(&source));
    let values: Vec<_> = output.lines().skip(9).map(number).collect();
    let ordinary: Vec<_> = values[..18].iter().step_by(2).copied().collect();
    let twisted: Vec<_> = values[1..18].iter().step_by(2).copied().collect();
    for i in 0..=8 {
        assert_eq!(
            ordinary[i] - twisted[i],
            if i == 4 { 81 * 25 * 25 } else { 0 }
        );
    }
    assert_eq!((twisted[3], twisted[4]), (-1000, 16_875));
    assert_eq!(alternating(&twisted), 543_051);
    assert_eq!(values[18], alternating(&twisted));
}

#[test]
fn rejects_invalid_orders_and_translations_and_respects_bounds() {
    let output = stdout(interpret(
        "rvGenKummerTranslateChar[340282366920938463463374607431768211455,1,2];
         rvGenKummerTranslateChar[2,2,340282366920938463463374607431768211455];",
    ));
    assert_eq!(output.lines().map(number).collect::<Vec<_>>(), [7, 0]);
    for (setup, statement, expected) in [
        (
            "",
            "rvGenKummerTranslateChar[0,1,2];",
            "length must be positive",
        ),
        ("", "rvGenKummerTranslateFixed[4,0];", "positive divisor"),
        ("", "rvGenKummerTranslateChar[4,3,100];", "positive divisor"),
        (
            "",
            "rvGenKummerTranslateChar[10000,10000,10000];",
            "coefficient storage limit",
        ),
        (
            "",
            "rvGenKummerTranslateFixed[1000000000000000000,1];",
            "unsigned 128-bit range",
        ),
        (
            EVEN,
            "rvGenKummerTranslateFrob[B,C,0,2](e2);",
            "extension degree must be positive",
        ),
        (
            EVEN,
            "rvGenKummerTranslatePoints[B,C,1](rvJacPoint[C](z,z));",
            "expected type RVJac(D)",
        ),
        (
            EVEN,
            "rvhyperelliptic X[F]{z,o,z,z,z,o}; rvGenKummerTranslateFrob[B,X,1,5](e2);",
            "not the declared base extension",
        ),
        (
            ODD,
            "rvGenKummerTranslateFrob[B,C,1,9](rvJacPoint[S](rvUnramElem[W]{},rvUnramElem[W]{}));",
            "requires a class killed by 3",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn number(line: &str) -> i128 {
    line.split_once(" = ")
        .unwrap()
        .1
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

fn alternating(values: &[i128]) -> i128 {
    values
        .iter()
        .enumerate()
        .map(|(i, &value)| if i.is_multiple_of(2) { value } else { -value })
        .sum()
}
