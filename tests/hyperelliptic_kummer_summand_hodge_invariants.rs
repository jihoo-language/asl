mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_cube_separates_divisors_quadratic_forms_and_unitary_skew_forms() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-summand-hodge-invariants")),
        concat!(
            "rvHilbLLVHodgeInvariant[3, 6, 15, real, 2][1, 1] = irrep=106 contribution=106\n",
            "rvHilbLLVHodgeInvariant[3, 6, 15, cm, 2][1, 1] = irrep=108 contribution=108\n",
            "rvHilbLLVHodgeInvariantDecomposition[3, 6, 15, cm, 2] = components=[([3], 1, 710, 710), ([1, 1], 1, 108, 108)] total=818 divisors=680 quadraticContractions=30 nonVerbitsky=108\n",
            "rvHilbLLVHodgeInvariant[6, 12, 3][2] = irrep=7 contribution=14\n",
        )
    );
    let result = stdout(interpret(
        r"
        rvHilbLLVHodgeInvariant[3,4,21]{1,1};
        rvHilbLLVHodgeInvariant[3,6,21]{1,1};
        rvHilbLLVHodgeInvariant[3,8,21]{1,1};
        rvHilbLLVHodgeInvariant[3,2,21]{1,1};
        rvHilbLLVHodgeInvariant[2,4,21]{1,1};
        ",
    ));
    assert_eq!(
        result,
        concat!(
            "rvHilbLLVHodgeInvariant[3, 4, 21][1, 1] = irrep=21 contribution=21\n",
            "rvHilbLLVHodgeInvariant[3, 6, 21][1, 1] = irrep=212 contribution=212\n",
            "rvHilbLLVHodgeInvariant[3, 8, 21][1, 1] = irrep=21 contribution=21\n",
            "rvHilbLLVHodgeInvariant[3, 2, 21][1, 1] = irrep=0 contribution=0\n",
            "rvHilbLLVHodgeInvariant[2, 4, 21][1, 1] = irrep=212 contribution=0\n",
        )
    );
}

#[test]
fn harmonic_young_restriction_subtracts_the_trace_part() {
    // S_(2,1) V = V times wedge^2 V minus wedge^3 V. Orthogonal
    // restriction then subtracts the copy of V from the contraction.
    assert_eq!(
        stdout(interpret("rvHilbLLVHodgeInvariantDecomposition[4,8,3];")),
        "rvHilbLLVHodgeInvariantDecomposition[4, 8, 3] = components=[([4], 1, 22, 22), ([2, 1], 1, 14, 14), ([2], 1, 7, 7), ([], 1, 1, 1)] total=44 divisors=15 quadraticContractions=7 nonVerbitsky=22\n"
    );
}

fn fields(line: &str) -> [u128; 4] {
    let fields = line.split_once(" total=").unwrap().1;
    fields
        .split(' ')
        .map(|field| field.rsplit('=').next().unwrap().parse().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

#[test]
fn exact_summands_and_splits_agree_with_the_whole_cohomology_counts() {
    for (rho, selector, command, index) in [
        (1, "", "rvHilbHodgeInvariant", ""),
        (21, "", "rvHilbHodgeInvariant", ""),
        (14, ",real,3", "rvHilbRealMultiplicationHodge", ",3"),
        (2, ",real,7", "rvHilbRealMultiplicationHodge", ",7"),
        (15, ",cm,2", "rvHilbComplexMultiplicationHodge", ",2"),
        (3, ",cm,5", "rvHilbComplexMultiplicationHodge", ",5"),
    ] {
        let mut program = String::new();
        for n in 2..=9 {
            for p in 0..=n {
                // The other half is equivalent by Poincare duality.
                let d = 2 * p;
                writeln!(
                    program,
                    "rvHilbLLVHodgeInvariantDecomposition[{n},{d},{rho}{selector}];"
                )
                .unwrap();
                writeln!(program, "{command}[{n},{d},{rho}{index}];").unwrap();
            }
        }
        let result = stdout(interpret(&program));
        let lines = result.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 2 * (2..=9).map(|n| n + 1).sum::<usize>());
        for pair in lines.chunks_exact(2) {
            let [total, divisors, quadratic, other] = fields(pair[0]);
            assert_eq!(total, divisors + quadratic + other, "{}", pair[0]);
            assert_eq!(
                total,
                pair[1]
                    .rsplit_once(" = ")
                    .unwrap()
                    .1
                    .parse::<u128>()
                    .unwrap(),
                "{}",
                pair[0]
            );
        }
    }
}

#[test]
fn sparse_boundary_queries_and_report_validation() {
    let report = stdout(interpret(
        r"
        rvHilbLLVHodgeInvariantDecomposition[0,0,3];
        rvHilbLLVHodgeInvariantDecomposition[3,9,3];
        rvHilbLLVHodgeInvariantDecomposition[3,14,3];
        rvHilbLLVHodgeInvariant[2,-2,3]{4};
        rvHilbLLVHodgeInvariant[1000000000000000000,2,15,cm,2]{1000000000000000000};
        rvHilbLLVHodgeInvariantDecomposition[1000000000000000000,0,3];
        rvHilbLLVHodgeInvariant[2,0,3]{2,1,0};
        ",
    ));
    let lines = report.lines().collect::<Vec<_>>();
    assert_eq!(fields(lines[0]), [1, 1, 0, 0]);
    assert_eq!(fields(lines[1]), [0, 0, 0, 0]);
    assert_eq!(fields(lines[2]), [0, 0, 0, 0]);
    assert!(lines[3].ends_with("irrep=3 contribution=0"));
    assert!(lines[4].ends_with("irrep=15 contribution=15"));
    assert_eq!(fields(lines[5]), [1, 1, 0, 0]);
    assert!(lines[6].ends_with("irrep=3 contribution=0"));

    for (code, expected) in [
        (
            "rvHilbLLVHodgeInvariantDecomposition[1,2,3];",
            "whose LLV type is D12",
        ),
        (
            "rvHilbLLVHodgeInvariantDecomposition[2,4,0];",
            "rank must be in 1..=21",
        ),
        (
            "rvHilbLLVHodgeInvariantDecomposition[2,4,20,real,3];",
            "at least three transcendental directions",
        ),
        (
            "rvHilbLLVHodgeInvariantDecomposition[2,4,2,cm,1];",
            "CM field degree 2e",
        ),
        (
            "rvHilbLLVHodgeInvariantDecomposition[2,4,2,other,1];",
            "real or cm",
        ),
        (
            "rvHilbLLVHodgeInvariant[2,4,3]{1,2};",
            "type B12 highest weight",
        ),
        (
            "rvHilbLLVHodgeInvariant[2,4,3]{1000000000000000000,1};",
            "polynomial computation limit",
        ),
    ] {
        assert!(diagnostic(interpret(code)).contains(expected), "{code}");
    }
}

#[test]
fn individual_torus_invariants_match_the_independent_bigraded_character() {
    let mut program = String::new();
    for (n, weight) in [
        (3, "1,1"),
        (4, "2,1"),
        (6, "2,2"),
        (7, "3,1"),
        (8, "3,2,1"),
        (9, "2,2,2"),
    ] {
        for p in 0..=2 * n {
            writeln!(
                program,
                "rvHilbLLVHodgeInvariant[{n},{},21]{{{weight}}};",
                2 * p
            )
            .unwrap();
            writeln!(
                program,
                "rvHilbLLVHodgeInvariant[{n},{},21,cm,1]{{{weight}}};",
                2 * p
            )
            .unwrap();
            writeln!(program, "rvHilbLLVHodge[{n},{p},{p}]{{{weight}}};").unwrap();
        }
    }
    let result = stdout(interpret(&program));
    let lines = result.lines().collect::<Vec<_>>();
    for triple in lines.chunks_exact(3) {
        let values =
            [triple[0], triple[1], triple[2]].map(|line| line.split_once(" = ").unwrap().1);
        assert_eq!(values[0], values[2], "{}", triple[0]);
        assert_eq!(values[1], values[2], "{}", triple[1]);
    }
}
