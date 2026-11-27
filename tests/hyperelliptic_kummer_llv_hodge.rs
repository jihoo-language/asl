mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn matrix(text: &str) -> Vec<Vec<u128>> {
    text.strip_prefix("[[")
        .unwrap()
        .strip_suffix("]]")
        .unwrap()
        .split("], [")
        .map(|row| {
            row.split(", ")
                .map(|value| value.parse().unwrap())
                .collect()
        })
        .collect()
}

fn goettsche(n: usize) -> Vec<Vec<u128>> {
    let diameter = 2 * n;
    let mut table = vec![vec![vec![0; diameter + 1]; diameter + 1]; n + 1];
    table[0][0][0] = 1;
    for m in 1..=n {
        for (p, q, multiplicity) in [(0, 0, 1), (2, 0, 1), (0, 2, 1), (1, 1, 20), (2, 2, 1)] {
            let a = p + m - 1;
            let b = q + m - 1;
            for _ in 0..multiplicity {
                for weight in m..=n {
                    for p in a..=diameter {
                        for q in b..=diameter {
                            table[weight][p][q] += table[weight - m][p - a][q - b];
                        }
                    }
                }
            }
        }
    }
    table.pop().unwrap()
}

fn expected_n3() -> [Vec<Vec<u128>>; 2] {
    let mut weights = vec![(0i32, 0i32); 21];
    weights.extend([(1, 0), (-1, 0), (0, 1), (0, -1)]);
    let mut exterior = vec![vec![0u128; 7]; 7];
    let mut symmetric = exterior.clone();
    let coords = |a: i32, b: i32| ((3 + a + b) as usize, (3 + a - b) as usize);
    for (i, &(a, b)) in weights.iter().enumerate() {
        for (j, &(c, d)) in weights.iter().enumerate().skip(i) {
            if i != j {
                let (p, q) = coords(a + c, b + d);
                exterior[p][q] += 1;
            }
            for &(e, f) in weights.iter().skip(j) {
                let (p, q) = coords(a + c + e, b + d + f);
                symmetric[p][q] += 1;
            }
        }
        let (p, q) = coords(a, b);
        symmetric[p][q] -= 1;
    }
    [exterior, symmetric]
}

#[test]
fn example_resolves_the_middle_bidegree() {
    let result = stdout(example("hyperelliptic-kummer-llv-hodge"));
    let lines = result.lines().collect::<Vec<_>>();
    assert_eq!(
        lines[0],
        "rvHilbLLVHodge[3, 3, 3][1, 1] = irrep=212 contribution=212"
    );
    assert_eq!(
        lines[1],
        "rvHilbLLVHodge[3, 3, 3][3] = irrep=1792 contribution=1792"
    );
    assert_eq!(lines[2],"rvHilbLLVHodgeDecomposition[3, 3, 3] = components=[([3], 1, 1792, 1792), ([1, 1], 1, 212, 212)] total=2004");
    assert_eq!(
        matrix(lines[3].split_once("hodge=").unwrap().1),
        goettsche(2)
    );
}

#[test]
fn every_bidegree_matches_the_independent_goettsche_hodge_product() {
    let mut source = String::new();
    for n in 2..=9 {
        writeln!(source, "rvHilbLLVHodgeTable[{n}];").unwrap();
    }
    let result = stdout(interpret(&source));
    assert_eq!(result.lines().count(), 8);
    for (n, line) in (2..=9).zip(result.lines()) {
        let got = matrix(line.split_once("hodge=").unwrap().1);
        let expected = goettsche(n);
        assert_eq!(got, expected, "n={n}");
        for p in 0..=2 * n {
            for q in 0..=2 * n {
                assert_eq!(
                    got[p][q], got[q][p],
                    "Hodge symmetry at n={n}, p={p}, q={q}"
                );
                assert_eq!(
                    got[p][q],
                    got[2 * n - p][2 * n - q],
                    "Poincaré duality at n={n}, p={p}, q={q}"
                );
            }
        }
    }
}

#[test]
fn the_hilbert_cube_irreps_match_exterior_and_harmonic_cube_characters() {
    let result = stdout(interpret(
        "rvHilbLLVHodgeTable[3]{1,1}; rvHilbLLVHodgeTable[3]{3};",
    ));
    let expected = expected_n3();
    assert_eq!(result.lines().count(), 2);
    for (line, expected) in result.lines().zip(expected) {
        let (_, table) = line.split_once("irrep=").unwrap();
        let (unweighted, weighted) = table.split_once(" contribution=").unwrap();
        assert_eq!(matrix(unweighted), expected);
        assert_eq!(matrix(weighted), expected);
    }
}

#[test]
fn weighted_tables_and_bidegree_decompositions_track_multiplicity() {
    let result = stdout(interpret(
        r"
        rvHilbLLVHodge[6,6,6]{2};
        rvHilbLLVHodgeDecomposition[6,6,6];
        rvHilbLLVHodgeTable[6]{2};
        rvHilbLLVHodgeTable[3]{2};
        rvHilbLLVHodge[3,2,4]{1,1};
        rvHilbLLVHodge[3,2,2]{1,1};
        rvHilbLLVHodge[3,3,1]{1,1};
        rvHilbLLVHodge[3,1,1]{1,1};
    ",
    ));
    let lines = result.lines().collect::<Vec<_>>();
    assert_eq!(
        lines[0],
        "rvHilbLLVHodge[6, 6, 6][2] = irrep=232 contribution=464"
    );
    assert!(lines[1].contains("([2], 2, 232, 464)"));
    assert!(lines[1].ends_with("total=481504"));
    for (line, factor) in [(lines[2], 2u128), (lines[3], 0)] {
        let (_, table) = line.split_once("irrep=").unwrap();
        let (raw, contribution) = table.split_once(" contribution=").unwrap();
        let expected = matrix(raw)
            .into_iter()
            .map(|row| row.into_iter().map(|value| factor * value).collect())
            .collect::<Vec<Vec<_>>>();
        assert_eq!(matrix(contribution), expected);
    }
    assert_eq!(
        &lines[4..],
        &[
            "rvHilbLLVHodge[3, 2, 4][1, 1] = irrep=21 contribution=21",
            "rvHilbLLVHodge[3, 2, 2][1, 1] = irrep=21 contribution=21",
            "rvHilbLLVHodge[3, 3, 1][1, 1] = irrep=1 contribution=1",
            "rvHilbLLVHodge[3, 1, 1][1, 1] = irrep=0 contribution=0",
        ]
    );
}

#[test]
fn parity_absent_formal_weights_and_vacuum_have_exact_zero_rules() {
    assert_eq!(
        stdout(interpret(
            r"
        rvHilbLLVHodge[0,0,0]{};
        rvHilbLLVHodge[3,3,2]{3};
        rvHilbLLVHodge[3,-1,-1]{3};
        rvHilbLLVHodge[3,7,7]{4};
        rvHilbLLVHodge[1000,1001,1000]{2};
        rvHilbLLVHodgeDecomposition[1000,-1,0];
        rvHilbLLVHodgeTable[0];
    "
        )),
        concat!(
            "rvHilbLLVHodge[0, 0, 0][] = irrep=1 contribution=1\n",
            "rvHilbLLVHodge[3, 3, 2][3] = irrep=0 contribution=0\n",
            "rvHilbLLVHodge[3, -1, -1][3] = irrep=0 contribution=0\n",
            "rvHilbLLVHodge[3, 7, 7][4] = irrep=1 contribution=0\n",
            "rvHilbLLVHodge[1000, 1001, 1000][2] = irrep=0 contribution=0\n",
            "rvHilbLLVHodgeDecomposition[1000, -1, 0] = components=[] total=0\n",
            "rvHilbLLVHodgeTable[0] = hodge=[[1]]\n",
        )
    );
}

#[test]
fn first_length_invalid_partitions_and_budgets_are_diagnosed() {
    for (input, expected) in [
        ("rvHilbLLVHodge[1,1,1]{1};", "whose LLV type is D12"),
        (
            "rvHilbLLVHodgeDecomposition[1,0,0];",
            "whose LLV type is D12",
        ),
        ("rvHilbLLVHodgeTable[1]{};", "whose LLV type is D12"),
        (
            "rvHilbLLVHodge[2,1,1]{1,2};",
            "highest weight must be weakly decreasing",
        ),
        (
            "rvHilbLLVHodgeTable[3]{10000000000};",
            "exceeds the polynomial computation limit",
        ),
        (
            "rvHilbLLVHodgeTable[10000000000];",
            "exceeds the polynomial computation limit",
        ),
    ] {
        let err = diagnostic(interpret(input));
        assert!(err.contains(expected), "{err}");
    }
}
