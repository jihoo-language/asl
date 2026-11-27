mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn number(line: &str) -> u128 {
    line.rsplit_once(" = ").unwrap().1.parse().unwrap()
}

fn binomial(n: u128, k: u128) -> u128 {
    if k > n {
        return 0;
    }
    (1..=k).fold(1, |out, i| out * (n + 1 - i) / i)
}

#[test]
fn example_counts_polarized_invariants_and_the_extra_orientation_line() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-hodge-invariants")),
        concat!(
            "rvHilbHodgeInvariant[2, 2, 3] = 3\n",
            "rvHilbHodgeInvariant[2, 4, 3] = 7\n",
            "rvHilbPrimitiveHodgeInvariant[2, 4, 3] = 4\n",
            "rvHilbHodgeInvariant[3, 4, 3] = 10\n",
            "rvHilbHodgeInvariantMultiplicity[21][1, 1] = 211\n",
        )
    );
}

#[test]
fn hilbert_square_and_cube_follow_classical_tensor_invariants_in_all_ranks() {
    let mut program = String::new();
    let mut expected = Vec::new();
    for rho in 1..=21 {
        // The traceless symmetric square has binomial(rho+1,2) invariants
        // after the metric line has been removed. In the cube's middle
        // degree add Sym^3(W), exterior^2(W), and one trivial line.
        let square = binomial(rho + 1, 2);
        let orientation = u128::from(rho == 21);
        let cube_middle = binomial(rho + 2, 3) + rho + binomial(rho, 2) + orientation + 1;
        for (n, d, all, primitive) in [
            (2, 0, 1, 1),
            (2, 2, rho, rho - 1),
            (2, 4, square + 1, square + 1 - rho),
            (2, 6, rho, 0),
            (2, 8, 1, 0),
            (3, 2, rho, rho - 1),
            (3, 4, square + rho + 1, square + 1),
            (3, 6, cube_middle, cube_middle - square - rho - 1),
            (3, 8, square + rho + 1, 0),
        ] {
            writeln!(program, "rvHilbHodgeInvariant[{n},{d},{rho}];").unwrap();
            writeln!(program, "rvHilbPrimitiveHodgeInvariant[{n},{d},{rho}];").unwrap();
            expected.extend([all, primitive]);
        }
    }
    let actual: Vec<_> = stdout(interpret(&program)).lines().map(number).collect();
    assert_eq!(actual, expected);
}

#[test]
fn exterior_powers_detect_the_connected_orientation_invariant() {
    let mut program = String::new();
    let mut expected = Vec::new();
    for rho in 1..=21 {
        let t = 23 - rho;
        for k in 1u128..=11 {
            let shape = vec!["1"; k as usize].join(",");
            writeln!(
                program,
                "rvHilbHodgeInvariantMultiplicity[{rho}]{{{shape}}};"
            )
            .unwrap();
            // In exterior^k(A + T), only exterior^0(T) and exterior^t(T)
            // are SO(T) invariant. Full O(T) would omit exterior^t(T).
            let orientation = k.checked_sub(t).map_or(0, |j| binomial(rho, j));
            expected.push(binomial(rho, k) + orientation);
        }
    }
    let actual: Vec<_> = stdout(interpret(&program)).lines().map(number).collect();
    assert_eq!(actual, expected);
}

#[test]
fn rank_two_invariants_equal_the_hodge_diagonal_and_its_primitive_difference() {
    for n in 2..=9 {
        let mut program = format!("rvHilbLLVHodgeTable[{n}];\n");
        for d in 0..=2 * n {
            writeln!(program, "rvHilbHodgeInvariant[{n},{},21];", 2 * d).unwrap();
            writeln!(program, "rvHilbPrimitiveHodgeInvariant[{n},{},21];", 2 * d).unwrap();
        }
        let reported = stdout(interpret(&program));
        let mut lines = reported.lines();
        let table = lines
            .next()
            .unwrap()
            .split_once("hodge=[[")
            .unwrap()
            .1
            .strip_suffix("]]")
            .unwrap();
        let mut older = 0;
        for (p, row) in table.split("], [").enumerate() {
            let diagonal = row.split(", ").nth(p).unwrap().parse::<u128>().unwrap();
            assert_eq!(number(lines.next().unwrap()), diagonal, "n={n}, p={p}");
            let primitive = if p <= n { diagonal - older } else { 0 };
            assert_eq!(number(lines.next().unwrap()), primitive, "n={n}, p={p}");
            older = diagonal;
        }
        assert_eq!(lines.next(), None);
    }
}

// For energy through nine, the first relation among Gram generators for
// T of dimension at least three costs more than nine. The first relation
// involving an orientation determinant for dim(T)=3 also costs more than
// nine, and no product of two such determinants or determinant for dim(T)>3
// fits. This gives an independent truncated invariant-theory product.
fn truncated_invariant_product(cap: usize, rho: usize) -> Vec<Vec<u128>> {
    let width = 2 * cap + 1;
    let mut count = vec![vec![0; width]; cap + 1];
    count[0][cap] = 1;
    let mut add_generator = |size: usize, charge: i32| {
        for n in size..=cap {
            for b in 0..width {
                let previous = b as i32 - charge;
                if (0..width as i32).contains(&previous) {
                    count[n][b] += count[n - size][previous as usize];
                }
            }
        }
    };
    for m in 1..=cap {
        add_generator(m, 1);
        add_generator(m, -1);
        for _ in 1..rho {
            add_generator(m, 0);
        }
        for l in m..=cap - m {
            add_generator(m + l, 0);
        }
    }
    if rho == 20 {
        let free = count.clone();
        for a in 1..=cap {
            for b in a + 1..=cap {
                for c in b + 1..=cap {
                    let size = a + b + c;
                    for n in size..=cap {
                        for (out, &base) in count[n].iter_mut().zip(&free[n - size]) {
                            *out += base;
                        }
                    }
                }
            }
        }
    }
    count
}

#[test]
fn generic_hodge_counts_match_the_independent_gram_and_orientation_product() {
    let cap = 9;
    for rho in [1, 2, 10, 19, 20] {
        let expected = truncated_invariant_product(cap, rho);
        let mut program = String::new();
        let mut dims = Vec::new();
        for (n, graded) in expected.iter().enumerate().skip(2) {
            for d in 0..=2 * n {
                let dimension = graded[cap - n + d];
                let primitive = if d == 0 {
                    dimension
                } else if d <= n {
                    dimension - graded[cap - n + d - 1]
                } else {
                    0
                };
                writeln!(program, "rvHilbHodgeInvariant[{n},{},{rho}];", 2 * d).unwrap();
                writeln!(
                    program,
                    "rvHilbPrimitiveHodgeInvariant[{n},{},{rho}];",
                    2 * d
                )
                .unwrap();
                dims.extend([dimension, primitive]);
            }
        }
        let reported = stdout(interpret(&program));
        assert_eq!(
            reported.lines().map(number).collect::<Vec<_>>(),
            dims,
            "rho={rho}"
        );
    }
}

#[test]
fn validates_geometry_domains_and_checks_polynomial_work() {
    let output = stdout(interpret(
        r"
        rvHilbHodgeInvariant[0,0,1];
        rvHilbPrimitiveHodgeInvariant[0,0,21];
        rvHilbHodgeInvariant[0,2,21];
        rvHilbHodgeInvariant[3,13,21];
        rvHilbPrimitiveHodgeInvariant[3,-2,1];
        rvHilbHodgeInvariant[1000000000000000000,7,21];
        rvHilbHodgeInvariantMultiplicity[1]{340282366920938463463374607431768211455};
        rvHilbHodgeInvariantMultiplicity[1]{340282366920938463463374607431768211455,1,0};
        ",
    ));
    assert_eq!(
        output.lines().map(number).collect::<Vec<_>>(),
        [1, 1, 0, 0, 0, 0, 1, 0]
    );

    for (program, message) in [
        (
            "rvHilbHodgeInvariant[2,2,0];",
            "polarization rank must be in 1..=21",
        ),
        (
            "rvHilbPrimitiveHodgeInvariant[2,2,22];",
            "polarization rank must be in 1..=21",
        ),
        ("rvHilbHodgeInvariant[1,2,1];", "whose LLV type is D12"),
        (
            "rvHilbHodgeInvariantMultiplicity[3]{1,2};",
            "type B11 highest weight",
        ),
        (
            "rvHilbHodgeInvariant[1000000000000000000,2,21];",
            "polynomial computation limit",
        ),
        (
            "rvHilbHodgeInvariantMultiplicity[2]{340282366920938463463374607431768211455};",
            "polynomial computation limit",
        ),
    ] {
        assert!(
            diagnostic(interpret(program)).contains(message),
            "{program}"
        );
    }
}
