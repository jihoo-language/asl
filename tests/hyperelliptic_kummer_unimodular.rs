mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zeroW=rvUnramElem[W]{};
    let oneW=rvUnramElem[W]{1};
    let a=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[L](zeroW,zeroW);
    let f1=rvJacPoint[L](oneW,zeroW);
    let e2=add(add(e1,f1),rvJacPoint[L](neg(oneW),zeroW));
    let f2=add(add(e1,f1),rvJacPoint[L](a,zeroW));
    rvjactorsionbasis B[L,2](neg(oneW)){e1,e2,f1,f2};
";

const PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

#[test]
fn actual_gram_matrix_is_even_integral_and_unimodular() {
    let mut source = format!("{MARKED} rvK3LatticeGram[B];");
    for i in 0..22 {
        writeln!(source, "let b{i}=rvK3LatticeBasis[B,{i}];").unwrap();
    }
    for i in 0..22 {
        for j in 0..22 {
            writeln!(source, "rvK3LatticePair(b{i},b{j});").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    let report = lines.next().unwrap();
    let gram = report
        .strip_prefix(
            "rvK3LatticeGram[B] = rank=22, signature=(3, 19), index=2048, determinant=-1, gram=",
        )
        .unwrap();
    let entries: Vec<i128> = gram
        .split(|c: char| c != '-' && !c.is_ascii_digit())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(entries.len(), 22 * 22);
    let mut matrix = [[0; 22]; 22];
    for i in 0..22 {
        for j in 0..22 {
            matrix[i][j] = entries[22 * i + j];
            let pairing = lines
                .next()
                .unwrap()
                .strip_prefix("rvK3LatticePair = ")
                .unwrap()
                .parse::<i128>()
                .unwrap();
            assert_eq!(matrix[i][j], pairing);
            assert_eq!(pairing, entries[22 * j + i]);
        }
        assert_eq!(matrix[i][i] % 2, 0);
    }
    assert!(lines.next().is_none());
    assert_eq!(determinant(matrix), -1);
}

#[test]
fn reflection_inverse_transvection_and_cross_pairings() {
    let output = stdout(example("hyperelliptic-kummer-unimodular"));
    assert!(
        output.starts_with(concat!(
            "rvK3LatticePair = -1\n",
            "rvK3LatticePair = -1\n",
            "rvK3LatticePair = -1\n",
            "rvK3LatticeGram[B] = rank=22, signature=(3, 19), index=2048, determinant=-1, gram=[[",
        )),
        "{output}"
    );
    let theta = literal([0, 2, 0, 0, 2, 0], [0; 16]);
    let output = stdout(interpret(&format!(
        "{MARKED} let theta={theta}; assert rvK3Transvection[15,{}](theta)==theta; rvK3LatticePair(theta,theta);",
        i128::MIN,
    )));
    assert_eq!(output, "rvK3LatticePair = 4\n");
}

#[test]
fn all_mask_actions_match_their_integral_exterior_square_and_node_map() {
    let w = [1, -2, 3, -4, 5, 6];
    let b = std::array::from_fn(|v| {
        let value: i128 = PAIRS
            .iter()
            .zip(w)
            .map(|(&(i, j), value)| (((v >> i) & (v >> j)) & 1) as i128 * value)
            .sum();
        (value & 1) + 2 * v as i128 - 8
    });
    let mut source = format!("{MARKED} let x={};", literal(w, b));
    for mask in 0..16 {
        for amount in [-2, -1, 1, 3] {
            let u: [i128; 4] = std::array::from_fn(|i| (mask >> i) & 1);
            let g = [u[2], u[3], -u[0], -u[1]];
            let inverse_transpose: [[i128; 4]; 4] = std::array::from_fn(|i| {
                std::array::from_fn(|j| i128::from(i == j) - amount * g[i] * u[j])
            });
            let shifted_w = PAIRS.map(|(i, j)| {
                PAIRS
                    .iter()
                    .zip(w)
                    .map(|(&(p, q), v)| {
                        v * (inverse_transpose[i][p] * inverse_transpose[j][q]
                            - inverse_transpose[i][q] * inverse_transpose[j][p])
                    })
                    .sum()
            });
            let mut shifted_b = [0; 16];
            for (node, value) in b.into_iter().enumerate() {
                let coordinates: [i128; 4] = std::array::from_fn(|i| ((node >> i) & 1) as i128);
                let contraction: i128 = coordinates.iter().zip(g).map(|(&x, y)| x * y).sum();
                let image = (0..4).fold(0, |v, i| {
                    v | ((((coordinates[i] + amount * u[i] * contraction) & 1) as usize) << i)
                });
                shifted_b[image] = value;
            }
            writeln!(
                source,
                "assert rvK3Transvection[{mask},{amount}](x)=={};",
                literal(shifted_w, shifted_b)
            )
            .unwrap();
            writeln!(
                source,
                "assert rvK3Transvection[{mask},{}](rvK3Transvection[{mask},{amount}](x))==x;",
                -amount
            )
            .unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn wide_pairing_and_reflection_keep_exact_cancellations() {
    let a = i128::MAX - 1;
    let left = literal([a, a, 0, 0, 0, 0], [0; 16]);
    let right = literal([0, 0, 0, 0, a, a], [0; 16]);
    let mut b = [0; 16];
    b[15] = 2;
    let root = literal([0; 6], b);
    b[15] = a;
    let x = literal([0; 6], b);
    assert_eq!(stdout(interpret(&format!(
        "{MARKED} let root={root}; let x={x}; assert rvK3Reflect(root,x)==neg(x); rvK3LatticePair({left},{right});"
    ))), "rvK3LatticePair = 0\n");
    assert!(
        diagnostic(interpret(&format!("{MARKED} rvK3LatticePair({x},{x});")))
            .contains("pairing exceeds the signed i128 range")
    );
}

#[test]
fn rejects_bad_gluing_bad_indices_nonroots_and_nominal_mismatches() {
    let surplus = literal([i128::MAX - 1, 0, 0, 0, 0, 0], [0; 16]);
    let increment = literal([2, 0, 0, 0, 0, 0], [0; 16]);
    for (source, expected) in [
        ("show rvK3Lattice[B]{[0],[0]};".into(), "requires six wedge coordinates"),
        ("show rvK3Lattice[B]{[0,0,0,0,0,0],[0]};".into(), "requires sixteen exceptional coordinates"),
        ("show rvK3Lattice[B]{[1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};".into(), "affine gluing condition"),
        ("show rvK3LatticeBasis[B,22];".into(), "between 0 and 21"),
        ("show rvK3Transvection[16,0](zero[RVK3Lattice(B)]);".into(), "between 0 and 15"),
        ("show rvK3Reflect(zero[RVK3Lattice(B)],rvK3LatticeBasis[B,0]);".into(), "requires a root of square"),
        ("rvjactorsionbasis D[L,2](neg(oneW)){e1,e2,f1,f2}; rvK3LatticePair(rvK3LatticeBasis[B,0],rvK3LatticeBasis[D,0]);".into(), "expected type RVK3Lattice(B)"),
        ("show rvK3Lattice[Missing]{[0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};".into(), "unknown Jacobian torsion basis"),
        (format!("show add({surplus},{increment});"), "coordinate exceeds the signed i128 range"),
    ] {
        let problem = diagnostic(interpret(&format!("{MARKED}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}

fn literal(w: [i128; 6], b: [i128; 16]) -> String {
    format!("rvK3Lattice[B]{{{w:?},{b:?}}}")
}

fn determinant(mut matrix: [[i128; 22]; 22]) -> i128 {
    let mut previous = 1;
    let mut sign = 1;
    for column in 0..21 {
        let Some(pivot) = (column..22).find(|&row| matrix[row][column] != 0) else {
            return 0;
        };
        if pivot != column {
            matrix.swap(pivot, column);
            sign = -sign;
        }
        let value = matrix[column][column];
        for r in column + 1..22 {
            for c in column + 1..22 {
                let updated = value * matrix[r][c] - matrix[r][column] * matrix[column][c];
                assert_eq!(updated % previous, 0);
                matrix[r][c] = updated / previous;
            }
        }
        previous = value;
    }
    sign * matrix[21][21]
}
