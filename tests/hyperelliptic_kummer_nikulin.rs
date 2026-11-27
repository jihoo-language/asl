mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const SHEAR: &str = r"
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

#[test]
fn shear_and_split_lattices_have_the_expected_integral_determinants() {
    let output = stdout(example("hyperelliptic-kummer-nikulin"));
    assert_eq!(
        output,
        concat!(
            "rvJacTorsionFrob[B, C, 1] = multiplier=1, matrix=[[1, 0, 0, 0], [0, 1, 0, 1], [0, 0, 1, 0], [0, 0, 0, 1]], characteristic=[1, 0, 0, 0, 1]\n",
            "rvKummerNodeLattice[B, C, 0] = orbits=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], rank=16, affine invariants=5, index=32, determinant=64, geometric rank=16, geometric determinant=64\n",
            "rvKummerNodeLattice[B, C, 1] = orbits=[1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2], rank=12, affine invariants=4, index=16, determinant=256, geometric rank=16, geometric determinant=64\n",
            "RVKummerDisc(B){5}\n",
            "PrimeScalar(2){1}\n",
        ),
    );
    let matrix_text = output
        .split_once("matrix=")
        .unwrap()
        .1
        .split_once(", characteristic")
        .unwrap()
        .0;
    let matrix: Vec<_> = matrix_text
        .bytes()
        .filter(|b| b.is_ascii_digit())
        .map(|b| b - b'0')
        .collect();
    assert_eq!(matrix.len(), 16);
    assert_eq!(fixed_gram_determinant(&matrix), 256);
}

#[test]
fn all_classes_have_the_expected_quadratic_and_polar_values() {
    let mut source = format!(
        r"{SHEAR}
        let q=\x:RVKummerDisc(B)=>rvKummerDiscQuadratic(x);
        let pair=\x:RVKummerDisc(B)=>\y:RVKummerDisc(B)=>rvKummerDiscPair(x,y);
        let polar=\x:RVKummerDisc(B)=>\y:RVKummerDisc(B)=>add(add(q(x),q(y)),q(add(x,y)));
        assert pair==polar;
        assert (\x:RVKummerDisc(B)=>neg(x))==(\x:RVKummerDisc(B)=>x);
        assert (\x:RVKummerDisc(B)=>add(x,x))==(\x:RVKummerDisc(B)=>zero[RVKummerDisc(B)]);
        assert (\x:RVKummerDisc(B)=>q(rvKummerDiscFrob[C,1](x)))==q;
        assert (\x:RVKummerDisc(B)=>rvKummerDiscFrob[L,7](x))==(\x:RVKummerDisc(B)=>x);
        show q;
        "
    );
    for bit in 0..6 {
        writeln!(
            source,
            r"show (\x:RVKummerDisc(B)=>rvKummerDiscPair(x,rvKummerDisc[B,{}]));",
            1 << bit
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    let q = lines.next().unwrap();
    let pair_columns: Vec<_> = lines.collect();
    assert_eq!(pair_columns.len(), 6);
    assert_eq!(q.matches("=> PrimeScalar(2){0}").count(), 36);
    assert_eq!(q.matches("=> PrimeScalar(2){1}").count(), 28);
    for mask in 0u8..64 {
        let weight = (0..16).filter(|&node| evaluates_to_one(mask, node)).count();
        assert!(q.contains(&format!(
            "RVKummerDisc(B){{{mask}}} => PrimeScalar(2){{{}}}",
            (weight / 2) % 2
        )));
        let mut nonzero_pairing = false;
        for (bit, output) in pair_columns.iter().enumerate() {
            let intersection = (0..16)
                .filter(|&node| evaluates_to_one(mask, node) && evaluates_to_one(1 << bit, node))
                .count();
            let value = intersection % 2;
            assert!(output.contains(&format!(
                "RVKummerDisc(B){{{mask}}} => PrimeScalar(2){{{value}}}"
            )));
            nonzero_pairing |= value == 1;
        }
        assert_eq!(nonzero_pairing, mask != 0);
    }
}

#[test]
fn order_three_frobenius_uses_the_inverse_permutation_on_functions() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            rvunram W[3,1]{1,2,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=neg(o);
            rvhyperelliptic C[F]{z,t,t,t,t,o};
            rvjacbase L[C,W];
            let zeroW=rvUnramElem[W]{};
            let oneW=rvUnramElem[W]{1};
            let a=rvUnramElem[W]{0,1};
            let e1=rvJacPoint[L](zeroW,zeroW);
            let f1=rvJacPoint[L](oneW,zeroW);
            let e2=add(add(e1,f1),rvJacPoint[L](a,zeroW));
            let f2=add(add(e1,f1),rvJacPoint[L](rvUnramFrob[1](a),zeroW));
            rvjactorsionbasis B[L,2](neg(oneW)){e1,e2,f1,f2};
            assert rvJacFrob[C,1](e2)==f2;
            assert rvJacFrob[C,1](f2)==add(e2,f2);
            let q=rvKummerDisc[B,1];
            assert rvKummerDiscFrob[C,1](q)==rvKummerDisc[B,5];
            assert rvKummerDiscFrob[C,2](q)==rvKummerDisc[B,4];
            let forward=\x:RVKummerDisc(B)=>rvKummerDiscFrob[C,1](x);
            let inverse=\x:RVKummerDisc(B)=>rvKummerDiscFrob[C,2](x);
            assert (\x:RVKummerDisc(B)=>inverse(forward(x)))==(\x:RVKummerDisc(B)=>x);
            assert (\x:RVKummerDisc(B)=>rvKummerDiscQuadratic(forward(x)))==(\x:RVKummerDisc(B)=>rvKummerDiscQuadratic(x));
            rvKummerNodeLattice[B,C,1];
            rvKummerNodeLattice[B,C,3];
            ",
        )),
        concat!(
            "rvKummerNodeLattice[B, C, 1] = orbits=[1, 1, 1, 1, 3, 3, 3, 3], rank=8, affine invariants=3, index=8, determinant=324, geometric rank=16, geometric determinant=64\n",
            "rvKummerNodeLattice[B, C, 3] = orbits=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], rank=16, affine invariants=5, index=32, determinant=64, geometric rank=16, geometric determinant=64\n",
        ),
    );
}

#[test]
fn validates_mask_nominal_basis_level_and_base_equation() {
    for (source, expected) in [
        ("show rvKummerDisc[B,64];", "between 0 and 63"),
        ("show rvKummerDisc[Missing,0];", "unknown Jacobian torsion basis"),
        ("rvjactorsionbasis Other[L,2](neg(oneW)){e1,e2,f1,f2}; show add(rvKummerDisc[B,0],rvKummerDisc[Other,0]);", "expected type RVKummerDisc(B)"),
        ("rvjactorsionbasis Other[L,2](neg(oneW)){e1,e2,f1,f2}; show rvKummerDiscPair(rvKummerDisc[B,0],rvKummerDisc[Other,0]);", "expected type RVKummerDisc(B)"),
        ("rvhyperelliptic D[F]{z,o,z,z,z,o}; show rvKummerDiscFrob[D,1](rvKummerDisc[B,1]);", "not the declared base extension"),
        ("rvhyperelliptic D[F]{z,o,z,z,z,o}; rvKummerNodeLattice[B,D,1];", "not the declared base extension"),
        ("rvhyperelliptic E[F]{z,neg(o),z,o}; rvjactorsionbasis Small[E,2](neg(o)){rvJacPoint[E](z,z),rvJacPoint[E](o,z)}; show zero[RVKummerDisc(Small)];", "certified genus-two level-2 basis"),
    ] {
        let problem = diagnostic(interpret(&format!("{SHEAR}{source}")));
        assert!(problem.contains(expected), "{problem}");
    }
}

fn evaluates_to_one(mask: u8, node: u8) -> bool {
    let mut bit = 0;
    let mut parity = 0;
    for i in 0..4 {
        for j in i + 1..4 {
            parity ^= ((mask >> bit) & (node >> i) & (node >> j)) & 1;
            bit += 1;
        }
    }
    parity == 1
}

fn fixed_gram_determinant(matrix: &[u8]) -> i128 {
    let image = |node: usize| {
        matrix
            .chunks_exact(4)
            .enumerate()
            .fold(0, |value, (i, row)| {
                let coordinate = row
                    .iter()
                    .enumerate()
                    .map(|(j, &entry)| ((node >> j) & 1) * usize::from(entry))
                    .sum::<usize>()
                    % 2;
                value | (coordinate << i)
            })
    };
    let mut lengths = Vec::new();
    let mut representatives = Vec::new();
    let mut seen = [false; 16];
    for node in 0..16 {
        if seen[node] {
            continue;
        }
        representatives.push(node);
        let mut length = 0;
        let mut next = node;
        while !seen[next] {
            seen[next] = true;
            length += 1;
            next = image(next);
        }
        lengths.push(length);
    }
    let rank = lengths.len();
    let mut echelon = vec![0usize; rank];
    for affine in 0usize..32 {
        let evaluate =
            |node: usize| ((affine & 1) as u32 + ((affine >> 1) & node).count_ones()) % 2;
        if (0..16).any(|node| evaluate(node) != evaluate(image(node))) {
            continue;
        }
        let mut word = representatives
            .iter()
            .enumerate()
            .fold(0, |word, (i, &r)| word | ((evaluate(r) as usize) << i));
        for (pivot, row) in echelon.iter_mut().enumerate() {
            if word & (1 << pivot) == 0 {
                continue;
            }
            if *row == 0 {
                *row = word;
                break;
            }
            word ^= *row;
        }
    }
    // Use the half-sum for each pivot and the original orbit sum elsewhere.
    let component = |i, k| {
        if echelon[i] == 0 {
            2 * i128::from(i == k)
        } else {
            ((echelon[i] >> k) & 1) as i128
        }
    };
    let mut gram = vec![vec![0i128; rank]; rank];
    for (i, row) in gram.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            let doubled = -(0..rank)
                .map(|k| lengths[k] * component(i, k) * component(j, k))
                .sum::<i128>();
            assert_eq!(doubled % 2, 0);
            *entry = doubled / 2;
        }
    }
    let mut previous = 1;
    for pivot in 0..rank - 1 {
        let value = gram[pivot][pivot];
        for i in pivot + 1..rank {
            for j in pivot + 1..rank {
                let numerator = value * gram[i][j] - gram[i][pivot] * gram[pivot][j];
                assert_eq!(numerator % previous, 0);
                gram[i][j] = numerator / previous;
            }
        }
        previous = value;
    }
    gram[rank - 1][rank - 1]
}
