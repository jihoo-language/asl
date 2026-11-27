mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const HUGE: u128 = u128::MAX;

fn numbers(output: &str) -> Vec<i128> {
    output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn example_reports_the_signed_square_and_normalizes_the_primitive_cyclic_part() {
    assert_eq!(
        stdout(example("witt-cyclic-lie")),
        concat!(
            "rvWittSpectrumGhost[1] = 2\n",
            "rvWittSpectrumGhost[1] = -1\n",
            "rvWittSpectrumCyclicGhost[6, 6, 1] = 18\n",
            "rvWittSpectrumLieGhost[6, 1] = 9\n",
            "rvWittSpectrumLieGhost[14, 1] = 1161\n",
        )
    );
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn cyclic_orbit_eigenvalues(length: u32, order: usize) -> [Vec<i128>; 2] {
    // The even eigenvalue 2 and odd eigenvalue 1 give a simple tensor basis.
    // Rotate actual words with the Koszul sign. On a signed orbit of length
    // L, the roots of the rotation satisfy z^L=sign and F is scalar.
    let all = 1u32 << length;
    let mut visited = vec![false; all as usize];
    let mut values = [Vec::new(), Vec::new()];
    let primitive_roots = (1..=order).filter(|&j| gcd(j, order) == 1).count();
    for word in 0..all {
        if visited[word as usize] {
            continue;
        }
        let odd_count = word.count_ones();
        let mut cycle_length: usize = 0;
        let mut negative = false;
        let mut image = word;
        loop {
            visited[image as usize] = true;
            negative ^= image & 1 != 0 && odd_count.is_multiple_of(2);
            image = (image >> 1) | ((image & 1) << (length - 1));
            cycle_length += 1;
            if image == word {
                break;
            }
        }
        let included = if negative {
            (2 * cycle_length).is_multiple_of(order) && !cycle_length.is_multiple_of(order)
        } else {
            cycle_length.is_multiple_of(order)
        };
        if included {
            let parity = (odd_count % 2) as usize;
            values[parity].extend(std::iter::repeat_n(
                1i128 << (length - odd_count),
                primitive_roots,
            ));
        }
    }
    values
}

fn characteristic(roots: &[i128]) -> Vec<i128> {
    let mut coeffs = vec![1];
    for &root in roots {
        coeffs.push(0);
        for j in (1..coeffs.len()).rev() {
            coeffs[j] -= root * coeffs[j - 1];
        }
    }
    coeffs
}

#[test]
fn signed_rotation_orbits_certify_all_projectors_and_their_direct_ghosts() {
    let through = 7u32;
    for n in [4u32, 6] {
        let mut source =
            String::from("let x = add(rvWittSpectrumLine[2],neg(rvWittSpectrumLine[1]));");
        let mut totals = vec![0; through as usize];
        let mut expected = Vec::new();
        for d in (1..=n).filter(|d| n.is_multiple_of(*d)) {
            let [even, odd] = cyclic_orbit_eigenvalues(n, d as usize);
            let p = characteristic(&even);
            let q = characteristic(&odd);
            let name = format!("c{d}");
            writeln!(source, "let {name} = rvWittSpectrumCyclic[{n},{d}](x); assert {name} == rvWittSpectrum{{{p:?},{q:?}}};").unwrap();
            for k in 1..=through {
                writeln!(
                    source,
                    "rvWittSpectrumGhost[{k}]({name}); rvWittSpectrumCyclicGhost[{n},{d},{k}](x);"
                )
                .unwrap();
                let contribution = |roots: &[i128]| roots.iter().map(|a| a.pow(k)).sum::<i128>();
                let ghost = contribution(&even) - contribution(&odd);
                totals[(k - 1) as usize] += ghost;
                expected.extend([ghost, ghost]);
            }
        }
        for k in 1..=through {
            assert_eq!(totals[(k - 1) as usize], (2i128.pow(k) - 1).pow(n));
        }
        writeln!(
            source,
            "let lie = rvWittSpectrumLie[{n}](x); assert c{n} == add(lie,lie);"
        )
        .unwrap();
        let primitive = &expected[expected.len() - 2 * through as usize..];
        let lie_ghosts: Vec<_> = primitive.iter().map(|ghost| ghost / 2).collect();
        for k in 1..=through {
            writeln!(
                source,
                "rvWittSpectrumLieGhost[{n},{k}](x); rvWittSpectrumGhost[{k}](lie);"
            )
            .unwrap();
        }
        expected.extend(lie_ghosts);
        assert_eq!(numbers(&stdout(interpret(&source))), expected);
    }
}

#[test]
fn all_cyclotomic_factors_add_to_the_tensor_power() {
    let output = stdout(interpret(
        r"
        let x = rvWittSpectrum{[1,-2,1]};
        let x2 = rvWittSpectrumTensor(x,x);
        let x3 = rvWittSpectrumTensor(x2,x);
        let y = add(rvWittSpectrumCyclic[6,1](x),rvWittSpectrumCyclic[6,2](x));
        let z = add(rvWittSpectrumCyclic[6,3](x),rvWittSpectrumCyclic[6,6](x));
        assert add(y,z) == rvWittSpectrumTensor(x3,x3);
    ",
    ));
    assert!(output.is_empty());
}

#[test]
fn quadratic_projectors_recover_symmetric_and_exterior_powers() {
    let output = stdout(interpret(
        r"
        let x = rvWittSpectrum{[1,-10,31,-30]};
        let ext = rvWittSpectrumExterior[2](x);
        let sym = add(rvWittSpectrumTensor(x,x),neg(ext));
        assert rvWittSpectrumCyclic[1,1](x) == x;
        assert rvWittSpectrumLie[1](x) == x;
        assert rvWittSpectrumCyclic[2,1](x) == sym;
        assert rvWittSpectrumCyclic[2,2](x) == ext;
        assert rvWittSpectrumLie[2](x) == ext;
        let odd = neg(rvWittSpectrumLine[2]);
        assert rvWittSpectrumCyclic[2,1](odd) == zero[RVWittSpectrum];
        assert rvWittSpectrumCyclic[2,2](odd) == rvWittSpectrumLine[4];
        assert rvWittSpectrumLie[2](odd) == rvWittSpectrumLine[4];
        assert rvWittSpectrumLie[3](odd) == zero[RVWittSpectrum];
        assert rvWittSpectrumCyclic[3,1](odd) == neg(rvWittSpectrumLine[8]);
        assert rvWittSpectrumCyclic[4,2](odd) == rvWittSpectrumLine[16];
        assert rvWittSpectrumCyclic[4,4](odd) == zero[RVWittSpectrum];
    ",
    ));
    assert!(output.is_empty());
}

#[test]
fn parity_ranks_satisfy_the_super_pbw_hilbert_series() {
    let mut source =
        String::from("let x = add(rvWittSpectrumLine[1],neg(rvWittSpectrumLine[-1]));");
    let max = 8;
    for n in 1..=max {
        writeln!(source, "rvWittSpectrumData(rvWittSpectrumLie[{n}](x));").unwrap();
    }
    let output = stdout(interpret(&source));
    let mut hilbert = vec![0; max + 1];
    hilbert[0] = 1;
    for (i, line) in output.lines().enumerate() {
        // Frobenius equals parity here, so the two determinants have no common roots.
        let fields: Vec<_> = line.split(", ").collect();
        let even: usize = fields[0].rsplit_once('=').unwrap().1.parse().unwrap();
        let odd: usize = fields[1].rsplit_once('=').unwrap().1.parse().unwrap();
        let n = i + 1;
        for _ in 0..even {
            for k in n..=max {
                hilbert[k] += hilbert[k - n];
            }
        }
        for _ in 0..odd {
            for k in (n..=max).rev() {
                hilbert[k] += hilbert[k - n];
            }
        }
    }
    assert_eq!(hilbert, (0..=max).map(|k| 1 << k).collect::<Vec<_>>());
}

#[test]
fn composite_mobius_ghosts_stay_available_above_the_determinant_bound() {
    let output = stdout(interpret(&format!(
        r"
        let binary = rvWittSpectrum{{[1,-2,1]}};
        rvWittSpectrumLieGhost[12,3](binary);
        rvWittSpectrumLieGhost[14,1](binary);
        rvWittSpectrumCyclicGhost[14,14,1](binary);
        let split = rvWittSpectrum{{[1,0,-1]}};
        rvWittSpectrumCyclicGhost[2,1,{HUGE}](split);
        rvWittSpectrumLieGhost[2,{HUGE}](split);
        let odd = neg(rvWittSpectrumLine[-1]);
        assert rvWittSpectrumLie[{HUGE}](rvWittSpectrumLine[100]) == zero[RVWittSpectrum];
        assert rvWittSpectrumCyclic[{HUGE},{HUGE}](rvWittSpectrumLine[100]) == zero[RVWittSpectrum];
        assert rvWittSpectrumCyclic[{HUGE},1](odd) == odd;
        rvWittSpectrumCyclicGhost[{HUGE},1,{HUGE}](odd);
    "
    )));
    assert_eq!(numbers(&output), [335, 1161, 6966, 1, -1, 1]);
}

#[test]
fn invalid_orders_nonpositive_indices_and_arithmetic_limits_are_reported() {
    for (input, expected) in [
        ("show rvWittSpectrumCyclic[4,3](rvWittSpectrumLine[1]);", "positive divisor"),
        ("show rvWittSpectrumCyclic[4,0](rvWittSpectrumLine[1]);", "positive divisor"),
        ("show rvWittSpectrumLie[0](rvWittSpectrumLine[1]);", "length must be positive"),
        ("rvWittSpectrumCyclicGhost[3,3,0](zero[RVWittSpectrum]);", "ghost index must be positive"),
        ("show rvWittSpectrumLie[2](S.hi);", "rvWittSpectrumLie"),
        ("rvWittSpectrumCyclicGhost[2,2,1](S.hi);", "RVWittSpectrum"),
        ("show rvWittSpectrumLie[14](rvWittSpectrum{[1,-2,1]});", "determinant degree limit"),
        ("show rvWittSpectrumCyclic[13,13](rvWittSpectrum{[1,-2,1]});", "determinant degree limit"),
        ("rvWittSpectrumLieGhost[800,1](rvWittSpectrum{[1,-2,1]});", "512 bits"),
        ("rvWittSpectrumLieGhost[170141183460469231731687303715884105727,1](rvWittSpectrum{[1,1,1]});", "trial division limit"),
    ] {
        let result = diagnostic(interpret(input));
        assert!(result.contains(expected), "{result}");
    }
}
