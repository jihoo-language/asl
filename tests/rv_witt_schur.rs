mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const HUGE: u128 = u128::MAX;

#[test]
fn example_exhibits_nontrivial_schur_multiplicities_and_the_degree_six_lie_spectrum() {
    assert_eq!(
        stdout(example("witt-schur-branching")),
        concat!(
            "rvWittSchurMajorIndex[[3, 2, 1]] = [2, 3, 3, 2, 3, 3]\n",
            "rvWittCyclicSchurMultiplicity[[3, 2, 1], 6] = single=3, rational=6\n",
            "rvWittLieSchurMultiplicity[[3, 2, 1]] = 3\n",
            "rvWittSpectrumSchurGhost[[4, 2], 1] = 3\n",
            "rvWittSpectrumGhost[1] = 9\n",
        )
    );
}

fn partitions(weight: usize) -> Vec<Vec<usize>> {
    fn fill(weight: usize, largest: usize, rows: &mut Vec<usize>, all: &mut Vec<Vec<usize>>) {
        if weight == 0 {
            all.push(rows.clone());
            return;
        }
        for row in (1..=weight.min(largest)).rev() {
            rows.push(row);
            fill(weight - row, row, rows, all);
            rows.pop();
        }
    }
    let mut out = Vec::new();
    fill(weight, weight, &mut Vec::new(), &mut out);
    out
}

fn standard_major_counts(rows: &[usize]) -> Vec<usize> {
    fn fill(
        rows: &[usize],
        lengths: &mut [usize],
        previous: usize,
        i: usize,
        major: usize,
        out: &mut [usize],
    ) {
        if i > out.len() {
            out[major % out.len()] += 1;
            return;
        }
        for (row, &width) in rows.iter().enumerate() {
            if lengths[row] == width || (row > 0 && lengths[row] == lengths[row - 1]) {
                continue;
            }
            lengths[row] += 1;
            let descent = if row > previous { i - 1 } else { 0 };
            fill(rows, lengths, row, i + 1, major + descent, out);
            lengths[row] -= 1;
        }
    }
    if rows.is_empty() {
        return vec![1];
    }
    let mut counts = vec![0; rows.iter().sum()];
    fill(
        rows,
        &mut vec![0; rows.len()],
        rows.len(),
        1,
        0,
        &mut counts,
    );
    counts
}

fn phi(n: usize) -> usize {
    (1..=n)
        .filter(|&d| {
            let (mut a, mut b) = (n, d);
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a == 1
        })
        .count()
}

#[test]
fn q_hook_residues_and_all_cyclotomic_multiplicities_match_standard_tableaux_through_seven() {
    let mut source = String::from("rvWittSchurMajorIndex[[]];\n");
    let mut expected = String::from("rvWittSchurMajorIndex[[]] = [1]\n");
    for n in 1..=7usize {
        for shape in partitions(n) {
            let counts = standard_major_counts(&shape);
            writeln!(source, "rvWittSchurMajorIndex[{shape:?}];").unwrap();
            writeln!(expected, "rvWittSchurMajorIndex[{shape:?}] = {counts:?}").unwrap();
            writeln!(source, "rvWittLieSchurMultiplicity[{shape:?}];").unwrap();
            writeln!(
                expected,
                "rvWittLieSchurMultiplicity[{shape:?}] = {}",
                counts[1 % n]
            )
            .unwrap();
            for d in (1..=n).filter(|&d| n.is_multiple_of(d)) {
                let single = counts[(n / d) % n];
                writeln!(source, "rvWittCyclicSchurMultiplicity[{shape:?},{d}];").unwrap();
                writeln!(
                    expected,
                    "rvWittCyclicSchurMultiplicity[{shape:?}, {d}] = single={single}, rational={}",
                    single * phi(d)
                )
                .unwrap();
            }
        }
    }
    assert_eq!(stdout(interpret(&source)), expected);
}

fn supertableau_eigenvalues(rows: &[usize]) -> [Vec<i128>; 2] {
    struct Tableaux {
        cells: Vec<(usize, usize)>,
        values: Vec<Vec<usize>>,
        spectra: [Vec<i128>; 2],
    }
    impl Tableaux {
        fn fill(&mut self, index: usize, eigenvalue: i128, parity: usize) {
            const EIGENVALUES: [i128; 3] = [2, 3, 5];
            let Some(&(i, j)) = self.cells.get(index) else {
                self.spectra[parity].push(eigenvalue);
                return;
            };
            for (letter, &root) in EIGENVALUES.iter().enumerate() {
                let odd = letter == 2;
                if j > 0 {
                    let left = self.values[i][j - 1];
                    if letter < left || letter == left && odd {
                        continue;
                    }
                }
                if i > 0 {
                    let up = self.values[i - 1][j];
                    if letter < up || letter == up && !odd {
                        continue;
                    }
                }
                self.values[i][j] = letter;
                self.fill(index + 1, eigenvalue * root, parity ^ usize::from(odd));
            }
        }
    }
    let mut tableaux = Tableaux {
        cells: rows
            .iter()
            .enumerate()
            .flat_map(|(i, &width)| (0..width).map(move |j| (i, j)))
            .collect(),
        values: rows.iter().map(|&width| vec![0; width]).collect(),
        spectra: [Vec::new(), Vec::new()],
    };
    tableaux.fill(0, 1, 0);
    tableaux.spectra
}

fn polynomial(roots: &[i128]) -> Vec<i128> {
    let mut coefficients = vec![1];
    for &root in roots {
        coefficients.push(0);
        for j in (1..coefficients.len()).rev() {
            coefficients[j] -= root * coefficients[j - 1];
        }
    }
    coefficients
}

fn transpose(rows: &[usize]) -> Vec<usize> {
    (1..=rows[0])
        .map(|j| rows.iter().filter(|&&width| width >= j).count())
        .collect()
}

#[test]
fn signed_schur_spectra_match_supertableaux_and_transposition() {
    let mut source = String::from("let x = rvWittSpectrum{[1,-5,6],[1,-5]};\n");
    let mut expected = Vec::new();
    for n in 1..=4 {
        for (i, rows) in partitions(n).into_iter().enumerate() {
            let [even, odd] = supertableau_eigenvalues(&rows);
            let (p, q) = (polynomial(&even), polynomial(&odd));
            let label = format!("t{n}x{i}");
            writeln!(source,"let {label} = rvWittSpectrumSchur[{rows:?}](x); assert {label} == rvWittSpectrum{{{p:?},{q:?}}};").unwrap();
            let transposed = transpose(&rows);
            let right = format!("rvWittSpectrumSchur[{transposed:?}](x)");
            let right = if n.is_multiple_of(2) {
                right
            } else {
                format!("neg({right})")
            };
            writeln!(
                source,
                "assert rvWittSpectrumSchur[{rows:?}](neg(x)) == {right};"
            )
            .unwrap();
            for k in 1..=3 {
                writeln!(
                    source,
                    "rvWittSpectrumSchurGhost[{rows:?},{k}](x); rvWittSpectrumGhost[{k}]({label});"
                )
                .unwrap();
                let sum = |roots: &[i128]| roots.iter().map(|a| a.pow(k)).sum::<i128>();
                let ghost = sum(&even) - sum(&odd);
                expected.extend([ghost, ghost]);
            }
        }
    }
    let output = stdout(interpret(&source));
    assert_eq!(
        output
            .lines()
            .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
            .collect::<Vec<i128>>(),
        expected
    );
}

fn sum_terms(terms: impl Iterator<Item = String>) -> String {
    terms
        .reduce(|a, b| format!("add({a},{b})"))
        .unwrap_or_else(|| "zero[RVWittSpectrum]".into())
}

#[test]
fn schur_branching_reconstructs_every_cyclic_factor_and_lie_character() {
    for signed in [false, true] {
        let input = if signed {
            "rvWittSpectrum{[1,-1],[1,1]}"
        } else {
            "rvWittSpectrum{[1,-2,1]}"
        };
        let mut source = format!("let x = {input};\n");
        for n in 1..=6usize {
            let shapes = partitions(n);
            let lift = |r, copies| {
                sum_terms(shapes.iter().flat_map(|shape| {
                    let count = standard_major_counts(shape)[r] * copies;
                    std::iter::repeat_n(format!("rvWittSpectrumSchur[{shape:?}](x)"), count)
                }))
            };
            let lie = lift(1 % n, 1);
            writeln!(source, "assert rvWittSpectrumLie[{n}](x) == {lie};").unwrap();
            for d in (1..=n).filter(|&d| n.is_multiple_of(d)) {
                let factor = lift((n / d) % n, phi(d));
                writeln!(
                    source,
                    "assert rvWittSpectrumCyclic[{n},{d}](x) == {factor};"
                )
                .unwrap();
            }
        }
        assert!(stdout(interpret(&source)).is_empty());
    }
}

#[test]
fn elementary_shapes_and_ghosts_avoid_unneeded_rank_and_exponent_work() {
    let shorter = HUGE - 1;
    let output = stdout(interpret(&format!(
        r"
        let x = rvWittSpectrumLine[-1];
        let split = rvWittSpectrum{{[1,0,-1]}};
        assert rvWittSpectrumSchur[[]](zero[RVWittSpectrum]) == rvWittSpectrumLine[1];
        assert rvWittSpectrumSchur[[{HUGE}]](x) == x;
        assert rvWittSpectrumSchur[[{shorter},1]](x) == zero[RVWittSpectrum];
        assert rvWittSpectrumSchur[[1,1]](neg(x)) == rvWittSpectrumLine[1];
        rvWittSpectrumSchurGhost[[],{HUGE}](x);
        rvWittSpectrumSchurGhost[[2],{HUGE}](split);
        rvWittSpectrumSchurGhost[[25],1](rvWittSpectrum{{[1,-3,3,-1]}});
        rvWittLieSchurMultiplicity[[{HUGE}]];
        rvWittCyclicSchurMultiplicity[[{HUGE}],{HUGE}];
    "
    )));
    assert_eq!(
        output
            .lines()
            .map(|l| l.split_once(" = ").unwrap().1)
            .collect::<Vec<_>>(),
        ["1", "1", "351", "0", "single=0, rational=0"]
    );
}

#[test]
fn invalid_shapes_orders_and_degree_bounds_are_diagnosed() {
    for (input, expected) in [
        (
            "show rvWittSpectrumSchur[[1,2]](rvWittSpectrumLine[1]);",
            "weakly decreasing",
        ),
        ("rvWittSchurMajorIndex[[2,0]];", "positive"),
        ("rvWittSchurMajorIndex[[257]];", "weight limit"),
        ("rvWittLieSchurMultiplicity[[]];", "nonempty shape"),
        (
            "rvWittCyclicSchurMultiplicity[[4,1],2];",
            "positive divisor",
        ),
        (
            "rvWittCyclicSchurMultiplicity[[],1];",
            "nonempty Schur shape",
        ),
        (
            "rvWittSpectrumSchurGhost[[2],0](zero[RVWittSpectrum]);",
            "ghost index must be positive",
        ),
        (
            "show rvWittSpectrumSchur[[1]](S.hi);",
            "rvWittSpectrumSchur",
        ),
        (
            "show rvWittSpectrumSchur[[25]](rvWittSpectrum{[1,-3,3,-1]});",
            "determinant degree limit",
        ),
        (
            "rvWittSchurMajorIndex[[340282366920938463463374607431768211455,1]];",
            "unsigned 128-bit range",
        ),
    ] {
        let error = diagnostic(interpret(input));
        assert!(error.contains(expected), "{error}");
    }
}
