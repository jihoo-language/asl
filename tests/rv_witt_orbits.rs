mod common;

use std::collections::HashMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const GEOMETRY: &str = include_str!("../examples/hyperelliptic-witt-orbits.asl");
const HUGE: u128 = u128::MAX;

fn geometry() -> &'static str {
    let start = GEOMETRY.find("rvunram F").unwrap();
    let end = GEOMETRY.find("rvK3FrobHilbOrbits").unwrap();
    &GEOMETRY[start..end]
}

fn numbers(output: &str) -> Vec<i128> {
    output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn example_connects_signed_necklaces_and_primitive_geometric_orbits() {
    assert_eq!(
        stdout(example("hyperelliptic-witt-orbits")),
        concat!(
            "RVWittSpectrum{[1, 0, 0, -1], [1]}\n",
            "rvWittSpectrumGhost[6] = 3\n",
            "rvWittSpectrumNecklace[3] = 1\n",
            "rvWittSpectrumNecklace[6] = 0\n",
            "rvK3FrobHilbOrbits[Idle, C, 1, 2] = 108\n",
            "rvK3FrobHilbOrbits[Moving, C, 1, 2] = 48\n",
            "rvK3FrobHilbOrbits[Moving, C, 2, 2] = 5172\n",
            "rvWittSpectrumData = numerator degree=24, denominator degree=0, virtual rank=24\n",
            "rvWittSpectrumNecklace[2] = 48\n",
        )
    );
}

#[test]
fn verschiebung_satisfies_frobenius_and_tensor_identities() {
    let input = r"
        let x = rvWittSpectrum{[1,-2],[1,1]};
        let y = rvWittSpectrumLine[3];
        let vx = rvWittSpectrumVerschiebung[4](x);
        let vy = rvWittSpectrumVerschiebung[6](y);
        assert vx == rvWittSpectrum{[1,0,0,0,-2],[1,0,0,0,1]};
        assert rvWittSpectrumVerschiebung[1](x) == x;
        assert rvWittSpectrumVerschiebung[3](rvWittSpectrumVerschiebung[2](x)) == rvWittSpectrumVerschiebung[6](x);
        assert rvWittSpectrumVerschiebung[4](add(x,y)) == add(vx,rvWittSpectrumVerschiebung[4](y));
        let a = rvWittSpectrumVerschiebung[2](rvWittSpectrumAdams[3](x));
        assert rvWittSpectrumAdams[6](vx) == add(a,a);
        assert rvWittSpectrumAdams[3](vx) == rvWittSpectrumVerschiebung[4](rvWittSpectrumAdams[3](x));
        let b = rvWittSpectrumVerschiebung[12](rvWittSpectrumTensor(rvWittSpectrumAdams[3](x),rvWittSpectrumAdams[2](y)));
        assert rvWittSpectrumTensor(vx,vy) == add(b,b);
        rvWittSpectrumGhost[3](vx);
        rvWittSpectrumGhost[4](vx);
        rvWittSpectrumGhost[8](vx);
        rvWittSpectrumGhost[12](vx);
    ";
    assert_eq!(numbers(&stdout(interpret(input))), [0, 12, 12, 36]);
}

#[test]
fn necklaces_recover_euler_products_and_can_be_negative() {
    let mut source = String::from(
        r"
        let sign = rvWittSpectrum{[1,0,-1],[1,-1]};
        let x = add(rvWittSpectrumVerschiebung[3](rvWittSpectrumLine[1]),neg(rvWittSpectrumVerschiebung[2](rvWittSpectrumLine[1])));
        rvWittSpectrumNecklace[1](sign);
        rvWittSpectrumNecklace[2](sign);
        rvWittSpectrumNecklace[1](x);
        rvWittSpectrumNecklace[2](x);
        rvWittSpectrumNecklace[3](x);
        rvWittSpectrumNecklace[6](x);
        let binary = rvWittSpectrumLine[2];
    ",
    );
    for k in [1, 2, 3, 4, 5, 6, 8, 12] {
        writeln!(source, "rvWittSpectrumNecklace[{k}](binary);").unwrap();
    }
    writeln!(source, "rvWittSpectrumNecklace[{HUGE}](rvWittSpectrumVerschiebung[3](rvWittSpectrumLine[1])); rvWittSpectrumNecklace[{HUGE}](zero[RVWittSpectrum]); assert rvWittSpectrumVerschiebung[{HUGE}](zero[RVWittSpectrum]) == zero[RVWittSpectrum];").unwrap();
    assert_eq!(
        numbers(&stdout(interpret(&source))),
        [-1, 1, 0, -1, 1, 0, 2, 1, 2, 3, 6, 9, 30, 335, 0, 0]
    );
}

fn permutation<'a>(map: &'a str, curve: &str) -> HashMap<&'a str, &'a str> {
    let prefix = format!("RVKummerK3({curve}){{");
    let points: Vec<_> = map
        .match_indices(&prefix)
        .map(|(start, _)| {
            let end = start + map[start..].find('}').unwrap() + 1;
            &map[start..end]
        })
        .collect();
    assert_eq!(points.len() % 2, 0);
    points
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect()
}

fn orbit_count(perm: &HashMap<&str, &str>, period: usize) -> usize {
    let primitive = perm
        .keys()
        .filter(|&&point| {
            let mut image = point;
            for step in 1..=period {
                image = perm[image];
                if image == point {
                    return step == period;
                }
            }
            false
        })
        .count();
    assert_eq!(primitive % period, 0);
    primitive / period
}

#[test]
fn ordinary_primitive_orbits_match_geometric_closed_points() {
    let output = stdout(interpret(&format!(
        r"
        {}
        rvunram E[3,1]{{1,2,0,1}};
        rvjacbase U[C,E];
        rvK3FrobHilbOrbits[Idle,C,1,1];
        rvK3FrobHilbOrbits[Idle,C,1,2];
        rvK3FrobHilbOrbits[Idle,C,1,3];
        show (\x:RVKummerK3(L)=>rvKummerK3Frob[C,1](x));
        card (\x:RVKummerK3(C)=>S.hi);
        card (\x:RVKummerK3(U)=>S.hi);
        assert (\x:RVKummerK3(U)=>rvKummerK3Frob[C,1](rvKummerK3Frob[C,1](rvKummerK3Frob[C,1](x)))) == (\x:RVKummerK3(U)=>x);
        rvk3autom MovingProject[B]{{project(8),project(10)}};
        let psi=\x:RVKummerK3(L)=>rvK3AutoPoint[MovingProject](rvKummerK3Frob[C,1](x));
        assert (\x:RVKummerK3(L)=>psi(psi(x))) == (\x:RVKummerK3(L)=>x);
        show psi;
        rvK3FrobHilbOrbits[MovingProject,C,1,1];
        rvK3FrobHilbOrbits[MovingProject,C,1,2];
    ",
        geometry()
    )));
    let lines: Vec<_> = output.lines().collect();
    let quadratic = permutation(lines[3], "L");
    // Frobenius swaps these projectors. The composite is conjugate to
    // ordinary Frobenius, so its first two periods also lie over F9.
    let moving = permutation(lines[6], "L");
    assert_eq!(quadratic.len(), 248);
    for period in [1, 2] {
        assert_eq!(
            numbers(lines[period - 1]),
            [orbit_count(&quadratic, period) as i128]
        );
        assert_eq!(
            numbers(lines[6 + period]),
            [orbit_count(&moving, period) as i128]
        );
    }
    let points = |line: &str| {
        let (yes, total) = line
            .strip_prefix("card = ")
            .unwrap()
            .split_once(" / ")
            .unwrap();
        assert_eq!(yes, total);
        yes.parse::<i128>().unwrap()
    };
    assert_eq!(
        numbers(lines[2]),
        [(points(lines[5]) - points(lines[4])) / 3]
    );
}

#[test]
fn hilbert_power_counts_agree_with_surface_cycles_and_full_spectrum() {
    let mut input = format!(
        r"
        {}
        rvhilbstratum Surface{{(1,1)}};
        let actual = rvK3FrobHilbSpectrum[Moving,C,1,1];
        let h2 = rvK3StratumSpectrum[Moving,C,Surface,1,2];
        assert actual == add(rvWittSpectrumLine[1],add(h2,rvWittSpectrumLine[9]));
        assert rvK3FrobHilbSpectrum[Moving,C,0,{HUGE}] == rvWittSpectrumLine[1];
        rvK3FrobHilbOrbits[Moving,C,0,1];
        rvK3FrobHilbOrbits[Moving,C,0,{HUGE}];
    ",
        geometry()
    );
    for period in 1..=4 {
        writeln!(input, "rvK3FrobHilbOrbits[Moving,C,1,{period}]; rvWittSpectrumNecklace[{period}](actual); rvK3FrobHilbOrbits[Moving,C,2,{period}];").unwrap();
    }
    let output = numbers(&stdout(interpret(&input)));
    assert_eq!(output[..2], [1, 0]);
    // The four nontrivial exceptional cycles have length four, while the
    // four cross terms in H² have determinant 1+8T+30T²+72T³+81T⁴.
    let cross = [1, 8, 30, 72, 81];
    let mut traces = [0i128; 9];
    for k in 1..=8 {
        traces[k] = -((k as i128) * cross.get(k).copied().unwrap_or(0)
            + (1..k.min(cross.len()))
                .map(|j| cross[j] * traces[k - j])
                .sum::<i128>());
    }
    let surface = |n: usize| {
        let q = 3i128.pow(n as u32);
        1 + q * q + (2 + 16 * i128::from(n.is_multiple_of(4))) * q + traces[n]
    };
    // An unordered pair, a degree-two orbit, or a tangent direction.
    let pairs =
        |n: usize| (surface(n).pow(2) + surface(2 * n)) / 2 + 3i128.pow(n as u32) * surface(n);
    for (i, reduced) in [0, 1, 1, 2].into_iter().enumerate() {
        let n = i + 1;
        let orbits = |counts: &dyn Fn(usize) -> i128| {
            let previous = if reduced == 0 { 0 } else { counts(reduced) };
            (counts(n) - previous) / n as i128
        };
        assert_eq!(
            output[2 + 3 * i..5 + 3 * i],
            [orbits(&surface), orbits(&surface), orbits(&pairs)]
        );
    }
}

#[test]
fn wide_fixed_point_counts_are_divided_before_narrowing() {
    let input = format!("{} rvK3FrobHilbOrbits[Moving,C,1,41];", geometry());
    assert_eq!(
        stdout(interpret(&input)),
        "rvK3FrobHilbOrbits[Moving, C, 1, 41] = 32445840603149105116422276603348101760\n"
    );
    let fixed = format!("{} rvK3FrobHilbFixed[Moving,C,1,41];", geometry());
    assert!(diagnostic(interpret(&fixed)).contains("unsigned 128-bit range"));
}

#[test]
fn invalid_types_periods_and_bounds_are_reported() {
    for (source, expected) in [
        ("show rvWittSpectrumVerschiebung[0](zero[RVWittSpectrum]);", "must be positive"),
        ("rvWittSpectrumNecklace[0](zero[RVWittSpectrum]);", "must be positive"),
        ("show rvWittSpectrumVerschiebung[4](S.hi);", "rvWittSpectrumVerschiebung"),
        ("rvWittSpectrumNecklace[1](S.hi);", "RVWittSpectrum"),
        ("show rvWittSpectrumVerschiebung[257](rvWittSpectrumLine[1]);", "determinant degree limit"),
        ("rvWittSpectrumNecklace[400](rvWittSpectrumLine[2]);", "signed 128-bit range"),
        ("rvWittSpectrumNecklace[170141183460469231731687303715884105727](rvWittSpectrum{[1,1,1]});", "trial division limit"),
    ] {
        let actual = diagnostic(interpret(source));
        assert!(actual.contains(expected), "{actual}");
    }
    for (query, expected) in [
        (
            "rvK3FrobHilbOrbits[Moving,C,0,0];",
            "primitive period must be positive",
        ),
        ("rvK3FrobHilbOrbits[Missing,C,0,1];", "automorphism"),
        (
            "show rvK3FrobHilbSpectrum[Moving,C,2,1];",
            "determinant degree limit",
        ),
        (
            "rvK3FrobHilbOrbits[Moving,C,2,340282366920938463463374607431768211455];",
            "512 bits",
        ),
        (
            "rvK3FrobHilbOrbits[Moving,C,1,50];",
            "unsigned 128-bit range",
        ),
    ] {
        let actual = diagnostic(interpret(&format!("{} {query}", geometry())));
        assert!(actual.contains(expected), "{actual}");
    }
}
