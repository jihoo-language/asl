mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/hyperelliptic-kummer-artin-tate.asl");
const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{o,z,z,z,z,o};
";

#[test]
fn artin_invariants_and_exact_brauer_orders() {
    let output = stdout(example("hyperelliptic-kummer-artin-tate"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 7);
    assert_eq!(lines[4], full_rank("N", 4, 81, 2, 81, 1));
    assert_eq!(lines[5], full_rank("N", 8, 6561, 2, 81, 81));
    assert_eq!(lines[6], full_rank("S", 1, 25, 1, 25, 1));
    assert!(lines[2].ends_with("Artin invariant=2, geometric NS discriminant=-81"));
    assert!(lines[3].ends_with("Artin invariant=2, geometric NS discriminant=-81"));
    assert!(!lines[0].contains("Artin invariant"));
    assert!(!lines[1].contains("Artin invariant"));
}

fn full_rank(
    name: &str,
    r: usize,
    q: u128,
    sigma: usize,
    discriminant: u128,
    brauer: u128,
) -> String {
    format!(
        "rvKummerK3ArtinTate[{name}, {r}] = field order={q}, rank=22, geometric rank=22, finite NS discriminant times Brauer order=-{q}, finite NS discriminant square class=[-{q}], Artin invariant={sigma}, geometric NS discriminant=-{discriminant}, Brauer order={brauer}"
    )
}

#[test]
fn direct_p2_leading_term_agrees_on_ordinary_and_nonscalar_rank_jumps() {
    let output = stdout(example("hyperelliptic-kummer-artin-tate"));
    let lines: Vec<_> = output.lines().collect();
    // These are Q times the cross quartics in x=QT. The two remaining
    // abelian eigenvalues give (1-x)^2 in every case.
    let cases = [
        ([3, 8, 11, 8, 3], vec![1, 1, 1, 1, 3, 3, 3, 3], 1, 20),
        ([27, -64, 74, -64, 27], vec![1, 1, 1, 1, 3, 3, 3, 3], 3, 20),
        ([3, 6, 6, 6, 3], vec![1, 1, 2, 4, 4, 4], 1, 22),
        ([9, 0, -18, 0, 9], vec![1, 1, 2, 4, 4, 4], 2, 22),
        ([81, -324, 486, -324, 81], vec![1, 1, 2, 4, 4, 4], 4, 22),
    ];
    for (line, (cross, orbits, r, geometric)) in lines.iter().zip(cases) {
        let mut p2 = product(&cross, &[1, -2, 1]);
        for d in orbits {
            let g = (1..=d).rev().find(|g| d % g == 0 && r % g == 0).unwrap();
            let mut factor = vec![0; d / g + 1];
            factor[0] = 1;
            factor[d / g] = -1;
            for _ in 0..g {
                p2 = product(&p2, &factor);
            }
        }
        let (rank, leading) = at_one(p2);
        let signed = if rank.is_multiple_of(2) {
            -leading
        } else {
            leading
        };
        assert!(
            line.contains(&format!("rank={rank}, geometric rank={geometric}")),
            "{line}"
        );
        assert!(line.contains(&format!("finite NS discriminant times Brauer order={signed}, finite NS discriminant square class=[{signed}]")), "{line}");
    }
}

fn at_one(mut polynomial: Vec<i128>) -> (usize, i128) {
    let mut multiplicity = 0;
    loop {
        let value = polynomial.iter().sum();
        if value != 0 {
            return (multiplicity, value);
        }
        // Division by 1-x in ascending coefficient order.
        for i in 1..polynomial.len() {
            polynomial[i] += polynomial[i - 1];
        }
        assert_eq!(polynomial.pop(), Some(0));
        multiplicity += 1;
    }
}

fn product(a: &[i128], b: &[i128]) -> Vec<i128> {
    let mut result = vec![0; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            result[i + j] += x * y;
        }
    }
    result
}

#[test]
fn field_orders_beyond_the_point_count_range_are_supported() {
    let output = stdout(interpret(&format!("{EXAMPLE}\nrvKummerK3ArtinTate[S,20];")));
    let expected = full_rank("S", 20, 25u128.pow(20), 1, 25, 5u128.pow(38));
    assert_eq!(output.lines().last(), Some(expected.as_str()));
}

#[test]
fn rejects_zero_excessive_degrees_and_inapplicable_curves() {
    for (source, expected) in [
        (
            "rvKummerK3ArtinTate[C,0];",
            "extension degree must be positive",
        ),
        (
            "rvKummerK3ArtinTate[C,1000];",
            "exceeds the checked 128-bit arithmetic range",
        ),
        (
            "rvKummerK3ArtinTate[C,340282366920938463463374607431768211455];",
            "exceeds the checked 128-bit arithmetic range",
        ),
        (
            "rvKummerK3ArtinTate[Missing,1];",
            "unknown residue hyperelliptic curve",
        ),
        (
            "rvhyperelliptic E[F]{o,o,z,o}; rvKummerK3ArtinTate[E,1];",
            "monic quintic",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F3}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
