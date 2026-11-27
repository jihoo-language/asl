mod common;

use std::collections::BTreeSet;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F23: &str = include_str!("../examples/hyperelliptic-k3-pencil-selmer.asl");
const ROOTS: [[i64; 5]; 3] = [[0, 12, 1, 20, 7], [0, 18, 20, 21, 1], [0, 21, 14, 9, 8]];
const PRIME: i64 = 23;

fn label(angular: i64, order: usize) -> u8 {
    let square = (1..PRIME).any(|i| (i * i - angular).rem_euclid(PRIME) == 0);
    (order % 2) as u8 + 2 * u8::from(!square)
}

fn symbol(x: &[i64], y: &[i64]) -> Option<u8> {
    x.iter()
        .zip(y)
        .enumerate()
        .find(|(_, (x, y))| (*x - *y).rem_euclid(PRIME) != 0)
        .map(|(order, (x, y))| label(x - y, order))
}

fn pair(symbols: [u8; 3]) -> u8 {
    symbols[0] + 4 * symbols[1]
}

fn exhaustive_image(parameter: Option<i64>) -> BTreeSet<u8> {
    let binomial = [
        [1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [1, 3, 3, 1, 0],
        [1, 4, 6, 4, 1],
    ];
    let chart = ROOTS.map(|root| {
        std::array::from_fn::<_, 5, _>(|j| {
            parameter.map_or(root[4 - j], |a| {
                (j..5)
                    .map(|k| root[k] * binomial[k][j] * a.pow((k - j) as u32))
                    .sum::<i64>()
                    .rem_euclid(PRIME)
            })
        })
    });
    let depth = [(0, 1), (0, 2), (1, 2)]
        .into_iter()
        .map(|(i, j)| (0..5).find(|&k| chart[i][k] != chart[j][k]).unwrap() + 1)
        .max()
        .unwrap();
    let mut image = BTreeSet::from([0]);
    // Enumerate the whole uniformizer quotient, including branches that
    // never meet a root. This does not share the implementation's trie.
    for mut index in 0..PRIME.pow(depth as u32) {
        let x: Vec<_> = (0..depth)
            .map(|_| {
                let digit = index % PRIME;
                index /= PRIME;
                digit
            })
            .collect();
        let symbols = chart.each_ref().map(|root| symbol(&x, root));
        if let [Some(a), Some(b), Some(c)] = symbols {
            if a ^ b ^ c == 0 {
                image.insert(pair([a, b, c]));
            }
        } else {
            // At this depth the ball of a root contains only that root.
            let root = symbols.iter().position(Option::is_none).unwrap();
            let mut symbols = symbols.map(|entry| entry.unwrap_or(0));
            symbols[root] = symbols[0] ^ symbols[1] ^ symbols[2];
            image.insert(pair(symbols));
        }
    }
    image
}

fn displayed_class(mask: usize, places: &[i64]) -> String {
    let width = places.len() + 1;
    let left = mask % (1 << width);
    let right = mask >> width;
    let coordinates = [left, right, left ^ right].map(|coordinate| {
        let mut support: Vec<_> = places
            .iter()
            .enumerate()
            .filter(|&(i, _)| coordinate & (2 << i) != 0)
            .map(|(_, &a)| [[(-a).rem_euclid(PRIME)], [1]])
            .collect();
        support.sort_unstable();
        let kind = if coordinate & 1 == 0 {
            "square"
        } else {
            "nonsquare"
        };
        format!("({kind}, {support:?})")
    });
    format!("RVK3PencilKummer2(P)[{}]", coordinates.join(", "))
}

#[test]
fn f23_local_images_and_all_selmer_classes_agree_with_independent_exhaustion() {
    let setup = F23.split_once("rvK3PencilSelmer[P]").unwrap().0;
    let places: Vec<_> = (0..PRIME)
        .filter(|&a| {
            let roots = ROOTS.map(|root| {
                root.into_iter()
                    .rev()
                    .fold(0, |acc, c| (acc * a + c).rem_euclid(PRIME))
            });
            roots[0] == roots[1] || roots[0] == roots[2] || roots[1] == roots[2]
        })
        .collect();
    assert_eq!(places, [0, 2, 10, 15, 16, 20, 21]);
    let local: Vec<_> = places
        .iter()
        .copied()
        .map(Some)
        .chain([None])
        .map(|place| {
            let image = exhaustive_image(place);
            assert_eq!(image.len(), 4);
            let mut symbols = vec![2];
            symbols.extend(places.iter().map(|&a| match place {
                Some(t) if t != a => label(t - a, 0),
                _ => 1,
            }));
            (image, symbols)
        })
        .collect();
    let width = places.len() + 1;
    let expected: BTreeSet<_> = (0..(1 << (2 * width)))
        .filter(|&mask| {
            local.iter().all(|(image, symbols)| {
                let chart = [0, width].map(|offset| {
                    symbols.iter().enumerate().fold(0, |acc, (j, &symbol)| {
                        acc ^ if mask & (1 << (offset + j)) != 0 {
                            symbol
                        } else {
                            0
                        }
                    })
                });
                image.contains(&(chart[0] + 4 * chart[1]))
            })
        })
        .map(|mask| displayed_class(mask, &places))
        .collect();
    assert_eq!(expected.len(), 64);

    let mut source = format!("{setup} rvK3PencilTorsion2[P]; rvK3PencilSelmer[P];");
    for index in 0..6 {
        writeln!(source, "let b{index}=rvK3PencilSelmerBasis[P,{index}];").unwrap();
    }
    for mask in 0..64 {
        let word = (0..6)
            .filter(|&i| mask & (1 << i) != 0)
            .fold("zero[RVK3PencilKummer2(P)]".to_owned(), |word, i| {
                format!("add({word},b{i})")
            });
        writeln!(
            source,
            "assert rvK3PencilSelmerContains({word})==S.hi; show {word};"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    assert_eq!(
        lines.next().unwrap(),
        format!(
            "rvK3PencilTorsion2[P] = order=4, affine x-roots={:?}",
            ROOTS.map(|root| root.map(|c| [c]))
        )
    );
    let finite: Vec<_> = places
        .iter()
        .zip(&local)
        .map(|(&a, (image, _))| ([a], image.iter().copied().collect::<Vec<_>>()))
        .collect();
    let infinity = local.last().unwrap().0.iter().copied().collect::<Vec<_>>();
    assert_eq!(
        lines.next().unwrap(),
        format!("rvK3PencilSelmer[P] = dim Selmer=6, Mordell-Weil rank=2, dim Sha[2]=2, finite local images (parameter, codes)={finite:?}, infinity codes={infinity:?}")
    );
    assert_eq!(lines.map(str::to_owned).collect::<BTreeSet<_>>(), expected);
    assert!(stdout(example("hyperelliptic-k3-pencil-selmer")).contains("dim Sha[2]=2"));
}

#[test]
fn rational_sections_and_two_torsion_map_into_the_space_over_f7_f13_and_f9() {
    let f7 = include_str!("../examples/hyperelliptic-k3-pencil-bsd.asl")
        .split_once("rvK3FiberEuler")
        .unwrap()
        .0;
    let f13 = include_str!("../examples/hyperelliptic-k3-pencil-two-descent.asl")
        .split_once("rvK3PencilTorsion2[P]")
        .unwrap()
        .0;
    let f9 = r"
        rvunram F[3,1]{1,0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let i=rvUnramElem[F]{0,1};
        rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
        let r0=rvJacPoint[C](z,z);
        let r1=rvJacPoint[C](o,z);
        let r2=rvJacPoint[C](neg(o),z);
        let r3=rvJacPoint[C](i,z);
        rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
        rvk3pencil P[B,0,1];
    ";
    for (setup, report) in [
        (f7, "dim Selmer=3, Mordell-Weil rank=1, dim Sha[2]=0"),
        (f13, "dim Selmer=4, Mordell-Weil rank=2, dim Sha[2]=0"),
        (
            f9,
            "dim Selmer=4, Mordell-Weil rank and dim Sha[2] unavailable in characteristic 3",
        ),
    ] {
        let mut source = format!(
            r"{setup}
            rvK3PencilSelmer[P];
            let t0=rvK3PencilDescent2(rvK3PencilTorsion2At[P,0]);
            let t1=rvK3PencilDescent2(rvK3PencilTorsion2At[P,1]);
            let b=rvK3PencilSelmerBasis[P,0];
            assert rvK3PencilSelmerContains(t0)==S.hi;
            assert rvK3PencilSelmerContains(t1)==S.hi;
            assert rvK3PencilSelmerContains(add(t0,t1))==S.hi;
            assert rvK3PencilSelmerContains(b)==S.hi;
            "
        );
        if setup == f13 {
            source.push_str(
                "let g=rvK3PencilSectionPoint[P]{[0],[1],[0],[6]}{[1]}{[0],[0],[7],[6],[3]}{[1]};\n\
                assert rvK3PencilSelmerContains(rvK3PencilDescent2(g))==S.hi;\n\
                assert rvK3PencilDivisible2(g)==S.lo;",
            );
        }
        assert!(stdout(interpret(&source)).contains(report));
    }
}

#[test]
fn invalid_index_operand_and_nominal_mixture_report_errors() {
    let setup = F23.split_once("rvK3PencilSelmer[P]").unwrap().0;
    for (operation, needle) in [
        ("show rvK3PencilSelmerBasis[P,6];", "outside dimension 6"),
        (
            "show rvK3PencilSelmerContains(rvK3PencilTorsion2At[P,0]);",
            "an RVK3PencilKummer2 class",
        ),
        (
            "rvk3pencil Q[B,14,0]; show add(rvK3PencilSelmerBasis[P,0],rvK3PencilSelmerBasis[Q,0]);",
            "expected type RVK3PencilKummer2(P)",
        ),
    ] {
        let failure = diagnostic(interpret(&format!("{setup}\n{operation}")));
        assert!(failure.contains(needle), "{failure}");
    }
}
