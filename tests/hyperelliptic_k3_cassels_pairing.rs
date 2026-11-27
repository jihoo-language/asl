mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/hyperelliptic-k3-cassels-pairing.asl");
const B_AND_SUM: &str = r"
rvk3cassels Second[Bcover]{
  0{[4],[12]}{[5],[3]}{[1]},
  1{}{[6],[3]}{[1]},
  2{[0],[8]}{[7],[21]}{[1]}
};
rvk3cassels Third[Sum]{
  0{[0],[3],[9]}{[9],[4],[22]}{[1],[1]},
  1{[2]}{[1],[17]}{[1]},
  2{[10]}{}{[1]}
};
";

fn setup() -> &'static str {
    EXAMPLE.split_once("rvK3CasselsTangents[T]").unwrap().0
}

fn reports(output: &str) -> impl Iterator<Item = &str> {
    output
        .lines()
        .filter(|line| line.starts_with("rvK3CasselsPair["))
}

fn sign(line: &str) -> i32 {
    line.split_once(" = sign=")
        .unwrap()
        .1
        .split_once(',')
        .unwrap()
        .0
        .parse()
        .unwrap()
}

fn numbers(text: &str) -> Vec<i64> {
    text.split(|c: char| !c.is_ascii_digit() && c != '-')
        .filter(|part| !part.is_empty())
        .map(|part| part.parse().unwrap())
        .collect()
}

fn character(angular: i64) -> i64 {
    assert_ne!(angular.rem_euclid(23), 0);
    if (1..23).any(|r| (r * r - angular).rem_euclid(23) == 0) {
        1
    } else {
        -1
    }
}

fn hilbert((v, a): (i64, i64), (w, b): (i64, i64)) -> i64 {
    // q=23, so -1 is nonsquare. This checks the reported local signs
    // of the actual tagged valuations and leading units independently.
    (-1i64).pow((v * w).rem_euclid(2) as u32)
        * character(a).pow(w.rem_euclid(2) as u32)
        * character(b).pow(v.rem_euclid(2) as u32)
}

#[test]
fn exact_local_hilbert_symbols_match_the_independent_laurent_oracle() {
    let output = stdout(example("hyperelliptic-k3-cassels-pairing"));
    assert!(output.contains("omitted indices [0, 1, 2]="));
    assert_eq!(reports(&output).map(sign).collect::<Vec<_>>(), [-1, 1, 1]);
    let row = reports(&output).next().unwrap();
    assert!(row.contains("sign=-1, value=1/2 in Q/Z"));
    let local = row
        .split_once("local factors=[")
        .unwrap()
        .1
        .trim_end_matches(']')
        .split("; ")
        .collect::<Vec<_>>();
    let places = [0, 2, 10, 15, 16, 20, 21].map(|x| format!("t=[{x}]"));
    let mut all_products = Vec::new();
    // The independent calculation with raw Q_i gives these local products.
    // Primitive normalization scales its three tangents by 1/5, 1/3, 19/f.
    let raw_products = [1, -1, -1, 1, 1, 1, 1, -1];
    for (i, factor) in local.iter().enumerate() {
        let label = places.get(i).map_or("infinity", String::as_str);
        assert!(factor.starts_with(label), "{factor}");
        let (left, tail) = factor
            .split_once("tangent (valuation, leading unit)=")
            .unwrap()
            .1
            .split_once(", right (valuation, leading unit)=")
            .unwrap();
        let (right, tail) = tail.split_once(", symbols=").unwrap();
        let (symbols, product) = tail.split_once(", product=").unwrap();
        let left = numbers(left);
        let right = numbers(right);
        let symbols = numbers(symbols);
        assert_eq!((left.len(), right.len(), symbols.len()), (6, 6, 3));
        for j in 0..3 {
            assert_eq!(
                symbols[j],
                hilbert(
                    (left[2 * j], left[2 * j + 1]),
                    (right[2 * j], right[2 * j + 1])
                )
            );
        }
        let product: i64 = product.parse().unwrap();
        assert_eq!(symbols.iter().product::<i64>(), product);
        let third_scalar = if i < 7 {
            let t = numbers(label)[0];
            let (order, leading) = if t == 0 {
                (1, -10)
            } else if t == 10 {
                (1, 10)
            } else {
                (0, t * (t - 10))
            };
            let angular = (1..23)
                .find(|&r| (r * leading - 19).rem_euclid(23) == 0)
                .unwrap();
            (-order, angular)
        } else {
            (2, 19)
        };
        let normalization = [(0, 14), (0, 8), third_scalar]
            .into_iter()
            .enumerate()
            .map(|(j, scale)| hilbert(scale, (right[2 * j], right[2 * j + 1])))
            .product::<i64>();
        assert_eq!(product, raw_products[i] * normalization, "{label}");
        if i == 7 {
            assert!(factor.contains("x-chart s^-2*[[1]]"));
            // Direct Laurent expansion with x=s^-6 and the normalized
            // tangents has cancellations through the negative orders.
            assert_eq!(left, [-5, 13, -4, 13, -4, 10]);
            assert_eq!(right, [-2, 5, -2, 1, -2, 5]);
            assert_eq!(symbols, [-1, 1, 1]);
        }
        all_products.push(product);
    }
    assert_eq!(all_products, [-1, 1, -1, 1, 1, 1, 1, -1]);
}

#[test]
fn bilinearity_and_the_global_mordell_weil_radical_are_verified() {
    let source = format!(
        r"{}{B_AND_SUM}
        let c=add(a,b);
        rvK3CasselsPair[T](a);
        rvK3CasselsPair[T](b);
        rvK3CasselsPair[T](c);
        rvK3CasselsPair[Second](a);
        rvK3CasselsPair[Second](b);
        rvK3CasselsPair[Second](c);
        rvK3CasselsPair[Third](a);
        rvK3CasselsPair[Third](b);
        rvK3CasselsPair[Third](c);
        let d0=rvK3PencilDescent2(rvK3PencilTorsion2At[P,0]);
        let d1=rvK3PencilDescent2(rvK3PencilTorsion2At[P,1]);
        let dq=rvK3PencilDescent2(q);
        rvK3CasselsPair[T](d0);
        rvK3CasselsPair[T](d1);
        rvK3CasselsPair[T](delta);
        rvK3CasselsPair[T](dq);
        rvK3CasselsPair[Second](d0);
        rvK3CasselsPair[Second](d1);
        rvK3CasselsPair[Second](delta);
        rvK3CasselsPair[Second](dq);
        rvK3CasselsPair[T](add(b,dq));
        rvk3cassels Global[RationalCover]{{
          0{{[0],[6],[18]}}{{[15],[19]}}{{[6]}},
          1{{[2],[1]}}{{[15],[19]}}{{[6]}},
          2{{[2],[1]}}{{[0],[6],[18]}}{{[6]}}
        }};
        rvK3CasselsPair[Global](a);
        rvK3CasselsPair[Global](b);
        rvk3pencilcover Trivial[P](zero[RVK3PencilKummer2(P)]);
        rvk3cassels Infinity[Trivial]{{0{{[1]}}{{[1]}}{{}},1{{[1]}}{{[1]}}{{}},2{{[1]}}{{[1]}}{{}}}};
        rvK3CasselsPair[Infinity](b);
    ",
        setup()
    );
    let result = stdout(interpret(&source));
    let signs = reports(&result).map(sign).collect::<Vec<_>>();
    assert_eq!(&signs[..9], [1, -1, -1, -1, 1, -1, -1, -1, 1]);
    assert_eq!(&signs[9..], [1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1]);
}

#[test]
fn homogeneous_rescaling_and_reordering_produce_identical_local_certificates() {
    let source = format!(
        r"{}
        rvk3cassels Scaled[A]{{
          2{{[0],[19],[20],[1],[17]}}{{[0],[16],[18],[1],[17]}}{{[0],[3],[2]}},
          0{{[5],[5]}}{{[0],[13],[17],[4]}}{{[1],[1]}},
          1{{[22],[0],[6]}}{{[0],[8],[16],[6]}}{{[2]}}
        }};
        rvK3CasselsTangents[T]; rvK3CasselsTangents[Scaled];
        rvK3CasselsPair[T](b); rvK3CasselsPair[Scaled](b);
        rvk3cassels Opposite[A]{{
          0{{[5]}}{{[0],[-13],[-4]}}{{[1]}},
          1{{[11],[0],[3]}}{{[0],[4],[8],[3]}}{{[1]}},
          2{{[0],[19],[20],[1],[17]}}{{[0],[16],[18],[1],[17]}}{{[0],[3],[2]}}
        }};
        rvK3CasselsPair[Opposite](b);
    ",
        setup()
    );
    let output = stdout(interpret(&source));
    let tangents: Vec<_> = output
        .lines()
        .filter(|line| line.starts_with("rvK3CasselsTangents["))
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(tangents.len(), 2);
    assert_eq!(tangents[0], tangents[1]);
    let pairings: Vec<_> = reports(&output).collect();
    assert_eq!(
        pairings[0].split_once(" = ").unwrap().1,
        pairings[1].split_once(" = ").unwrap().1
    );
    assert_eq!(sign(pairings[2]), -1);
}

#[test]
fn odd_residue_hilbert_symbols_use_the_actual_extension_field() {
    let output = stdout(interpret(
        r"
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
        rvk3pencilcover T[P](zero[RVK3PencilKummer2(P)]);
        rvk3cassels Local[T]{0{[1,0]}{[1,0]}{},1{[1,0]}{[1,0]}{},2{[1,0]}{[1,0]}{}};
        rvK3CasselsPair[Local](rvK3PencilDescent2(rvK3PencilTorsion2At[P,1]));
        ",
    ));
    let row = reports(&output).next().unwrap();
    assert_eq!(sign(row), 1);
    assert!(row.contains("t=[0, 1]"));
    assert!(row.contains("infinity: x-chart s^-2*[[1, 0]]"));
}

#[test]
fn invalid_or_incomplete_conic_witnesses_and_wrong_classes_are_rejected() {
    let source = format!(
        r"{}
        rvk3pencil N[B,14,0];
        let localFail=rvK3PencilKummer2[P]{{[5]}}{{[5]}}{{[1]}};
        rvk3pencilcover Fail[P](localFail);
    ",
        setup()
    );
    for (operation, expected) in [
        (
            "rvk3cassels C[A]{};",
            "each omitted index 0, 1 and 2 exactly once",
        ),
        ("rvk3cassels C[A]{3{[1]}{[1]}{[1]}};", "each omitted index"),
        (
            "rvk3cassels C[A]{0{[5]}{[0],[13],[4]}{[1]},0{[5]}{[0],[13],[4]}{[1]}};",
            "each omitted index",
        ),
        (
            "rvk3cassels C[A]{0{}{}{}};",
            "must be a nonzero projective point",
        ),
        (
            "rvk3cassels C[A]{0{[1]}{[1]}{[1]}};",
            "does not lie on its conic",
        ),
        (
            "rvk3cassels C[Fail]{};",
            "first covering must be everywhere locally soluble",
        ),
        (
            "rvK3CasselsPair[T](localFail);",
            "second class must be everywhere locally soluble",
        ),
        (
            "rvK3CasselsPair[T](zero[RVK3PencilKummer2(N)]);",
            "expected type RVK3PencilKummer2(P)",
        ),
    ] {
        let failure = diagnostic(interpret(&format!("{source}\n{operation}")));
        assert!(failure.contains(expected), "{operation}: {failure}");
    }
}
