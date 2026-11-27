mod common;

use std::collections::BTreeSet;

use common::{diagnostic, example, interpret, stdout};

const F13: &str = include_str!("../examples/hyperelliptic-k3-pencil-regulator.asl");
const F23: &str = include_str!("../examples/hyperelliptic-k3-pencil-two-covering.asl");

fn f13() -> &'static str {
    F13.split_once("rvK3PencilSectionHeight[P]").unwrap().0
}

fn f23() -> &'static str {
    F23.split_once("rvk3pencilmw2").unwrap().0
}

#[test]
fn true_shioda_heights_and_pairing_certify_rank_two_f23() {
    let output = stdout(interpret(&format!(
        r"{}
        let torsion=rvK3PencilTorsion2At[P,0];
        rvK3PencilSectionHeight[P](p);
        rvK3PencilSectionHeight[P](q);
        rvK3PencilSectionHeight[P](add(p,q));
        rvK3PencilSectionHeight[P](add(p,neg(q)));
        rvK3PencilSectionHeight[P](add(p,p));
        rvK3PencilSectionHeight[P](torsion);
        rvK3PencilSectionHeight[P](add(p,torsion));
        rvK3PencilSectionPair[P](p,q);
        rvK3PencilSectionPair[P](add(p,p),q);
        rvK3PencilSectionPair[P](p,torsion);
        rvK3PencilRegulator[P]{{p,q}};
        rvK3PencilRegulator[P]{{add(p,q),add(q,torsion)}};
        ",
        f23(),
    )));
    let mut lines = output.lines();
    for height in [1, 2, 1, 5, 4, 0, 1] {
        let line = lines.next().unwrap();
        assert!(line.starts_with(&format!("rvK3PencilSectionHeight[P] = {height}, component multiplier=4, minimal discriminant degree=24")), "{line}");
    }
    for pair in [-1, -2, 0] {
        let line = lines.next().unwrap();
        assert!(
            line.starts_with(&format!(
                "rvK3PencilSectionPair[P] = {pair}, component multiplier=4"
            )),
            "{line}"
        );
    }
    for expected_gram in ["[[1, -1], [-1, 2]]", "[[1, 1], [1, 2]]"] {
        let line = lines.next().unwrap();
        assert!(line.contains(&format!("Gram={expected_gram}, subgroup determinant=1, D^r*det=16, mod-2 image dimension=4/4")), "{line}");
        assert!(line.contains("arithmetic index candidates=[1], good-fiber sieve=[], exact index=1, regulator=1, Sha order=4, Sha=(Z/2)^2"), "{line}");
    }
    assert_eq!(lines.next(), None);
}

const PRIME: i64 = 13;
type FP = Option<(i64, i64)>;

fn field(value: i64) -> i64 {
    value.rem_euclid(PRIME)
}

fn poly_at(coefficients: &[i64], t: i64) -> i64 {
    coefficients.iter().rev().fold(0, |v, &c| field(t * v + c))
}

fn add(x: FP, y: FP, a2: i64, a4: i64) -> FP {
    let (Some((x1, y1)), Some((x2, y2))) = (x, y) else {
        return x.or(y);
    };
    let (numerator, denominator) = if x1 == x2 {
        if field(y1 + y2) == 0 {
            return None;
        }
        (3 * x1 * x1 + 2 * a2 * x1 + a4, 2 * y1)
    } else {
        (y2 - y1, x2 - x1)
    };
    let slope = (0..PRIME)
        .find(|&u| field(u * denominator - numerator) == 0)
        .unwrap();
    let x3 = field(slope * slope - x1 - x2 - a2);
    Some((x3, field(slope * (x1 - x3) - y1)))
}

fn local_kernel(t: i64) -> BTreeSet<(i64, i64)> {
    let roots = [[0, 1, 9, 6], [0, 6, 8, 10], [0, 8, 6, 9]].map(|e| poly_at(&e, t));
    let [a, b, c] = roots;
    assert!(a != b && b != c && c != a);
    let a2 = field(-a - b - c);
    let a4 = field(a * b + b * c + c * a);
    let add = |x, y| add(x, y, a2, a4);
    let mult = |n, x| (0..n).fold(None, |y, _| add(y, x));
    let mut points = vec![None];
    for x in 0..PRIME {
        for y in 0..PRIME {
            if field(y * y - (x - a) * (x - b) * (x - c)) == 0 {
                points.push(Some((x, y)));
            }
        }
    }
    assert_eq!(points.len(), 12);
    let triples: BTreeSet<_> = points.into_iter().map(|x| mult(3, x)).collect();
    let p = Some((poly_at(&[0, 1, 0, 6], t), poly_at(&[0, 0, 7, 6, 3], t)));
    let q = Some((poly_at(&[12, 0, 12, 6], t), poly_at(&[8, 8, 10, 3, 6], t)));
    if t == 1 {
        // At this witness place the section 3P has denominator poles and
        // reduces to O, a case tested in the unsaturated run as well.
        assert_ne!(p, None);
        assert_eq!(mult(3, p), None);
    }
    (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .filter(|&(i, j)| triples.contains(&add(mult(i, p), mult(j, q))))
        .collect()
}

#[test]
fn independently_enumerated_good_fibers_remove_the_f13_index_three() {
    let first = local_kernel(1);
    let second = local_kernel(7);
    assert_eq!(first, BTreeSet::from([(0, 0), (1, 2), (2, 1)]));
    assert_eq!(
        first.intersection(&second).copied().collect::<Vec<_>>(),
        [(0, 0)]
    );

    let output = stdout(example("hyperelliptic-k3-pencil-regulator"));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].starts_with("rvK3PencilSectionHeight[P] = 1, component multiplier=2"));
    assert!(lines[1].starts_with("rvK3PencilSectionHeight[P] = 5/2, component multiplier=2"));
    assert!(lines[2].starts_with("rvK3PencilSectionPair[P] = 1/2, component multiplier=2"));
    assert!(lines[3].contains("Gram=[[1, 1/2], [1/2, 5/2]], subgroup determinant=9/4, D^r*det=9, mod-2 image dimension=4/4, arithmetic index candidates=[1, 3]"));
    assert!(lines[3].contains("ell=3, reductions (parameter, group order, common kernel dimension)=[([1], 12, 1), ([7], 12, 0)], excluded"));
    assert!(lines[3].contains("exact index=1, regulator=9/4, Sha order=1, Sha=0"));
}

#[test]
fn unsaturated_multiples_keep_only_certified_conclusions() {
    let output = stdout(interpret(&format!(
        "{} let triplePoint=add(add(p,p),p); rvK3PencilRegulator[P]{{triplePoint,q}};",
        f13(),
    )));
    assert!(output.contains("subgroup determinant=81/4, D^r*det=81, mod-2 image dimension=4/4, arithmetic index candidates=[3, 9]"), "{output}");
    assert!(
        output.contains("ell=3") && output.contains("inconclusive (kernel dimension 1)"),
        "{output}"
    );
    assert!(
        output.contains("unresolved (index, regulator, Sha order)=[(3, 9/4, 1), (9, 1/4, 9)]"),
        "{output}"
    );

    let output = stdout(interpret(&format!(
        "{} rvK3PencilRegulator[P]{{add(p,p),q}};",
        f23()
    )));
    assert!(output.contains("subgroup determinant=4, D^r*det=64, mod-2 image dimension=3/4, arithmetic index candidates=[2, 4, 8]"), "{output}");
    assert!(
        output.contains(
            "unresolved (index, regulator, Sha order)=[(2, 1, 4), (4, 1/4, 16), (8, 1/16, 64)]"
        ),
        "{output}"
    );
}

#[test]
fn dependent_rank_and_mixed_nominal_inputs_are_rejected() {
    for (case, expected) in [
        ("rvK3PencilRegulator[P]{p};", "needs 2 sections"),
        (
            "rvK3PencilRegulator[P]{p,add(p,rvK3PencilTorsion2At[P,0])};",
            "dependent modulo torsion",
        ),
        (
            "rvK3PencilSectionHeight[P](rvK3PencilDescent2(p));",
            "expected type RVK3PencilSection(P)",
        ),
        (
            "rvk3pencil Q[B,1,0]; rvK3PencilSectionPair[P](p,zero[RVK3PencilSection(Q)]);",
            "expected type RVK3PencilSection(P)",
        ),
    ] {
        let message = diagnostic(interpret(&format!("{}\n{case}", f13())));
        assert!(message.contains(expected), "{case}: {message}");
    }
}
