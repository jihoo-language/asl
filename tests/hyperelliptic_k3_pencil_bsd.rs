mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup(prime: i64, extension: bool) -> String {
    let mut source = format!(
        r"
        rvunram K[{prime},1]{{0,1}};
        let zk=rvUnramElem[K]{{}};
        let ok=rvUnramElem[K]{{1}};
        let c1=rvUnramElem[K]{{24}};
        let c2=rvUnramElem[K]{{-50}};
        let c3=rvUnramElem[K]{{35}};
        let c4=rvUnramElem[K]{{-10}};
        rvhyperelliptic C[K]{{zk,c1,c2,c3,c4,ok}};
        "
    );
    let (curve, field) = if extension {
        source.push_str("rvunram F[5,1]{2,0,1}; rvjacbase L[C,F];");
        ("L", "F")
    } else {
        ("C", "K")
    };
    write!(
        source,
        r"
        let z=rvUnramElem[{field}]{{}};
        let o=rvUnramElem[{field}]{{1}};
        let two=rvUnramElem[{field}]{{2}};
        let three=rvUnramElem[{field}]{{3}};
        let r0=rvJacPoint[{curve}](z,z);
        let r1=rvJacPoint[{curve}](o,z);
        let r2=rvJacPoint[{curve}](two,z);
        let r3=rvJacPoint[{curve}](three,z);
        rvjactorsionbasis B[{curve},2](neg(o)){{r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)}};
        rvk3pencil P[B,0,1];
        "
    )
    .unwrap();
    source
}

#[test]
fn exact_special_values_agree_with_independent_tamagawa_and_artin_tate_products() {
    for (prime, extension, expected) in [
        (5, false, "q=5, rank=2, root number=1, L*=4, Tamagawa product=256, component orbit product=4, arithmetic trivial |det|=1024, signed NS discriminant times Brauer order=-80, BSD quotient=5/64, Taylor coefficients=[0, 0, 4, -4, 1]"),
        (5, true, "q=25, rank=4, root number=1, L*=1, Tamagawa product=1024, component orbit product=1, arithmetic trivial |det|=1024, signed NS discriminant times Brauer order=-25, BSD quotient=25/1024, Taylor coefficients=[0, 0, 0, 0, 1]"),
        (13, false, "q=13, rank=2, root number=1, L*=144/13, Tamagawa product=1024, component orbit product=1, arithmetic trivial |det|=1024, signed NS discriminant times Brauer order=-144, BSD quotient=9/64, Taylor coefficients=[0, 0, 144/13, -288/13, 232/13, -88/13, 1]"),
    ] {
        let mut source = setup(prime, extension);
        let (curve, field) = if extension { ("L", "F") } else { ("C", "K") };
        writeln!(source, "rvK3PencilBSD[P]; rvKummerK3ArtinTate[{curve},1];").unwrap();
        source.push_str("rvK3FiberEuler[P](rvK3PencilParam[P](z,o));");
        for a in 0..prime {
            for b in 0..if extension { prime } else { 1 } {
                let scalar = if extension { format!("{a},{b}") } else { a.to_string() };
                writeln!(
                    source,
                    "rvK3FiberEuler[P](rvK3PencilParam[P](o,rvUnramElem[{field}]{{{scalar}}}));"
                ).unwrap();
            }
        }
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        let bsd = lines.next().unwrap();
        assert_eq!(bsd, format!("rvK3PencilBSD[P] = {expected}"));
        let artin_tate = lines.next().unwrap();
        let product = value(artin_tate, "finite NS discriminant times Brauer order=");
        assert_eq!(product, value(bsd, "signed NS discriminant times Brauer order="));
        let mut tamagawa = 1;
        let mut orbit_product = 1;
        for relation in lines {
            tamagawa *= value(relation, "Tamagawa=");
            let cycles = relation.split_once("component cycles=[").unwrap().1.split(']').next().unwrap();
            if !cycles.is_empty() {
                orbit_product *= cycles.split(',').map(|n| n.trim().parse::<i128>().unwrap()).product::<i128>();
            }
        }
        assert_eq!(tamagawa, value(bsd, "Tamagawa product="));
        assert_eq!(orbit_product, value(bsd, "component orbit product="));
        let trivial = value(bsd, "arithmetic trivial |det|=");
        assert_eq!(trivial, tamagawa * orbit_product);
        let (numerator, denominator) = fraction(bsd, "BSD quotient=");
        assert_eq!(numerator * trivial, product.abs() * denominator);
    }
}

#[test]
fn runnable_odd_rank_example_checks_functional_sign_and_a_finite_height() {
    let output = stdout(example("hyperelliptic-k3-pencil-bsd"));
    assert!(output.contains("geometric Kodaira=I4"));
    assert!(output.contains("split=no"));
    assert!(output.contains("split=yes"));
    let bsd = output
        .lines()
        .find(|line| line.starts_with("rvK3PencilBSD"))
        .unwrap();
    assert_eq!(value(bsd, "rank=") % 2, 1);
    assert_eq!(value(bsd, "root number="), -1);
    assert_eq!(
        bsd,
        "rvK3PencilBSD[P] = q=7, rank=1, root number=-1, L*=32/7, Tamagawa product=512, component orbit product=2, arithmetic trivial |det|=1024, signed NS discriminant times Brauer order=-64, BSD quotient=1/16, Taylor coefficients=[0, 32/7, -48/7, 30/7, -1]"
    );
    assert!(output.ends_with("rvK3MordellHeight[M] = 1\n"));
    let artin = output
        .lines()
        .find(|line| line.starts_with("rvKummerK3ArtinTate"))
        .unwrap();
    assert_eq!(
        value(artin, "finite NS discriminant times Brauer order="),
        value(bsd, "signed NS discriminant times Brauer order=")
    );
}

#[test]
fn unknown_pencil_has_a_nominal_diagnostic() {
    let error = diagnostic(interpret("rvK3PencilBSD[Missing];"));
    assert!(error.contains("unknown K3 pencil"), "{error}");
}

fn tail<'a>(report: &'a str, label: &str) -> &'a str {
    report
        .split_once(label)
        .unwrap()
        .1
        .split(',')
        .next()
        .unwrap()
}

fn value(report: &str, label: &str) -> i128 {
    tail(report, label).parse().unwrap()
}

fn fraction(report: &str, label: &str) -> (i128, i128) {
    let text = tail(report, label);
    if let Some((numerator, denominator)) = text.split_once('/') {
        (numerator.parse().unwrap(), denominator.parse().unwrap())
    } else {
        (text.parse().unwrap(), 1)
    }
}
