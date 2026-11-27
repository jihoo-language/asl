mod common;

use common::{diagnostic, example, interpret, stdout};

const SOURCE: &str = include_str!("../examples/hyperelliptic-k3-iwasawa.asl");

fn base() -> &'static str {
    SOURCE
        .split_once("rvKummerK3Iwasawa[C,2,2,1000000]")
        .unwrap()
        .0
}

#[test]
fn f23_newton_breakpoints_and_absolute_sha_are_certified_from_base_data() {
    let output = stdout(example("hyperelliptic-k3-iwasawa"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 4);
    assert!(lines[0].contains("NS rank=20, geometric NS rank=20"));
    assert!(lines[0].contains("zero valuations=0, breakpoint k0=0, positive Newton valuations=[4: 2], determinant valuations through breakpoint=[8], mu=0, lambda=2, relative v_ell(Brauer) from breakpoint=2000000, relative v_ell(Brauer) from base=2000000"), "{}",lines[0]);
    assert!(lines[1].contains("zero valuations=2, breakpoint k0=0, positive Newton valuations=[], determinant valuations through breakpoint=[0], mu=0, lambda=0, relative v_ell(Sha) from breakpoint=0, relative v_ell(Sha) from base=0"));
    assert!(lines[2].contains("positive Newton valuations=[2: 2], determinant valuations through breakpoint=[4], mu=0, lambda=2, relative v_ell(Brauer) from breakpoint=2000000"));
    assert!(lines[3].contains("NS rank=18, geometric NS rank=20"));
    assert!(lines[3].contains("zero valuations=2, breakpoint k0=0, positive Newton valuations=[1: 2], determinant valuations through breakpoint=[2], mu=0, lambda=2, relative v_ell(Sha) from breakpoint=2000000"));
    assert!(lines[3].ends_with("certified Sha order over F_q=4, v_ell(Sha)=2000002"));
}

#[test]
fn finite_tower_ratios_and_the_certified_growth_law_agree() {
    let output = stdout(interpret(&format!(
        r"{}
        rvK3PencilArtinTateRatio[P,2,4]; rvKummerK3Iwasawa[C,2,2,1];
        rvK3PencilArtinTateRatio[P,2,6]; rvK3PencilShaIwasawa[P,2,3,1];
        rvK3PencilArtinTateRatio[P,4,12]; rvKummerK3Iwasawa[C,4,3,1];
        rvK3PencilShaPower[P,3]{{p,q}}; rvK3PencilShaIwasawa[P,3,5,0]{{p,q}};
    ",
        base()
    )));
    let lines: Vec<_> = output.lines().collect();
    for (i, prime) in [2u64, 3, 3].into_iter().enumerate() {
        let ratio: u64 = lines[2 * i]
            .split_once("order ratio=")
            .unwrap()
            .1
            .parse()
            .unwrap();
        let expected = std::iter::successors(Some(ratio), |&n| {
            n.is_multiple_of(prime).then_some(n / prime)
        })
        .count()
            - 1;
        assert!(
            lines[2 * i + 1].ends_with(&format!("from base={expected}")),
            "{}",
            lines[2 * i + 1]
        );
    }
    assert!(lines[6].contains("extension Sha order=100"));
    assert!(lines[7].ends_with("certified Sha order over F_q=4, v_ell(Sha)=2"));
}

#[test]
fn large_base_and_level_and_prime_to_ell_algebraic_factors_are_supported() {
    let output = stdout(interpret(&format!(
        r"{}
        rvKummerK3Iwasawa[C,85070591730234615865843651857942052864,2,340282366920938463463374607431768211455];
        rvKummerK3Iwasawa[C,1,3,1000];
        rvK3PencilShaIwasawa[P,2,2,1000]{{p,q}};
    ",
        base()
    )));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].contains(
        "positive Newton valuations=[129: 2], determinant valuations through breakpoint=[258]"
    ));
    assert!(lines[0]
        .ends_with("relative v_ell(Brauer) from base=680564733841876926926749214863536422910"));
    assert!(lines[1].contains("NS rank=18, geometric NS rank=20"));
    assert!(lines[1].contains("zero valuations=4, breakpoint k0=0, positive Newton valuations=[]"));
    assert!(lines[2].ends_with(
        "certified Sha order over F_q=4, NS rank over F_q=18, absolute tower valuation unresolved"
    ));
}

#[test]
fn actual_nonsplit_node_cycles_and_a_fully_algebraic_surface() {
    let output = stdout(interpret(
        r"
        rvunram F[3,1]{0,1}; let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1}; let two=rvUnramElem[F]{2};
        rvhyperelliptic C[F]{z,two,z,two,o,o};
        rvhyperelliptic N[F]{o,z,z,z,z,o};
        rvKummerK3Iwasawa[C,1,2,100];
        rvKummerK3Iwasawa[N,4,2,100];
    ",
    ));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].contains("NS rank=10, geometric NS rank=20, integral NS lattice fixed throughout tower, transcendental degree=2, zero valuations=12"), "{}",lines[0]);
    assert!(lines[1].contains("NS rank=22, geometric NS rank=22, integral NS lattice fixed throughout tower, transcendental degree=0, zero valuations=0"));
    assert!(lines
        .iter()
        .all(|line| line.ends_with("relative v_ell(Brauer) from base=0")));
}

#[test]
fn rank_jumps_bad_primes_and_wrong_sections_are_diagnosed() {
    for (operation, expected) in [
        ("rvKummerK3Iwasawa[C,1,2,100];", "rank jumps from 18 to 20 at level k=1"),
        ("rvK3PencilShaIwasawa[P,0,2,1];", "base degree must be positive"),
        ("rvKummerK3Iwasawa[C,2,4,100];", "Iwasawa prime must be prime"),
        ("rvKummerK3Iwasawa[C,2,9,100];", "9 is not prime"),
        ("rvKummerK3Iwasawa[C,2,4294967296,100];", "prime must fit in u32"),
        ("rvK3PencilShaIwasawa[P,3,5,1]{add(p,p),q};", "requires an exact base Sha order"),
        ("rvk3pencil Other[B,14,0]; rvK3PencilShaIwasawa[P,3,5,1]{p,zero[RVK3PencilSection(Other)]};", "expected type RVK3PencilSection(P)"),
    ] {
        let error = diagnostic(interpret(&format!("{}\n{operation}",base())));
        assert!(error.contains(expected), "{operation}: {error}");
    }
}
