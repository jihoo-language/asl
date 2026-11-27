mod common;

use common::{diagnostic, example, interpret, stdout};

const F23: &str = include_str!("../examples/hyperelliptic-k3-iwasawa.asl");
const DEMO: &str = include_str!("../examples/hyperelliptic-k3-crystalline-iwasawa.asl");

fn f23() -> &'static str {
    F23.split_once("rvKummerK3Iwasawa[C,2,2,1000000]")
        .unwrap()
        .0
}

#[test]
fn supersingular_mu_and_absolute_orders_match_independent_artin_tate_calculation() {
    let output = stdout(example("hyperelliptic-k3-crystalline-iwasawa"));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].contains("Artin invariant=1, geometric NS discriminant=-25, Brauer order=1"));
    assert!(
        lines[1].contains("Artin invariant=1, geometric NS discriminant=-25, Brauer order=390625")
    );
    for (level, value) in [
        (0, "0"),
        (1, "8"),
        (2, "48"),
        (1_000_000, "2*5^1000000 - 2"),
    ] {
        let prefix = format!("rvKummerK3Iwasawa[S, 2, 5, {level}]");
        let row = lines.iter().find(|line| line.starts_with(&prefix)).unwrap();
        assert!(row.contains("NS rank=22, geometric NS rank=22"));
        assert!(row.contains("formal Brauer height=infinite, nonzero v_p(u) at base=[], negative Newton valuations=[], positive determinant valuations through breakpoint=[0], mu=2, lambda=0"), "{row}");
        assert!(row.ends_with(&format!("Artin invariant=1, NS discriminant=-25, base v_ell(Brauer)=0, v_ell(Brauer)={value}")), "{row}");
    }
    let direct = DEMO.split_once("rvKummerK3ArtinTate[S,2]").unwrap().0;
    let output = stdout(interpret(&format!("{direct} rvKummerK3ArtinTate[S,50];")));
    let order: u128 = output
        .trim()
        .rsplit_once(", Brauer order=")
        .unwrap()
        .1
        .parse()
        .unwrap();
    assert_eq!(order, 5u128.pow(48));
}

#[test]
fn ordinary_and_actual_height_two_slopes_force_mu_zero() {
    let output = stdout(example("hyperelliptic-k3-crystalline-iwasawa"));
    let ordinary = output
        .lines()
        .find(|line| line.starts_with("rvKummerK3Iwasawa[C"))
        .unwrap();
    assert!(ordinary.contains("zero valuations=1, breakpoint k0=0, positive Newton valuations=[]"));
    assert!(ordinary.contains("formal Brauer height=1, nonzero v_p(u) at base=[-2: 1, 2: 1], negative Newton valuations=[-2: 1], positive determinant valuations through breakpoint=[0], mu=0, lambda=0"),"{ordinary}");
    assert!(ordinary.ends_with("relative v_ell(Brauer) from base=0"));
    let height_two = output
        .lines()
        .find(|line| line.starts_with("rvKummerK3Iwasawa[D"))
        .unwrap();
    assert!(output.contains("rvKummerK3Height[D] = 2"));
    assert!(height_two.contains("formal Brauer height=2, nonzero v_p(u) at base=[-1/2: 2, 1/2: 2], negative Newton valuations=[-1/2: 2], positive determinant valuations through breakpoint=[0], mu=0, lambda=0"), "{height_two}");
    assert!(height_two.ends_with("relative v_ell(Brauer) from base=0"));
}

#[test]
fn inherited_sha_certificate_and_rank_jump_scope_are_correct_in_characteristic() {
    let output = stdout(interpret(&format!(
        r"{}
        rvK3PencilShaIwasawa[P,1,23,1000]{{p,q}};
        rvK3PencilShaIwasawa[P,2,23,1000]{{p,q}};
        rvK3PencilShaIwasawa[P,2,23,1000];
        ",
        f23()
    )));
    let lines: Vec<_> = output.lines().collect();
    assert!(lines[0].contains("NS rank=18, geometric NS rank=20"));
    assert!(lines[0].ends_with("certified Sha order over F_q=4, v_ell(Sha)=0"));
    assert!(lines[1].ends_with(
        "certified Sha order over F_q=4, NS rank over F_q=18, absolute tower valuation unresolved"
    ));
    assert!(lines[2].ends_with("relative v_ell(Sha) from base=0"));
}

#[test]
fn precision_guards_and_future_p_rank_jumps_do_not_produce_a_certificate() {
    let failure = diagnostic(interpret(&format!(
        "{} rvKummerK3Iwasawa[C,1000,23,100];",
        f23()
    )));
    assert!(
        failure.contains("checked ell-adic precision after clearing Frobenius denominators"),
        "{failure}"
    );
    let failure = diagnostic(interpret(
        r"
        rvunram F[3,1]{0,1}; let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1}; let two=rvUnramElem[F]{2};
        rvhyperelliptic C[F]{z,two,z,two,o,o};
        rvKummerK3Iwasawa[C,1,3,100];
    ",
    ));
    assert!(
        failure.contains("rank jumps from 10 to 20 at level k=1"),
        "{failure}"
    );
}
