mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn signed_leibniz_and_filtration_shift_have_nonzero_examples() {
    assert_eq!(
        stdout(example("bp-page-products")),
        "BPBock(3, 1, 0, 3, 12, 0){2*g0}\nBPBock(2, 2, 4, 1, 8, 5){g0}\n"
    );
}

#[test]
fn page_products_are_associative_and_next_is_multiplicative_on_cycles() {
    let program = r"
        let t = bpt[3,1,1,1,1];
        let x = bpBockClass[0,4,0](t);
        let y = bpBockClass[0,8,0](bpPow[2](t));
        assert bpCup(bpCup(y,x),x) == bpCup(y,bpCup(x,x));
        assert bpBockD(bpCup(y,x)) == bpCup(bpBockD(y),x);
        let unit = bpBockgen[3,1,0,0,0,0,0];
        assert bpCup(unit,x) == x;
        assert bpCup(x,unit) == x;
        assert bpBockNext(bpCup(x,x)) == bpCup(bpBockNext(x),bpBockNext(x));
        let v = bpBockClass[3,8,0](bpPow[4](bpv[2,2,0,1]));
        let p = bpBockgen[2,2,3,0,0,1,0];
        assert bpBockNext(bpCup(p,v)) == bpCup(bpBockNext(p),bpBockNext(v));
        let large = bpBockgen[3,1,1,0,0,1000,0];
        let u = bpBockgen[3,1,1,0,4,1000,0];
        assert bpCup(large,u) == bpBockgen[3,1,1,0,4,2000,0];
        assert bpBockD(bpCup(large,u)) == bpCup(large,bpBockD(u));
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn novikov_products_transport_the_base_in_the_graded_filtration() {
    let program = r"
        let v = bpv[3,1,0,1];
        let t = bpt[3,1,1,1,1];
        let x = bpNovClass[0,4,0](t);
        let y = bpNovClass[0,4,1](v);
        let p_t_squared = bpMul(bpint[3,1,1,3],bpPow[2](t));
        let difference = bpNovClass[0,8,1](p_t_squared);
        assert add(bpCup(x,y),pscale[2](bpCup(y,x))) == difference;
        assert difference == bpNovgen[3,1,0,1,8,1,1];
        assert bpNovD(bpCup(x,y)) == add(bpCup(bpNovD(x),y),pscale[2](bpCup(x,bpNovD(y))));
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn page_cup_checks_profile_page_and_grade_sums() {
    let maximum = usize::MAX;
    let cases = [
        (
            "bpCup(zero[BPBock(3,1,0,0,4,0)],zero[BPNov(3,1,0,0,4,0)])".to_owned(),
            "same BP page family",
        ),
        (
            "bpCup(zero[BPBock(3,1,0,0,4,0)],zero[BPBock(3,1,1,0,4,0)])".to_owned(),
            "same BP page family",
        ),
        (
            "bpCup(zero[BPBock(3,1,0,0,4,0)],zero[BPBock(3,2,0,0,4,0)])".to_owned(),
            "same prime, height",
        ),
        (
            format!("bpCup(zero[BPNov(3,1,1,0,{maximum},0)],zero[BPNov(3,1,1,0,1,0)])"),
            "internal degree is out of range",
        ),
        (
            format!("bpCup(zero[BPBock(3,1,1,0,0,{maximum})],zero[BPBock(3,1,1,0,0,1)])"),
            "filtration is out of range",
        ),
    ];
    for (term, reason) in cases {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
}
