mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn p_adic_differentials_and_corrected_transitions_use_integral_lifts() {
    assert!(stdout(example("bp-pages")).contains("bpBockdim[2, 2, 5, 0, 8, 0] = 0"));
    let program = r"
        let x = bpBockClass[3,8,0](bpPow[4](bpv[2,2,0,1]));
        assert bpBockD(x) == zero[BPBock(2,2,3,1,8,3)];
        let y = bpBockNext(x);
        assert bpBockNext(bpBockPrev(y)) == y;
        assert bpBockClass[4,8,0](bpBockRep(y)) == y;
        assert y == bpBockgen[2,2,4,0,8,0,0];
        let d4 = bpBockD(y);
        assert bpBockD(d4) == zero[BPBock(2,2,4,2,8,8)];
        assert bpBockNext(d4) == zero[BPBock(2,2,5,1,8,4)];
        assert add(d4,pscale[3](d4)) == zero[BPBock(2,2,4,1,8,4)];
        let height_one = bpBockClass[3,8,0](bpPow[4](bpv[2,1,0,1]));
        assert bpBockD(height_one) == bpBockgen[2,1,3,1,8,3,0];
        let odd = bpBockClass[2,12,0](bpPow[3](bpv[3,2,0,1]));
        let odd_d2 = bpBockD(odd);
        assert odd_d2 == bpBockgen[3,2,2,1,12,2,0];
        assert bpBockD(pscale[2](odd)) == pscale[2](odd_d2);
        assert bpBockD(bpBockgen[3,1,1,0,4,1000,0]) == bpBockgen[3,1,1,1,4,1001,0];
        betti BPBock(3,2,2,0,12,0);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(BPBock(3, 2, 2, 0, 12, 0)) = [3]\n"
    );
}

#[test]
fn the_novikov_page_two_transgression_has_a_nonzero_integral_representative() {
    let program = r"
        let v1 = bpv[2,2,1,1];
        let v2 = bpv[2,2,1,2];
        let t1 = bpt[2,2,1,1,1];
        let t2 = bpt[2,2,1,1,2];
        let term1 = bpMul(bpPow[2](v1),bpPow[6](t1));
        let term2 = bpNeg(bpMul(bpPow[2](v1),bpPow[2](t2)));
        let term3 = bpMul(bpPow[2](v2),bpPow[2](t1));
        let term4 = bpMul(bpint[2,2,1,-2],bpPow[8](t1));
        let term5 = bpMul(bpint[2,2,1,-4],bpMul(bpPow[2](t1),bpPow[2](t2)));
        let representative = add(add(term1,term2),add(term3,add(term4,term5)));
        let source = bpNovClass[2,16,1](representative);
        assert source == bpNovgen[2,2,2,1,16,1,0];
        let image = bpNovD(source);
        assert image == bpNovgen[2,2,2,2,16,3,0];
        assert image == bpNovClass[2,16,3](bpD(representative));
        assert bpNovD(image) == zero[BPNov(2,2,2,3,16,5)];
        assert bpNovNext(bpNovPrev(source)) == source;
        assert bpNovClass[2,16,3](bpNovRep(image)) == image;
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn page_zero_and_invalid_classes_have_explicit_witness_conditions() {
    let program = r"
        let v = bpv[2,1,0,1];
        let x = bpBockClass[0,2,0](v);
        assert x == bpBockgen[2,1,0,0,2,0,0];
        assert bpBockD(x) == zero[BPBock(2,1,0,1,2,0)];
        assert bpBockNext(x) == bpBockgen[2,1,1,0,2,0,0];
        assert bpNovClass[0,2,0](v) == zero[BPNov(2,1,0,0,2,0)];
        assert bpNovClass[0,2,1](v) == bpNovgen[2,1,0,0,2,1,0];
    ";
    assert_eq!(stdout(interpret(program)), "");
    for (term, reason) in [
        (
            "bpBockClass[2,2,0](bpv[2,1,0,1])",
            "requires a filtered cycle",
        ),
        (
            "bpBockClass[0,2,1](bpv[2,1,0,1])",
            "requires a filtered cycle",
        ),
        (
            "bpBockClass[0,4,0](bpv[2,1,0,1])",
            "outside its normalized cobar basis",
        ),
        (
            "bpBockClass[0,2,0](bpv[2,1,1,1])",
            "outside its normalized cobar basis",
        ),
        (
            "bpBockClass[0,2,0](bpReduce[1](bpv[2,1,0,1]))",
            "an integral BP cochain",
        ),
        (
            "bpBockNext(bpBockgen[2,1,1,0,2,0,0])",
            "requires a cycle for the page differential",
        ),
        (
            "bpNovPrev(zero[BPNov(2,1,0,0,2,0)])",
            "page zero has no previous page",
        ),
        ("bpNovD(zero[BPBock(2,1,0,0,2,0)])", "BPNov"),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
    let filtration = usize::MAX;
    let program = format!("show bpBockD(zero[BPBock(2,1,1,0,2,{filtration})]);");
    assert!(diagnostic(interpret(&program)).contains("differential filtration is out of range"));
}
