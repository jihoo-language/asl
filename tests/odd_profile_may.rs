mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn restriction_kills_the_full_fifth_differential_and_the_lift_survives() {
    assert_eq!(
        stdout(example("odd-profile-may")),
        "PMayAn(3, 1, 6, 3, 15, 6){g0}\npmayandim[3, 1, 6, 3, 15, 6] = 1\n"
    );
}

#[test]
fn odd_finite_may_products_restriction_and_passage_respect_prime_coefficients() {
    let source = "
        let a = mayClass[1, 0](pword[5]{([], [0])});
        let x = mayClass[1, 7](pword[5]{([1], [])});
        let v = mayClass[0, 6](pword[5]{([], [1])});
        let fa = resA[1](a);
        let fx = resA[1](x);
        let fv = resA[1](v);
        assert mayD(fv) == zero[PMayAn(5, 1, 0, 2, 9, 6)];
        let page1 = mayNext(fv);
        assert mayPrev(pscale[3](page1)) == pscale[3](fv);
        assert mayD(page1) == mayMul(fx, fa);
        assert mayNext(mayD(page1)) == zero[PMayAn(5, 1, 2, 2, 9, 7)];
        assert resA[1](mayMul(x, a)) == mayMul(fx, fa);
        assert mayMul(pscale[2](fa), pscale[3](fa)) == mayMul(fa, fa);
        assert mayMul(fx, fx) == zero[PMayAn(5, 1, 1, 2, 16, 14)];
        assert resA[0](resA[1](a)) == resA[0](a);
        assert mayNext(resA[0](a)) == resA[0](mayNext(a));
        betti PMayAn(5, 0, 0, 2, 2, 0);
        pmayandim[5, 0, 10, 2, 2, 1];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "betti(PMayAn(5, 0, 0, 2, 2, 0)) = [5]\npmayandim[5, 0, 10, 2, 2, 1] = 0\n"
    );
}

#[test]
fn odd_finite_pages_check_profile_prime_and_cycle_conditions() {
    for (expression, expected) in [
        (
            "resA[2](zero[PMayAn(3, 1, 0, 0, 0, 0)])",
            "smaller Steenrod algebra",
        ),
        (
            "mayMul(zero[PMay(3, 0, 0, 0, 0)], zero[PMayAn(3, 1, 0, 0, 0, 0)])",
            "same Steenrod algebra profile",
        ),
        (
            "mayMul(zero[PMayAn(3, 1, 0, 0, 0, 0)], zero[PMayAn(5, 1, 0, 0, 0, 0)])",
            "same prime",
        ),
        ("mayNext(pmayangen[3, 1, 1, 1, 5, 2, 0])", "cycle for d_1"),
        (
            "mayClass[2, 2](resA[1](pword[3]{([], [1])}))",
            "cycle condition",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(expected), "{expression}\n{error}");
    }
}
