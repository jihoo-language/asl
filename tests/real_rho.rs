mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn truncated_real_cobar_has_an_actual_nonzero_internal_bockstein() {
    assert_eq!(
        stdout(example("real-rho-quotient")),
        concat!(
            "RMRhoExt(2, 2, 4, 0){g0 + g1}\n",
            "rmrhoextdim[32, 0, 0, -32] = 1\n",
            "rmrhoextdim[33, 0, 0, -32] = 0\n",
            "rhoorder = 32\n",
            "RMExt(1, 32, 0){g0}\n",
            "RMRhoExt(32, 1, 32, 0){g0}\n",
        )
    );
}

#[test]
fn high_scalar_connecting_representatives_have_exact_orders_for_all_lower_exponents() {
    for e in 1..=32 {
        let source = format!(
            "let tau = rmcword[0,32]{{}};
             let cycle = extClass(rhoReduce[{e}](tau));
             let conn = rhoConnectRep(cycle);
             assert conn == rhoShift[{}](rmcword[0,0]{{([],[32])}});
             assert rhoShift[{e}](conn) == cobarD(tau);
             rhoorder(conn);",
            32 - e,
        );
        assert_eq!(stdout(interpret(&source)), format!("rhoorder = {e}\n"));
    }
}

#[test]
fn real_reduction_and_connecting_respect_cutoffs_and_the_les() {
    let source = r"
        let x = extClass(rmcword[0,0]{([],[2])});
        let r = rhoReduce[2](x);
        assert rhoConnect(r) == zero[RMExt(2,4,2)];
        assert realComplex(r) == realComplex(x);
        assert rhoReduce[1](rhoReduce[2](x)) == rhoReduce[1](x);
        assert rhoReduce[2](rhoReduce[2](x)) == r;
        assert rhoShift[2](r) == zero[RMRhoExt(2,1,0,-2)];

        let t = rmcword[0,2]{};
        let z = extClass(rhoReduce[2](t));
        assert rhoConnect(z) == x;
        assert rhoShift[2](rhoConnect(z)) == zero[RMExt(1,0,-2)];
        assert rhoConnect(rhoReduce[1](z)) == rhoShift[1](x);
        assert rhoBeta(z) == r;
        assert rhoBeta(r) == zero[RMRhoExt(2,2,4,2)];

        let h = rmcword[0,0]{([],[1])};
        assert rhoReduce[2](concat(h,t)) == concat(rhoReduce[2](h),rhoReduce[2](t));
        assert rhoReduce[2](concat(t,h)) == concat(rhoReduce[2](t),rhoReduce[2](h));
        assert concat(rhoReduce[2](rmcword[1,0]{}),rhoReduce[2](rmcword[1,0]{})) == zero[RMRhoCobar(2,0,-2,-2)];
        card(\q : RMRhoExt(2,0,0,-2) => S.hi);
        card(\q : RMRhoCobar(2,0,0,-2) => S.hi);
        rmrhoextdim[2,0,0,-2];
        rmbarRank[0,2,0,0];
        rmbarRank[1,2,2,0];
        rmbarRank[1,2,2,2];
    ";
    let result = stdout(interpret(source));
    assert_eq!(
        result,
        concat!(
            "card = 2 / 2\n",
            "card = 2 / 2\n",
            "rmrhoextdim[2, 0, 0, -2] = 1\n",
            "rmbarRank[0, 2, 0, 0] = 0\n",
            "rmbarRank[1, 2, 2, 0] = 2\n",
            "rmbarRank[1, 2, 2, 2] = 1\n",
        )
    );
}

#[test]
fn rho_truncation_checks_types_cutoffs_and_closure() {
    for (source, reason) in [
        (
            "show zero[RMRhoCobar(0,0,0,0)];",
            "truncation exponent must be positive",
        ),
        (
            "rmrhoextdim[0,0,0,0];",
            "truncation exponent must be positive",
        ),
        (
            "show rhoReduce[2](rhoReduce[1](rmcword[0,0]{}));",
            "cannot increase the truncation exponent",
        ),
        (
            "show rhoReduce[2](mcword[1]{});",
            "complex-motivic value when reducing to rho^1",
        ),
        (
            "show extClass(rhoReduce[33](rmcword[0,32]{}));",
            "requires a closed real-motivic cobar cochain",
        ),
        (
            "show rhoConnect(rhoReduce[1](rmcword[0,0]{}));",
            "an RMRhoExt class",
        ),
        (
            "show yoneda(rmrhoextgen[1,0,0,0,0],rmrhoextgen[2,0,0,0,0]);",
            "same truncation exponent",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
