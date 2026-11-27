mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_comparison_recovers_non_dyadic_bocksteins_and_sparse_high_powers() {
    assert_eq!(
        stdout(example("rho-minimal")),
        concat!(
            "RhoMin(2, 2, -3){rho^6*g0 + rho^7*g0}\n",
            "RhoBock(6, 2, 2, -3, 6){g0}\n",
            "RhoBock(7, 3, 4, -2, 7){g0}\n",
            "RhoBock(7, 4, 3, -4, 10){g0}\n",
            "rhoorder = 32\n",
            "RhoMin(1, 0, -32){rho^32*g0}\n",
            "MExt(1, 32, 0){g0}\n",
        ),
    );
}

#[test]
fn scalar_differentials_and_rho_linearity() {
    for n in [1, 2, 4, 8] {
        let source = format!(
            "let a = rhoMin[0](mextgen[0,0,-{n},0]);
             let x = rhoMinD(a);
             let c = extClass(realComplex(rmcword[0,0]{{([],[{n}])}}));
             assert x == rhoMin[{n}](c);
             assert rhoMinI(a) == rmcword[0,{n}]{{}};
             assert rhoMinI(x) == cobarD(rhoMinI(a));
             assert rhoMinD(rhoShift[3](a)) == rhoShift[3](x);
             assert rhoMinI(rhoShift[3](x)) == rhoShift[3](rhoMinI(x));
             assert rhoMinP(rhoShift[3](rhoMinI(x))) == rhoShift[3](x);
             assert rhoMinCoef[0](x) == zero[MExt(1,0,-{n})];
             assert rhoMinCoef[{n}](x) == c;
             assert rhoMinD(x) == zero[RhoMin(2,0,-{n})];"
        );
        assert!(stdout(interpret(&source)).is_empty());
    }
}

#[test]
fn polynomial_operations_and_corrected_contraction_on_noncycles() {
    let source = r"
        let u = rmcword[0,0]{([0,1],[]),([],[1])};
        let du = cobarD(u);
        let px = rhoMinP(u);
        let hi = rhoMinH(u);
        assert add(cobarD(hi),rhoMinH(du)) == add(u,rhoMinI(px));
        assert rhoMinD(px) == rhoMinP(du);
        assert rhoMinP(rhoMinI(px)) == px;
        assert rhoMinH(rhoMinI(px)) == add(hi,hi);
        assert rhoMinH(rhoShift[2](u)) == rhoShift[2](hi);
        assert rhoMinP(rhoShift[2](u)) == rhoShift[2](px);
        assert rhoMinP(hi) == add(rhoMinP(hi),rhoMinP(hi));
        assert rhoMinH(hi) == add(rhoMinH(hi),rhoMinH(hi));
        let scalar = \x : RhoMin(0,-2,-3) => rhoMinI(x);
        card(\x : RhoMin(0,-2,-3) => S.hi);
        assert scalar(rhoMin[2](mextgen[0,0,-1,0])) == rmcword[2,1]{};
        card(\x : RhoMin(2,0,1) => S.hi);
    ";
    assert_eq!(stdout(interpret(source)), "card = 2 / 2\ncard = 1 / 1\n");
}

#[test]
fn rho_minimal_rejects_invalid_grade_exponents_and_operands() {
    for (source, reason) in [
        ("show rhoMinD(rmcword[0,0]{});", "a RhoMin polynomial"),
        ("show rhoMin[0](rmextgen[0,0,0,0]);", "an MExt class"),
        ("show rhoMinH(rmcword[0,0]{});", "positive cobar length"),
        (
            "show rhoShift[-1](zero[RhoMin(0,0,0)]);",
            "unsigned exponent",
        ),
        (
            "show rhoMinCoef[1](zero[RhoMin(0,0,9223372036854775807)]);",
            "coefficient grade exceeds",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
