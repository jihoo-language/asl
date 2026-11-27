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

#[test]
fn high_internal_bockstein_square_has_an_exact_witness_without_ambient_targets() {
    for e in 17..=32 {
        let primitive_check = if matches!(e, 17 | 31 | 32) {
            "assert cobarD(extPrimitive(raw))==raw;"
        } else {
            ""
        };
        let source = format!(
            "let a=extClass(rhoReduce[{e}](rmcword[0,32]{{}}));
             let b=rhoBeta(a);
             assert b==rmrhoextgen[{e},1,{e},{},0];
             assert rhoBeta(b)==zero[RMRhoExt({e},2,{},{})];
             let raw=rhoBetaRep(b);
             let u=rhoBetaNull(b);
             assert raw==rhoReduce[{e}](rhoConnectRep(b));
             assert cobarD(u)==raw;
             {primitive_check}
             show raw;
             show u;",
            e - 32,
            2 * e,
            2 * e - 32,
        );
        let output = stdout(interpret(&source));
        let mut lines = output.lines();
        let (raw, witness) = (lines.next().unwrap(), lines.next().unwrap());
        if e != 32 {
            for (s, text) in [(2, raw), (1, witness)] {
                assert!(
                    text.starts_with(&format!(
                        "RMRhoCobar({e}, {s}, {}, {}){{",
                        2 * e,
                        2 * e - 32
                    )) && !text.ends_with("{0}"),
                    "{text}",
                );
            }
        }
        assert!(lines.next().is_none());
    }
}

#[test]
fn internal_bockstein_witnesses_preserve_naturality_and_diagnose_nonboundaries() {
    let source = r"
        let z=extClass(rhoReduce[2](add(rmcword[0,3]{([1],[])},
            add(rmcword[1,2]{([1],[1])},rmcword[1,2]{([],[0,1])}))));
        let b=rhoBeta(z);
        assert extClass(rhoBetaRep(z))==b;
        assert rhoBeta(b)==zero[RMRhoExt(2,3,6,2)];
        let raw=rhoBetaRep(b);
        assert cobarD(rhoBetaNull(b))==raw;
        show raw;

        let a=extClass(rhoReduce[4](rmcword[0,4]{}));
        let a3=rhoReduce[3](a);
        assert rhoBeta(a3)==rhoShift[1](rhoReduce[3](rhoBeta(a)));
        assert cobarD(rhoBetaNull(rhoBeta(a3)))==rhoBetaRep(rhoBeta(a3));
        assert rhoBetaNull(zero[RMRhoExt(2,0,0,-2)])==zero[RMRhoCobar(2,0,2,0)];
    ";
    let output = stdout(interpret(source));
    assert!(output.starts_with("RMRhoCobar(2, 3, 6, 2){"), "{output}");
    assert!(!output.ends_with("{0}\n"), "{output}");

    let input = "let a=extClass(rhoReduce[1](rmcword[0,1]{}));";
    for expression in ["rhoBetaNull(a)", "extPrimitive(rhoBetaRep(a))"] {
        let error = diagnostic(interpret(&format!("{input} show {expression};")));
        assert!(
            error.contains("requires a rho-truncated cobar boundary"),
            "{error}"
        );
    }
    let error = diagnostic(interpret("show rhoBetaRep(rhoReduce[1](rmcword[0,1]{}));"));
    assert!(error.contains("an RMRhoExt class"), "{error}");
}

#[test]
fn internal_bockstein_null_example_returns_a_real_cochain_identity() {
    assert_eq!(
        stdout(example("real-rho-bockstein-null")),
        concat!(
            "RMRhoExt(17, 1, 17, -15){g0}\n",
            "RMRhoCobar(17, 2, 34, 2){rmcword[14, 8]{([8], []), ([16], [])}}\n",
            "RMRhoCobar(17, 1, 34, 2){rmcword[14, 8]{([24], [])} + rmcword[14, 8]{([0, 8], [])}}\n",
        ),
    );
}
