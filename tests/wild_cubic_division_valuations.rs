mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

const BASE: &str = "
    rvunram B[3,7]{0,1};
    rvcyclotomiclocal K[B,1];
    let lambda = rvEisensteinPi[K,14];
    let lambda2 = rvEisensteinMul(lambda,lambda);
    let lambda3 = rvEisensteinMul(lambda,lambda2);
    let one = rvEisensteinElem[K,14]([1]);
    let z = add(one,lambda);
    let u2 = add(one,lambda2);
    let u3 = add(one,lambda3);
";

#[test]
fn both_invariants_have_actual_graded_commutators_and_expected_cardinalities() {
    let output = stdout(example("wild-cubic-division-valuations"));
    for expected in [
        "invariant=1/3, frobenius_power=1, ramification_index=3, residue_degree=3, residue_size=27, precision=4, coefficient_precisions=[2, 1, 1], order_size=531441, units=511758, Pi=v",
        "rvCubicDivCommutatorCheck = layers=[1, 1], grade=2, leading=[1, 0, 0], valuation=2",
        "rvCubicDivVal = valuation=3, leading=[1, 0, 0], norm_precision=4, reduced_norm=[0, 6]",
        "invariant=2/3, frobenius_power=2, ramification_index=3, residue_degree=3, residue_size=27, precision=7, coefficient_precisions=[3, 2, 2], order_size=10460353203, units=10072932714, Pi=v^2/lambda",
        "rvCubicDivCommutatorCheck = layers=[1, 1], grade=2, leading=[2, 0, 0], valuation=2",
        "rvCubicDivVal = valuation=1, leading=[1, 0, 0], norm_precision=3, reduced_norm=[0, 1]",
        "RVEisenstein(K,3){[[0], [1]]}",
    ] {
        assert!(output.contains(expected), "{expected}\n{output}");
    }
}

#[test]
fn canonical_symbol_homomorphisms_and_reduced_norms_are_compatible() {
    let mut program = BASE.to_owned();
    for (t, v_param) in [(1, "lambda"), (2, "lambda2")] {
        writeln!(
            program,
            "rvcubicdiv D{t}[K,{t},8]; rvcubicalgebra C{t}[K,8](u3,{v_param});"
        )
        .unwrap();
        writeln!(program, "
            let e{t} = rvCubicDivElem[D{t}]([1]);
            let l{t} = rvCubicDivElem[D{t}]([0,1]);
            let p{t} = rvCubicDivPi[D{t}];
            let t{t} = rvCubicDivTheta[D{t}];
            let pp{t} = rvCubicDivMul(p{t},p{t});
            let th{t} = rvCubicDivMul(t{t},t{t});
            assert rvCubicDivMul(pp{t},p{t}) == l{t};
            assert rvCubicDivMul(th{t},t{t}) == add(
                add(e{t},neg(rvCubicDivMul(rvCubicDivElem[D{t}]([2,1]),t{t}))),
                rvCubicDivMul(rvCubicDivElem[D{t}]([3,1]),th{t}));
            assert rvCubicDivConjugatePi(rvCubicDivConjugatePi(rvCubicDivConjugatePi(t{t}))) == t{t};
            assert rvCubicDivMul(p{t},t{t}) == rvCubicDivMul(rvCubicDivConjugatePi(t{t}),p{t});
            assert rvCubicDivIsUnit(add(rvCubicDivConjugatePi(t{t}),neg(t{t}))) == S.hi;
            assert rvCubicDivNorm[3](t{t}) == rvEisensteinElem[K,3]([1]);
            assert rvCubicDivNorm[4](p{t}) == rvEisensteinPi[K,4];
            assert rvCubicDivNorm[3](add(e{t},p{t})) == rvEisensteinElem[K,3]([1],[1]);
            let su{t} = rvCubicAlgBasis[C{t},1,0];
            let v{t} = rvCubicAlgBasis[C{t},0,1];
            let x{t} = add(rvCubicAlgBasis[C{t},0,0],add(su{t},rvCubicAlgBasis[C{t},2,2]));
            let y{t} = add(su{t},v{t});
            let ud{t} = rvCubicDivFromSymbol[D{t}](su{t});
            let vd{t} = rvCubicDivFromSymbol[D{t}](v{t});
            let xd{t} = rvCubicDivFromSymbol[D{t}](x{t});
            let yd{t} = rvCubicDivFromSymbol[D{t}](y{t});
            assert ud{t} == add(e{t},rvCubicDivMul(l{t},t{t}));
            assert vd{t} == {};
            assert rvCubicDivMul(vd{t},ud{t}) ==
              rvCubicDivMul(rvCubicDivElem[D{t}]([1,1]),rvCubicDivMul(ud{t},vd{t}));
            assert rvCubicDivFromSymbol[D{t}](rvCubicAlgMul(x{t},y{t})) == rvCubicDivMul(xd{t},yd{t});
            assert rvCubicDivNorm[3](rvCubicDivMul(xd{t},yd{t})) ==
              rvEisensteinMul(rvCubicDivNorm[3](xd{t}),rvCubicDivNorm[3](yd{t}));
        ", if t == 1 { format!("p{t}") } else { format!("pp{t}") }).unwrap();
        if t == 2 {
            writeln!(
                program,
                "assert rvCubicDivMul(vd{t},vd{t}) == rvCubicDivMul(l{t},p{t});"
            )
            .unwrap();
        }
    }
    writeln!(
        program,
        "show rvCubicReducedNorm(x1); show rvCubicDivNorm[3](xd1);
         show rvCubicReducedNorm(x2); show rvCubicDivNorm[3](xd2);"
    )
    .unwrap();
    let output = stdout(interpret(&program));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 4, "{output}");
    for pair in lines.chunks_exact(2) {
        let expected: Vec<u64> = pair[0]
            .split_once("{[[")
            .unwrap()
            .1
            .trim_end_matches("]]}")
            .split("], [")
            .map(|digit| digit.parse().unwrap())
            .collect();
        assert_eq!(
            pair[1],
            format!(
                "RVEisenstein(K,3){{[[{}], [{}]]}}",
                expected[0] % 9,
                expected[1] % 3
            ),
        );
    }
}

#[test]
fn residue_field_units_truncation_and_cancellation() {
    let mut program = format!(
        "{BASE} rvcubicdiv D[K,2,8]; rvcubicdiv F[K,2,4]; rvcubicdiv R[K,2,1];
         let e = rvCubicDivElem[R]([1]); let p = rvCubicDivPi[D]; let t = rvCubicDivTheta[D];
         let x = add(rvCubicDivElem[D]([1]),p);
         assert rvCubicDivPi[R] == zero[RVCubicDiv(R)];
         assert rvCubicDivReduce[R](p) == zero[RVCubicDiv(R)];
         assert rvCubicDivReduce[F](rvCubicDivInverse(x)) == rvCubicDivInverse(rvCubicDivReduce[F](x));
         assert rvCubicDivReduce[F](rvCubicDivMul(t,x)) ==
           rvCubicDivMul(rvCubicDivReduce[F](t),rvCubicDivReduce[F](x));
         rvCubicDivVal(add(x,neg(x)));
         rvCubicDivVal(rvCubicDivElem[D]([0,0],[],[],[0,0],[],[],[],[],[3,0]));
         rvCubicDivCommutatorCheck(x,x);
         let high = add(rvCubicDivElem[F]([1]),rvCubicDivElem[F]([0,1]));
         rvCubicDivCommutatorCheck(high,high);
         rvcubicdiv Wide[K,2,42];
         rvCubicDivVal(rvCubicDivElem[Wide]([0],[],[],[],[],[],[0,729]));
         rvCubicDivVal(zero[RVCubicDiv(Wide)]);
         rvunram Bsmall[3,1]{{0,1}};
         rvcyclotomiclocal Ksmall[Bsmall,1];
         rvcubicdiv Small[Ksmall,1,6];
         rvCubicDivVal(rvCubicDivPi[Small]);"
    );
    for a in 0..3 {
        for b in 0..3 {
            for c in 0..3 {
                let digit = format!("r{a}{b}{c}");
                writeln!(program, "let {digit}=rvCubicDivElem[R]([{a}],[{b}],[{c}]);").unwrap();
                if (a, b, c) == (0, 0, 0) {
                    writeln!(program, "assert rvCubicDivIsUnit({digit}) == S.lo;").unwrap();
                } else {
                    writeln!(
                        program,
                        "assert rvCubicDivIsUnit({digit}) == S.hi;
                        assert rvCubicDivMul(rvCubicDivInverse({digit}),{digit}) == e;"
                    )
                    .unwrap();
                }
            }
        }
    }
    let output = stdout(interpret(&program));
    assert_eq!(
        output
            .matches("at_least=8, leading=unresolved, norm_precision=8, reduced_norm=[0, 0]")
            .count(),
        2,
        "{output}"
    );
    for expected in [
        "layers=[1, 1], grade=2, leading=[0, 0, 0], at_least=8",
        "layers=[3, 3], at_least=4, grade=unresolved",
        "valuation=41, leading=[1, 0, 0], norm_precision=14, reduced_norm=[0, 0]",
        "at_least=42, leading=unresolved, norm_precision=14, reduced_norm=[0, 0]",
        "valuation=1, leading=[1, 0, 0], norm_precision=2, reduced_norm=[0, 1]",
    ] {
        assert!(output.contains(expected), "{expected}\n{output}");
    }
}

#[test]
fn uncertified_norms_nominality_and_inappropriate_symbol_charts_fail() {
    let fixtures = format!(
        "{BASE} rvcubicdiv D[K,2,4]; rvcubicdiv E[K,1,4]; rvcubicdiv Low[K,2,1];
         rvcubicalgebra C[K,8](u3,lambda); rvcubicalgebra Bad[K,8](z,u2);
         rvcubicalgebra Tiny[K,2](u3,lambda2);
         let p = rvCubicDivPi[D]; let e = rvCubicDivElem[D]([1]);"
    );
    for (statement, expected) in [
        ("rvcubicdiv N[K,0,4];", "numerator must be 1 or 2"),
        ("rvcubicdiv D[K,2,4];", "already defined"),
        (
            "rvcubicdiv N[K,2,43];",
            "precision must be positive and supported",
        ),
        (
            "show rvCubicDivInverse(p);",
            "maximal order requires a unit",
        ),
        ("show rvCubicDivNorm[3](e);", "lambda precision at most 2"),
        ("show rvCubicDivMul(p,rvCubicDivPi[E]);", "RVCubicDiv(D)"),
        (
            "show rvCubicDivReduce[D](rvCubicDivPi[E]);",
            "same named base and invariant",
        ),
        (
            "show rvCubicDivReduce[D](rvCubicDivPi[Low]);",
            "sufficient source precision",
        ),
        (
            "show rvCubicDivFromSymbol[D](rvCubicAlgBasis[C,1,0]);",
            "canonical division invariant",
        ),
        (
            "show rvCubicDivFromSymbol[D](rvCubicAlgBasis[Bad,1,0]);",
            "parameters U_3 and lambda^t",
        ),
        (
            "show rvCubicDivFromSymbol[D](rvCubicAlgBasis[Tiny,1,0]);",
            "sufficient source precision",
        ),
        (
            "rvCubicDivCommutatorCheck(e,neg(e));",
            "requires principal units",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{fixtures} {statement}")));
        assert!(error.contains(expected), "{statement}\n{error}");
    }
}
