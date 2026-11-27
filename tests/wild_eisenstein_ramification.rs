mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_ramification_polygons_and_wild_cyclotomic_filtrations() {
    let output = stdout(example("wild-eisenstein-ramification"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        &lines[..5],
        [
            "rvRamificationPolygon[K] = valuations=[3, 3, 0], vertices=[(1, 3), (3, 0)], embedding_distances=[(5, 2, 2)], different=5",
            "rvEisensteinVal = 5",
            "rvEisensteinTrace = precision=4, modulus=81, coordinates=[0]",
            "rvEisensteinNorm = precision=4, modulus=81, coordinates=[3]",
            "rvEisensteinPolynomial[C] = p=3, degree=6, coefficients=[3, 9, 18, 21, 15, 6, 1]",
        ]
    );
    assert!(lines[5].contains("vertices=[(1, 4), (3, 0), (6, 0)], embedding_distances=[(1, 1, 3), (3, 1, 2)], different=9"), "{output}");
    assert!(
        lines[6].contains("lower_groups=[[1, 2, 4, 5, 7, 8], [1, 4, 7], [1, 4, 7], [1]]"),
        "{output}"
    );
    assert!(
        lines[6].ends_with(
            "sizes=[6, 3, 3, 1], lower_breaks=[0, 2], upper_breaks=[(0, 1), (1, 1)], different=9"
        ),
        "{output}"
    );
    assert!(lines[7].ends_with("sizes=[4, 4, 2, 2, 1], lower_breaks=[1, 3], upper_breaks=[(1, 1), (2, 1)], different=8"), "{output}");
    assert!(lines[8].ends_with("sizes=[18, 9, 9, 3, 3, 3, 3, 3, 3, 1], lower_breaks=[0, 2, 8], upper_breaks=[(0, 1), (1, 1), (2, 1)], different=45"), "{output}");
}

#[test]
fn all_small_cyclotomic_automorphisms_compose_on_actual_elements() {
    let mut program = String::new();
    for (key, p, n, order, precision, coefficients) in [
        ("C", 3, 2, 9, 18, "3,9,18,21,15,6,1"),
        ("D", 2, 3, 8, 12, "2,4,6,4,1"),
    ] {
        writeln!(
            program,
            "rvunram B{key}[{p},3]{{0,1}}; rvcyclotomiclocal {key}[B{key},{n}];"
        )
        .unwrap();
        writeln!(program, "let p{key} = rvEisensteinPi[{key},{precision}]; let q{key} = rvEisensteinElem[{key},{precision}]([1],[2],[3]);").unwrap();
        let units: Vec<_> = (1..order).filter(|a| a % p != 0).collect();
        for &a in &units {
            writeln!(program, "assert rvEisensteinEval(rvEisensteinSigma[{a}](p{key}),rvTatePoly[{p},3]({coefficients})) == zero[RVEisenstein({key},{precision})];").unwrap();
            for &b in &units {
                writeln!(program, "assert rvEisensteinSigma[{a}](rvEisensteinSigma[{b}](q{key})) == rvEisensteinSigma[{}](q{key});", a*b%order).unwrap();
            }
        }
    }
    let output = stdout(interpret(&program));
    assert!(output.is_empty());
}

#[test]
fn general_eisenstein_relation_trace_and_norm_over_unramified_coefficients() {
    let output = stdout(interpret(
        "
        rvunram U[3,4]{1,0,1};
        rveisenstein K[U]{-3,-6,3,1};
        let pi = rvEisensteinPi[K,12];
        assert rvEisensteinMul(rvEisensteinMul(pi,pi),pi) == rvEisensteinElem[K,12]([3],[6],[-3]);
        let x = rvEisensteinElem[K,12]([1],[2],[1]);
        assert add(x,neg(x)) == zero[RVEisenstein(K,12)];
        rvEisensteinTrace(x); rvEisensteinNorm(x);
        let zpi = rvEisensteinElem[K,12]([],[0,1]);
        rvEisensteinTrace(zpi); rvEisensteinNorm(zpi);
        rvEisensteinTrace(rvEisensteinPi[K,2]);
        rvEisensteinNorm(rvEisensteinPi[K,5]);

        let mixed = rvEisensteinElem[K,5]([10,11],[10,11],[5,4]);
        show mixed;
        rvEisensteinVal(add(mixed,neg(mixed)));
        rveisenstein L[U]{-3,1};
        show rvEisensteinPi[L,4];
        rvRamificationPolygon[L];
        rvcyclotomiclocal C[U,2];
        let fixed = rvEisensteinElem[C,18]([1,1]);
        assert rvEisensteinSigma[4](fixed) == fixed;
        let cp = rvEisensteinPi[C,18];
        let zeta = add(cp,rvEisensteinElem[C,18]([1]));
        rvEisensteinTrace(zeta); rvEisensteinNorm(zeta);
        rvEisensteinVal(rvEisensteinEval(cp,rvTatePoly[3,3](9,36,63,60,30,6)));",
    ));
    // Independent resultants give the characteristic polynomial of x as
    // X^3-18X^2+81X-25, and Norm(z*pi) = -3z when z^2 = -1.
    assert_eq!(
        output,
        concat!(
            "rvEisensteinTrace = precision=4, modulus=81, coordinates=[18, 0]\n",
            "rvEisensteinNorm = precision=4, modulus=81, coordinates=[25, 0]\n",
            "rvEisensteinTrace = precision=4, modulus=81, coordinates=[0, 78]\n",
            "rvEisensteinNorm = precision=4, modulus=81, coordinates=[0, 78]\n",
            "rvEisensteinTrace = precision=0, unresolved\n",
            "rvEisensteinNorm = precision=1, modulus=3, coordinates=[0, 0]\n",
            "RVEisenstein(K,5){[[1, 2], [1, 2], [2, 1]]}\n",
            "rvEisensteinVal = at_least=5\n",
            "RVEisenstein(L,4){[[3, 0]]}\n",
            "rvRamificationPolygon[L] = valuations=[0], vertices=[(1, 0)], embedding_distances=[], different=0\n",
            "rvEisensteinTrace = precision=3, modulus=27, coordinates=[0, 0]\n",
            "rvEisensteinNorm = precision=3, modulus=27, coordinates=[1, 0]\n",
            "rvEisensteinVal = 9\n",
        )
    );
}

#[test]
fn finite_precision_distinguishes_unresolved_and_bad_inputs() {
    let header = "rvunram B[3,3]{0,1}; rveisenstein K[B]{-3,0,0,1}; rvcyclotomiclocal C[B,2];";
    assert_eq!(
        stdout(interpret(&format!("{header} rvHerbrand[C,3];"))),
        "rvHerbrand[C,3] = distances=[(2, 1), (5, 1), (8, 1)], valuation_at_least=3, unresolved=[4, 7]\n",
    );
    for (source, expected) in [
        ("rveisenstein Bad[B]{-9,0,0,1};", "v_p(a_0)=1"),
        (
            "rveisenstein Bad[B]{-3,1,0,1};",
            "p to divide every lower coefficient",
        ),
        ("rveisenstein Bad[B]{-3,0,0,-1};", "must be monic"),
        ("rvcyclotomiclocal D[B,0];", "level must be positive"),
        ("rvcyclotomiclocal D[B,4];", "degree must be in 1..=32"),
        ("rvcyclotomiclocal K[B,2];", "is already defined"),
        (
            "rvHerbrand[K,9];",
            "requires a named p-power cyclotomic extension",
        ),
        ("rvHerbrand[C,19];", "uniformizer precision"),
        (
            "show rvEisensteinSigma[3](rvEisensteinPi[C,9]);",
            "must be a unit",
        ),
        (
            "show rvEisensteinSigma[2](rvEisensteinPi[K,9]);",
            "require a named p-power cyclotomic extension",
        ),
        (
            "show rvEisensteinMul(rvEisensteinPi[K,9],rvEisensteinPi[C,9]);",
            "RVEisenstein(K,9)",
        ),
        (
            "show rvEisensteinEval(rvEisensteinPi[K,9],rvTatePoly[3,2](-3,0,0,1));",
            "precision is insufficient",
        ),
        (
            "show rvEisensteinEval(rvEisensteinPi[K,9],rvTatePoly[5,3](-3,0,0,1));",
            "primes do not match",
        ),
    ] {
        let message = diagnostic(interpret(&format!("{header} {source}")));
        assert!(message.contains(expected), "{message}");
    }
}
