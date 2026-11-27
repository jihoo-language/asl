mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn wild_uniformizers_breaks_and_absolute_polynomials_are_computed() {
    let output = stdout(example("cyclotomic-kummer-wild-towers"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        &lines[..7],
        [
            "rvKummerLocalData[First,24] = b=1, w_order=1, uniformizer=(w^1)/lambda^0, distances=[3, 3], lower_break=2, upper_break=2, relative_different=6, absolute_different=9",
            "rvKummerLocalPolynomial[First,24] = degree=6, precision=4, modulus=81, coefficients=[3, 9, 18, 21, 15, 6, 1], base_trace=[75], base_norm=[3], derivative_valuation=9",
            "rvKummerLocalData[Second,24] = b=2, w_order=2, uniformizer=(w^2)/lambda^1, distances=[2, 2], lower_break=1, upper_break=1, relative_different=4, absolute_different=7",
            "rvKummerLocalPolynomial[Second,24] = degree=6, precision=4, modulus=81, coefficients=[3, 18, 66, 30, 18, 3, 1], base_trace=[78], base_norm=[3], derivative_valuation=7",
            "rvKummerLocalTrace = base=K, lambda_precision=8, coordinates=[[78], [80]]",
            "rvKummerLocalNorm = base=K, lambda_precision=8, coordinates=[[0], [1]]",
            "rvKummerLocalNorm = base=K, lambda_precision=8, coordinates=[[79], [78]]",
        ]
    );
    assert!(lines[7].contains("trace=1, kind=unramified, residue_degree=3, different=0, residue_polynomial=[[2], [2], [0], [1]], residual_roots=[]"), "{output}");
}

#[test]
fn galois_actions_use_the_actual_kummer_generator_and_nominal_base() {
    let output = stdout(interpret(
        "
        rvunram F[5,6]{0,1};
        rvcyclotomiclocal K[F,1];
        rvkummerlocal L[K,4]([2],[1]);
        let l0 = rvEisensteinPi[K,3];
        let o0 = rvEisensteinElem[K,3]([1]);
        let l = rvKummerLocalBase[L,12](l0);
        let one = rvKummerLocalBase[L,12](o0);
        let zeta = add(one,l);
        let y = rvKummerLocalY[L,12];
        let pi = rvKummerLocalPi[L,12];
        let a = add(y,pi);
        assert add(a,neg(a)) == zero[RVKummerLocal(L,12)];
        assert rvKummerLocalSigma[1](y) == rvKummerLocalMul(zeta,y);
        assert rvKummerLocalSigma[2](rvKummerLocalSigma[3](a)) == a;
        assert rvKummerLocalSigma[3](rvKummerLocalSigma[-1](a)) == rvKummerLocalSigma[2](a);
        assert rvKummerLocalSigma[1](rvKummerLocalMul(y,pi)) ==
            rvKummerLocalMul(rvKummerLocalSigma[1](y),rvKummerLocalSigma[1](pi));
        assert rvKummerLocalSigma[4](l) == l;
        rvKummerLocalData[L,12];
        rvKummerLocalVal(pi);
        rvKummerLocalTrace(y);",
    ));
    assert_eq!(
        output,
        concat!(
            "rvKummerLocalData[L,12] = b=4, w_order=4, uniformizer=(w^4)/lambda^3, distances=[2, 2, 2, 2], lower_break=1, upper_break=1, relative_different=8, absolute_different=23\n",
            "rvKummerLocalVal = 1\n",
            "rvKummerLocalTrace = base=K, lambda_precision=3, coordinates=[[0], [0], [0], [0]]\n",
        )
    );
}

#[test]
fn boundary_uses_absolute_residue_trace_and_hensel_lifts_split_roots() {
    let output = stdout(interpret("
        rvunram F[3,5]{1,0,1};
        rvcyclotomiclocal K[F,1];
        rvkummerlocal Unramified[K,3]([1]);
        rvkummerlocal Split[K,3]([0,1]);
        rvKummerLocalBoundary[Unramified];
        rvKummerLocalBoundary[Split];
        let theta = rvKummerLocalTheta[Split,1];
        let u = rvKummerLocalBase[Split,1](rvEisensteinElem[K,1]([0,1]));
        let one = rvKummerLocalElem[Split,1](([1]));
        let th2 = rvKummerLocalMul(theta,theta);
        assert rvKummerLocalMul(th2,theta) == add(theta,u);
        assert rvKummerLocalSigma[1](theta) == add(theta,one);
        assert rvKummerLocalSigma[1](rvKummerLocalSigma[2](theta)) == theta;
        let y = rvKummerLocalY[Split,4];
        let lp = rvKummerLocalBase[Split,4](rvEisensteinPi[K,4]);
        assert y == add(rvKummerLocalElem[Split,4](([1])),rvKummerLocalMul(lp,rvKummerLocalTheta[Split,4]));
        rvunram Low[3,1]{0,1}; rvcyclotomiclocal H[Low,1]; rvkummerlocal R[H,3]([1]);
        rvKummerLocalBoundary[R];"
    ));
    let lines: Vec<_> = output.lines().collect();
    assert!(
        lines[0].contains("trace=2, kind=unramified, residue_degree=3, different=0"),
        "{output}"
    );
    assert!(lines[0].ends_with("split_lifts=[]"), "{output}");
    assert!(lines[1].contains("trace=0, kind=split, residue_polynomial=[[0, 2], [2, 0], [0, 0], [1, 0]], residual_roots=[[0, 1], [2, 1], [1, 1]], theta_precision=8"), "{output}");
    assert!(
        lines[1].ends_with(
            "split_lifts=[[[39, 49], [40, 49]], [[44, 49], [2, 0]], [[1, 64], [40, 32]]]"
        ),
        "{output}"
    );
    assert!(lines[2].ends_with("theta_precision=0"), "{output}");
}

#[test]
fn insufficient_precision_bad_units_and_other_towers_are_diagnosed() {
    let setup = "
        rvunram Q[3,5]{0,1}; rvcyclotomiclocal K[Q,1]; rvcyclotomiclocal H[Q,1];
        rvcyclotomiclocal C[Q,2]; rvkummerlocal L[K,2]([1]); rvkummerlocal B[K,3]([1]);";
    assert_eq!(
        stdout(interpret(&format!("{setup} rvKummerLocalData[L,2]; rvKummerLocalPolynomial[L,3];"))),
        concat!(
            "rvKummerLocalData[L,2] = b=2, w_order=2, uniformizer=(w^2)/lambda^1, distances_at_least=2, unresolved\n",
            "rvKummerLocalPolynomial[L,3] = precision=0, unresolved\n",
        )
    );
    for (source, expected) in [
        ("rvkummerlocal T[K,0]([1]);", "unit level must lie in 1..=p"),
        (
            "rvkummerlocal T[K,2]([0],[1]);",
            "parameter u must be a unit",
        ),
        (
            "rvkummerlocal T[C,2]([1]);",
            "base must be a named zeta_p extension",
        ),
        (
            "show rvKummerLocalPi[L,25];",
            "exceeds the base digits after clearing",
        ),
        (
            "show rvKummerLocalTheta[L,5];",
            "theta generator requires b=p",
        ),
        (
            "show rvKummerLocalPi[B,3];",
            "requires a totally ramified unit level b<p",
        ),
        ("rvKummerLocalBoundary[L];", "boundary report requires b=p"),
        (
            "show rvKummerLocalBase[L,9](rvEisensteinPi[K,2]);",
            "RVEisenstein(K,3)",
        ),
        (
            "show rvKummerLocalBase[L,9](rvEisensteinPi[H,3]);",
            "RVEisenstein(K,3)",
        ),
        (
            "show rvKummerLocalEval(rvKummerLocalPi[L,24],rvTatePoly[3,3](1));",
            "precision is insufficient",
        ),
    ] {
        let result = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(result.contains(expected), "{result}");
    }
}
