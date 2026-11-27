mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

const BASE: &str = "
    rvunram B[3,7]{0,1};
    rvcyclotomiclocal K[B,1];
    let l = rvEisensteinPi[K,14];
    let one = rvEisensteinElem[K,14]([1]);
    let z = add(one,l);
    let l2 = rvEisensteinMul(l,l);
    let u2 = add(one,l2);
    let u3 = add(one,rvEisensteinMul(l2,l));
";

#[test]
fn division_example_has_the_expected_reduced_polynomial_and_unit_inverse() {
    let output = stdout(example("wild-cubic-central-simple-algebras"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        &lines[..5],
        [
            "rvCubicAlgData[D] = base=K, precision=8, parameters=[[1, 1], [79, 78]], shifts=[0, 0], invariant=2/3, kind=division",
            "rvCubicReducedChar = precision=8, coefficients=[[1, 0], [78, 0], [3, 0], [64, 9]]",
            "RVEisenstein(K,8){[[17], [72]]}",
            "rvCubicAlgData[A] = base=K, precision=8, parameters=[[1, 1], [0, 80]], shifts=[0, 0], invariant=0, kind=split",
            "rvCubicSplitData[S] = algebra=A, precision=8, norm_c=[0, 80], relations=verified",
        ]
    );
    assert!(
        lines[5].starts_with("RVCubicSplitMatrix(S){[[[0, 0]"),
        "{output}"
    );

    let program = format!(
        "{BASE}
         rvcubicalgebra A[K,8](z,u2);
         let u = rvCubicAlgBasis[A,1,0]; let v = rvCubicAlgBasis[A,0,1];
         let zz = rvEisensteinElem[K,8]([1],[1]);
         assert rvCubicAlgCoeff[1,1](rvCubicAlgMul(v,u)) == zz;
         assert rvCubicAlgMul(v,u) == rvCubicAlgTerm[A,1,1](zz);
         let x = rvCubicAlgElem[A]([1],[2],[],[1],[],[],[],[],[1]);
         let y = add(u,v);
         assert rvCubicReducedTrace(x) == rvEisensteinElem[K,8]([3]);
         let nx = rvCubicReducedNorm(x);
         assert rvCubicReducedNorm(rvCubicAlgMul(x,y)) == rvEisensteinMul(nx,rvCubicReducedNorm(y));
         assert rvCubicAlgMul(rvCubicAdj(x),x) == rvCubicAlgTerm[A,0,0](nx);
         assert rvCubicAlgMul(rvCubicInverse(x),x) == rvCubicAlgBasis[A,0,0];
         assert add(x,neg(x)) == zero[RVCubicAlg(A)];
         assert rvCubicAlgIsSplit[A] == S.lo;
         rvcubicalgebra L[K,1](z,u2);
         let small = rvCubicAlgElem[L]([1],[2],[],[1],[],[],[],[],[1]);
         rvCubicReducedChar(small);"
    );
    assert_eq!(
        stdout(interpret(&program)),
        "rvCubicReducedChar = precision=1, coefficients=[[1, 0], [0, 0], [0, 0], [1, 0]]
"
    );
}

#[test]
fn steinberg_witness_gives_actual_matrix_products_and_nonzero_zero_divisors() {
    let mut program = BASE.to_owned();
    for (i, a) in ["z", "u2", "u3", "l"].into_iter().enumerate() {
        writeln!(
            program,
            "rvcubicalgebra A{i}[K,8]({a},add(one,neg({a}))); rvcubicsplit S{i}[A{i}]([1],[-1]);"
        )
        .unwrap();
        writeln!(program, "let e{i}=rvCubicAlgBasis[A{i},0,0]; let u{i}a=rvCubicAlgBasis[A{i},1,0]; let v{i}=rvCubicAlgBasis[A{i},0,1];").unwrap();
        writeln!(
            program,
            "let x{i}=add(e{i},neg(add(u{i}a,v{i}))); let adj{i}=rvCubicAdj(x{i});"
        )
        .unwrap();
        writeln!(program, "assert rvCubicAlgIsSplit[A{i}] == S.hi;").unwrap();
        writeln!(
            program,
            "assert rvCubicAlgMul(x{i},adj{i}) == zero[RVCubicAlg(A{i})];"
        )
        .unwrap();
        writeln!(
            program,
            "assert rvCubicReducedNorm(x{i}) == zero[RVEisenstein(K,8)];"
        )
        .unwrap();
        writeln!(program, "show adj{i};").unwrap();
        writeln!(program, "let m{i}=rvCubicSplitMat[S{i}](x{i}); let n{i}=rvCubicSplitMat[S{i}](add(v{i},u{i}a));").unwrap();
        writeln!(program, "assert m{i} == rvCubicSplitMat[S{i}](x{i});").unwrap();
        writeln!(program, "assert rvCubicSplitMul(m{i},n{i}) == rvCubicSplitMat[S{i}](rvCubicAlgMul(x{i},add(v{i},u{i}a)));").unwrap();
        writeln!(
            program,
            "assert rvCubicSplitDet(n{i}) == rvCubicReducedNorm(add(v{i},u{i}a));"
        )
        .unwrap();
        for row in 0..3 {
            writeln!(
                program,
                "assert rvCubicSplitEntry[{row},0](m{i}) == zero[RVEisenstein(K,8)];"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&program));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 4);
    for line in lines {
        assert!(
            !line.ends_with(
                "{[[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]}"
            ),
            "{output}"
        );
    }
}

#[test]
fn signed_valuations_change_generators_by_actual_cube_scales() {
    let program = format!(
        "{BASE}
        let a = rvCubicLocal[K,10,-5](z);
        let b = rvCubicLocal[K,10,-2](u2);
        rvcubicalgebra A[K,8](a,b);
        rvCubicAlgData[A];
        let u = rvCubicAlgBasis[A,1,0]; let v = rvCubicAlgBasis[A,0,1];
        let u_cubed = rvCubicAlgMul(rvCubicAlgMul(u,u),u);
        let v3 = rvCubicAlgMul(rvCubicAlgMul(v,v),v);
        assert rvCubicAlgCoeff[0,0](u_cubed) == rvEisensteinElem[K,8]([-3],[-2]);
        assert rvCubicAlgCoeff[0,0](v3) == rvEisensteinElem[K,8]([9],[7]);
        rvCubicHilbertInvariant(rvCubicClass(a),rvCubicClass(b));"
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvCubicAlgData[A] = base=K, precision=8, parameters=[[78, 79], [9, 7]], shifts=[2, 1], invariant=2/3, kind=division
",
            "rvCubicHilbertInvariant = exponent=2, invariant=2/3
",
        )
    );
}

#[test]
fn chart_precision_and_nominal_parameters_cannot_forge_inverses_or_witnesses() {
    let prefix = format!(
        "{BASE}
         rvcubicalgebra D[K,8](z,u2);
         rvcubicalgebra A[K,8](l,add(one,neg(l))); rvcubicsplit S[A]([1],[-1]);
         rvcubicalgebra Low[K,1](l,add(one,neg(l)));
         rvcyclotomiclocal J[B,1];"
    );
    for (source, expected) in [
        (
            "rvcubicalgebra Bad[K,8](rvEisensteinElem[J,12]([1]),z);",
            "named cubic local base 'K'",
        ),
        (
            "rvcubicalgebra Bad[K,8](rvCubicLocal[K,6,0](z),u2);",
            "insufficient relative precision",
        ),
        (
            "rvcubicalgebra Bad[K,1](rvCubicLocal[K,3,0](z),u2);",
            "insufficient relative precision",
        ),
        (
            "rvcubicalgebra Bad[K,8](zero[RVEisenstein(K,12)],z);",
            "cannot certify a nonzero",
        ),
        (
            "rvcubicsplit Bad[D]([1],[-1]);",
            "division algebra cannot admit",
        ),
        ("rvcubicsplit Bad[A]([1]);", "norm differs from b"),
        (
            "show rvCubicInverse(rvCubicAlgBasis[A,1,0]);",
            "positive valuation",
        ),
        (
            "show rvCubicInverse(rvCubicAlgBasis[Low,1,0]);",
            "norm vanishes at the available precision",
        ),
        (
            "show rvCubicAlgMul(rvCubicAlgBasis[A,0,1],rvCubicAlgBasis[D,0,1]);",
            "RVCubicAlg(A)",
        ),
        (
            "show rvCubicSplitMat[S](rvCubicAlgBasis[D,0,1]);",
            "RVCubicAlg(A)",
        ),
        ("show rvCubicAlgBasis[A,3,0];", "indices must be in 0..=2"),
        (
            "show rvCubicAlgTerm[A,0,0](rvEisensteinElem[K,7]([1]));",
            "RVEisenstein(K,8)",
        ),
    ] {
        let report = diagnostic(interpret(&format!("{prefix}{source}")));
        assert!(report.contains(expected), "{report}");
    }
}
