mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

const BASE: &str = "
    rvunram B[3,7]{0,1};
    rvcyclotomiclocal K[B,1];
    let one = rvEisensteinElem[K,14]([1]);
    let l = rvEisensteinPi[K,14];
    let z = add(one,l);
    let l2 = rvEisensteinMul(l,l);
    let l3 = rvEisensteinMul(l2,l);
    let u2 = add(one,l2);
    let u3 = add(one,l3);
";

#[test]
fn wild_gram_signs_are_derived_and_export_typed_roots() {
    assert_eq!(
        stdout(example("wild-cubic-hilbert-symbol")),
        concat!(
            "rvCubicHilbertCertificate[K] = basis=[lambda, zeta, U2, U3], gram=[[0, 0, 0, 2], [0, 0, 2, 0], [0, 1, 0, 0], [1, 0, 0, 0]], determinant_mod3=1, boundary_frobenius=1, uniformizer_norms_mod_lambda5=[[0, 1], [0, 1]], norm_U2_mod27=13, zeta9_action=7\n",
            "rvCubicHilbertInvariant = exponent=2, invariant=2/3\n",
            "rvCubicNormSubgroup = functional=[1, 0, 0, 0], basis=[[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], order=27\n",
            "RVCubicClass(K){[1, 1, 0, 0]}\n",
            "RVCubicClass(K){[1, 1, 2, 1]}\n",
        )
    );

    let mut program = BASE.to_owned();
    let names = ["l", "z", "u2", "u3"];
    let matrix = [[0, 0, 0, 2], [0, 0, 2, 0], [0, 1, 0, 0], [1, 0, 0, 0]];
    for (i, &a) in names.iter().enumerate() {
        for (j, &b) in names.iter().enumerate() {
            writeln!(
                program,
                "assert rvCubicHilbert(rvCubicClass({a}),rvCubicClass({b})) == rvCyclo[3,{}];",
                matrix[i][j]
            )
            .unwrap();
        }
    }
    assert!(stdout(interpret(&program)).is_empty());
}

#[test]
fn cube_representatives_bilinearity_and_norm_kernels_are_composable() {
    let source = format!(
        r"{BASE}
         let a = rvCubicLocal[K,8,-7](u2);
         let b = rvCubicLocal[K,8,2](z);
         let p = rvCubicLocalMul(a,b);
         let b2 = rvCubicLocalMul(b,b);
         let cube = rvCubicLocalMul(b2,b);
         assert rvCubicClass(p) == add(rvCubicClass(a),rvCubicClass(b));
         assert rvCubicClass(cube) == zero[RVCubicClass(K)];
         assert rvCubicClass(rvCubicLocalMul(a,cube)) == rvCubicClass(a);
         assert rvCubicLocalInv(rvCubicLocalInv(a)) == a;
         assert rvCubicClass(rvCubicLocalInv(a)) == neg(rvCubicClass(a));
         assert rvCubicClass(add(one,rvEisensteinMul(l2,l2))) == zero[RVCubicClass(K)];
         assert rvCubicClass(neg(u3)) == rvCubicClass(u3);

         let l0 = rvCubicLocalFrom[8](l);
         assert l0 == rvCubicLocal[K,8,1](one);
         assert rvCubicLocalMul(l0,rvCubicLocalInv(l0)) == rvCubicLocal[K,8,0](one);

         let A = rvCubicClass(a); let V = rvCubicClass(u3); let W = rvCubicClass(b);
         assert rvCubicHilbert(add(A,V),W) == rvCycloMul(rvCubicHilbert(A,W),rvCubicHilbert(V,W));
         assert rvCubicHilbert(W,add(A,V)) == rvCycloMul(rvCubicHilbert(W,A),rvCubicHilbert(W,V));
         assert rvCycloMul(rvCubicHilbert(A,V),rvCubicHilbert(V,A)) == rvCyclo[3,0];
         assert rvCubicHilbert(A,A) == rvCyclo[3,0];
         assert rvCubicIsNorm(rvCubicClass(z),rvCubicClass(l)) == S.hi;
         assert rvCubicIsNorm(V,rvCubicClass(l)) == S.lo;
         card (\x : RVCubicClass(K) => rvCubicIsNorm(A,x));
         card (\x : RVCubicClass(K) => rvCubicAlgebraSplits(zero[RVCubicClass(K)],x));
         rvCubicNormSubgroup(zero[RVCubicClass(K)]);
         assert rvCubicClass[K](-2,3,7,-6) == rvCubicClass[K](1,0,1,0);"
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "card = 27 / 81\n",
            "card = 81 / 81\n",
            "rvCubicNormSubgroup = functional=[0, 0, 0, 0], basis=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], order=81\n",
        )
    );
}

#[test]
fn actual_steinberg_pairs_cover_integral_and_negative_valuations() {
    let mut program = BASE.to_owned();
    // Nine unit digits allow up to five digits of cancellation in 1-a.
    for a in -12..=12 {
        for b in -12..=12 {
            if [(0, 0), (1, 0)].contains(&(a, b)) {
                continue;
            }
            let value = format!("rvEisensteinElem[K,14]([{a}],[{b}])");
            writeln!(
                program,
                "assert rvCubicHilbert(rvCubicClass({value}),rvCubicClass(rvCubicOneMinus[4](rvCubicLocalFrom[9]({value})))) == rvCyclo[3,0];"
            )
            .unwrap();
        }
    }
    for valuation in -8..=-1 {
        for unit in ["one", "z", "u2", "u3"] {
            let value = format!("rvCubicLocal[K,8,{valuation}]({unit})");
            writeln!(
                program,
                "assert rvCubicIsNorm(rvCubicClass({value}),rvCubicClass(rvCubicOneMinus[4]({value}))) == S.hi;"
            )
            .unwrap();
        }
    }
    assert!(stdout(interpret(&program)).is_empty());
}

#[test]
fn two_named_bases_precision_and_nonzero_guards() {
    let header = format!(
        "{BASE} rvunram C[3,6]{{0,1}}; rvcyclotomiclocal H[C,1];
         rvcyclotomiclocal Tower[C,2]; rvunram E[3,6]{{1,0,1}}; rvcyclotomiclocal Wide[E,1];"
    );
    let source = format!(
        "{header}
         let p = rvEisensteinPi[H,7]; let t = rvEisensteinElem[H,7]([1]);
         let t3 = add(t,rvEisensteinMul(p,rvEisensteinMul(p,p)));
         assert rvCubicHilbert(rvCubicClass(t3),rvCubicClass(p)) ==
             rvCubicHilbert(rvCubicClass(u3),rvCubicClass(l));
         let huge = rvCubicLocal[K,4,-170141183460469231731687303715884105728](one);
         assert rvCubicIsNorm(rvCubicClass(huge),rvCubicClass(rvCubicOneMinus[4](huge))) == S.hi;"
    );
    assert!(stdout(interpret(&source)).is_empty());
    for (source, expected) in [
        ("rvCubicHilbertCertificate[Tower];", "requires a named Q_3(zeta_3)"),
        ("rvCubicHilbertCertificate[Wide];", "requires a named Q_3(zeta_3)"),
        ("show rvCubicClass(rvEisensteinPi[K,4]);", "insufficient relative precision"),
        ("show rvCubicClass(rvCubicLocal[K,3,0](one));", "at least four relative lambda digits"),
        ("show rvCubicClass(zero[RVEisenstein(K,8)]);", "cannot certify a nonzero"),
        ("show rvCubicLocal[K,4,-2](l);", "normalized cubic local factor must be a unit"),
        ("show rvCubicOneMinus[4](rvCubicLocal[K,4,0](one));", "cannot certify a nonzero"),
        ("show rvCubicHilbert(rvCubicClass[K](1,0,0,0),rvCubicClass[H](1,0,0,0));", "RVCubicClass(K)"),
        ("show rvCubicClass[K](1,0,0);", "requires four basis coordinates"),
        ("show rvCubicLocalInv(rvCubicLocal[K,4,-170141183460469231731687303715884105728](one));", "signed cubic local valuation exceeds"),
    ] {
        let result = diagnostic(interpret(&format!("{header}{source}")));
        assert!(result.contains(expected), "{result}");
    }
}
