mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const HUGE: u128 = u128::MAX;
const SYMBOLS: &str = r"
    let s1 = rvSchurCharacter[[1]];
    let s2 = rvSchurCharacter[[2]];
    let s11 = rvSchurCharacter[[1,1]];
    let one = rvSchurCharacter[[]];
";

#[test]
fn runnable_example_prints_stable_and_virtual_characters() {
    assert_eq!(
        stdout(example("witt-schur-plethysm")),
        concat!(
            "RVSchurCharacter(4){([4],1), ([2, 2],1)}\n",
            "rvSchurCharacterCoeff[[2, 2]] = 1\n",
            "RVSchurCharacter(4){([4],1), ([3, 1],-1), ([2, 2],1)}\n",
            "rvWittSpectrumGhost[1] = 21\n",
        )
    );
}

#[test]
fn stable_degree_four_and_higher_identities_hold() {
    let source = format!(
        r"{SYMBOLS}
        let s4 = rvSchurCharacter[[4]];
        let s31 = rvSchurCharacter[[3,1]];
        let s22 = rvSchurCharacter[[2,2]];
        let s211 = rvSchurCharacter[[2,1,1]];
        let s1111 = rvSchurCharacter[[1,1,1,1]];
        assert rvSchurPlethysm(s2,s2) == add(s4,s22);
        assert rvSchurPlethysm(s11,s2) == s31;
        assert rvSchurPlethysm(s2,s11) == add(s22,s1111);
        assert rvSchurPlethysm(s11,s11) == s211;
        let full = add(add(s4,s31),add(rvSchurCharacterScale[2](s22),add(s211,s1111)));
        assert rvSchurPlethysm(s2,add(s2,s11)) == full;
        assert rvSchurPlethysm(s2,add(s2,neg(s11))) == add(add(s4,neg(s31)),s22);
        assert rvSchurPlethysm(rvSchurCharacter[[3]],s2)
            == add(rvSchurCharacter[[6]],add(rvSchurCharacter[[4,2]],rvSchurCharacter[[2,2,2]]));
        assert rvSchurPlethysm(rvSchurCharacter[[2,1]],s2)
            == add(rvSchurCharacter[[5,1]],add(rvSchurCharacter[[4,2]],rvSchurCharacter[[3,2,1]]));
        rvSchurCharacterCoeff[[1,1,1,1]](rvSchurPlethysm(s2,s2));
        rvSchurCharacterCoeff[[3,2,1]](rvSchurPlethysm(rvSchurCharacter[[2,1]],s2));
    "
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvSchurCharacterCoeff[[1, 1, 1, 1]] = 0\nrvSchurCharacterCoeff[[3, 2, 1]] = 1\n"
    );
}

fn power_sum(n: usize) -> String {
    (0..n)
        .map(|j| {
            let mut hook = vec![1; j + 1];
            hook[0] = n - j;
            let term = format!("rvSchurCharacter[{hook:?}]");
            if j.is_multiple_of(2) {
                term
            } else {
                format!("neg({term})")
            }
        })
        .reduce(|a, b| format!("add({a},{b})"))
        .unwrap()
}

#[test]
fn power_sum_plethysm_composes_indices_with_signed_rim_hooks() {
    let mut source = String::new();
    for (a, b) in [(2, 2), (3, 4), (4, 4), (2, 7)] {
        writeln!(
            source,
            "assert rvSchurPlethysm({},{}) == {};",
            power_sum(a),
            power_sum(b),
            power_sum(a * b)
        )
        .unwrap();
    }
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn associativity_outer_linearity_and_native_witt_evaluation_agree() {
    let source = format!(
        r"{SYMBOLS}
        let f = add(s2,neg(s11));
        let g = add(s2,rvSchurCharacterScale[-2](s11));
        assert rvSchurPlethysm(rvSchurPlethysm(f,g),s11) == rvSchurPlethysm(f,rvSchurPlethysm(g,s11));
        assert rvSchurPlethysm(rvSchurPlethysm(s2,g),s11) == rvSchurPlethysm(s2,rvSchurPlethysm(g,s11));
        assert rvSchurPlethysm(f,g) == add(rvSchurPlethysm(s2,g),neg(rvSchurPlethysm(s11,g)));
        assert rvSchurPlethysm(rvSchurCharacterScale[-3](f),g) == rvSchurCharacterScale[-3](rvSchurPlethysm(f,g));
        let x = rvWittSpectrum{{[1,-3,2]}};
        let evaluated = rvSchurCharacterEvaluate(g,x);
        assert rvSchurCharacterEvaluate(rvSchurPlethysm(f,g),x) == rvSchurCharacterEvaluate(f,evaluated);
        assert rvSchurCharacterEvaluate(rvSchurPlethysm(f,g),x) == rvWittSpectrumAdams[2](evaluated);
        let h = rvSchurCharacter[[2,1]];
        let w = rvWittSpectrum{{[1,-2,1]}};
        assert rvSchurCharacterEvaluate(rvSchurPlethysm(h,h),w)
            == rvSchurCharacterEvaluate(h,rvSchurCharacterEvaluate(h,w));
        assert rvSchurCharacterEvaluate(rvSchurPlethysm(s2,g),w)
            == rvWittSpectrumSchur[[2]](rvSchurCharacterEvaluate(g,w));
    "
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn constants_and_homogeneous_zero_use_the_integral_binomial_lambda_ring() {
    let source = format!(
        r"{SYMBOLS}
        let blank = zero[RVSchurCharacter(0)];
        let minus3 = rvSchurCharacterScale[-3](one);
        assert rvSchurPlethysm(one,zero[RVSchurCharacter(91)]) == one;
        assert rvSchurPlethysm(blank,s2) == blank;
        assert rvSchurPlethysm(s2,blank) == blank;
        assert rvSchurPlethysm(s2,zero[RVSchurCharacter(8)]) == zero[RVSchurCharacter(16)];
        assert rvSchurPlethysm(zero[RVSchurCharacter(17)],s2) == zero[RVSchurCharacter(34)];
        assert rvSchurPlethysm(s1,zero[RVSchurCharacter({HUGE})]) == zero[RVSchurCharacter({HUGE})];
        assert rvSchurPlethysm(zero[RVSchurCharacter({HUGE})],s1) == zero[RVSchurCharacter({HUGE})];
        assert rvSchurPlethysm(s2,minus3) == rvSchurCharacterScale[3](one);
        assert rvSchurPlethysm(s11,minus3) == rvSchurCharacterScale[6](one);
        assert rvSchurPlethysm(rvSchurCharacter[[2,1]],rvSchurCharacterScale[3](one)) == rvSchurCharacterScale[8](one);
        assert rvSchurPlethysm(rvSchurCharacter[[1,1,1]],neg(one)) == neg(one);
        let huge = rvSchurCharacter[[{HUGE}]];
        assert rvSchurPlethysm(huge,one) == one;
        assert rvSchurPlethysm(huge,neg(one)) == blank;
        assert rvSchurPlethysm(s1,huge) == huge;
        assert rvSchurCharacterScale[0](huge) == zero[RVSchurCharacter({HUGE})];
        assert add(huge,neg(huge)) == zero[RVSchurCharacter({HUGE})];
        assert rvSchurCharacterEvaluate(one,zero[RVWittSpectrum]) == rvWittSpectrumLine[1];
        assert rvSchurCharacterEvaluate(minus3,zero[RVWittSpectrum]) == neg(rvWittSpectrum{{[1,-3,3,-1]}});
    "
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn incompatible_types_shapes_and_bounds_have_direct_diagnostics() {
    for (query, expected) in [
        ("show add(s1,s11);", "RVSchurCharacter"),
        ("rvSchurCharacterCoeff[[1,1]](s1);", "RVSchurCharacter(2)"),
        ("show rvSchurCharacter[[1,3]];", "weakly decreasing"),
        ("show rvSchurPlethysm(s1,S.hi);", "requires an RVSchurCharacter"),
        ("show rvSchurCharacterEvaluate(s1,S.hi);", "RVWittSpectrum"),
        ("show rvSchurCharacterScale[4](S.hi);", "requires an RVSchurCharacter"),
        ("show rvSchurPlethysm(s2,rvSchurCharacter[[9]]);", "weight limit of 16 boxes"),
        ("show rvSchurPlethysm(rvSchurCharacter[[340282366920938463463374607431768211455]],s2);", "weight exceeds unsigned 128-bit range"),
        ("show rvSchurPlethysm(s2,rvSchurCharacterScale[170141183460469231731687303715884105727](s1));", "coefficient exceeds signed 128-bit range"),
        ("show neg(rvSchurCharacterScale[-170141183460469231731687303715884105728](s1));", "coefficient exceeds signed 128-bit range"),
        ("show rvSchurCharacterEvaluate(rvSchurCharacterScale[257](s1),rvWittSpectrumLine[1]);", "determinant degree limit"),
    ] {
        let actual = diagnostic(interpret(&format!("{SYMBOLS} {query}")));
        assert!(actual.contains(expected), "{actual}");
    }
}
