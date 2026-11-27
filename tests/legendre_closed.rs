mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn degree_two_unit_and_actual_semilinear_factors() {
    assert_eq!(
        stdout(example("legendre-closed-points")),
        concat!(
            "RVUnram(W)[3039, 49]\n",
            "RVUnram(W)[3039, 3076]\n",
            "RVUnram(W)[2823, 0]\n",
            "rvLegendreTrace = q=25, trace=-2, points=28\n",
            "rvLegendreHensel mod 3125 = q=25, trace=-2, unit=2823\n",
        )
    );
}

#[test]
fn degree_three_factors_use_the_signed_ratio_and_the_true_frobenius() {
    let source = r"
        rvunram W[3,4]{1,-1,0,1};
        let x = rvUnramElem[W]{0,1};
        let z = rvLegendreParameter(x);
        assert z == rvUnramElem[W]{54,49,21};
        assert rvLegendreParameter(add(x,rvUnramElem[W]{3})) == z;
        assert rvLegendreFactor[0](x) == rvUnramElem[W]{23,23,42};
        assert rvLegendreFactor[1](x) == rvUnramElem[W]{78,80,0};
        assert rvLegendreFactor[2](x) == rvUnramElem[W]{52,59,39};
        assert rvLegendreFactor[3](x) == rvLegendreFactor[0](x);
        assert rvLegendreFactor[0](rvUnramFrob[1](x)) == rvLegendreFactor[1](x);
        assert rvLegendreUnit(x) == rvUnramElem[W]{58};
        rvLegendreTrace(x);
        rvLegendreHensel(x);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvLegendreTrace = q=27, trace=4, points=24\n",
            "rvLegendreHensel mod 81 = q=27, trace=4, unit=58\n",
        )
    );
}

#[test]
fn degree_one_and_proper_subfield_points_give_the_expected_power() {
    let source = r"
        rvunram K[5,5]{-2,1};
        let lambda = rvUnramElem[K]{2};
        assert rvUnramFrob[1](lambda) == lambda;
        assert rvLegendreUnit(lambda) == rvUnramElem[K]{988};

        rvunram Precise[5,8]{-2,1};
        assert rvLegendreUnit(rvUnramElem[Precise]{2}) == rvUnramElem[Precise]{169738};
        rvLegendreHensel(lambda);
        rvunram W[5,5]{2,0,1};
        let ext = rvUnramElem[W]{2};
        assert rvLegendreFactor[0](ext) == rvUnramElem[W]{988};
        assert rvLegendreFactor[0](ext) == rvLegendreFactor[1](ext);
        assert rvLegendreUnit(ext) == rvUnramElem[W]{1144};
        assert rvLegendreUnit(rvUnramElem[W]{0,2}) == rvUnramElem[W]{302};
        rvLegendreHensel(ext);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvLegendreHensel mod 3125 = q=5, trace=-2, unit=988\n",
            "rvLegendreHensel mod 3125 = q=25, trace=-6, unit=1144\n",
        )
    );
}

#[test]
fn actual_unramified_ring_arithmetic_and_finite_values() {
    let source = r"
        rvunram W[5,5]{2,-1,1};
        let x = rvUnramElem[W]{0,1};
        assert add(x,neg(x)) == zero[RVUnram(W)];
        assert rvUnramTrace(x) == rvUnramElem[W]{1};
        assert rvUnramNorm(x) == rvUnramElem[W]{2};
        assert rvUnramMul(x,rvUnramInv(x)) == rvUnramElem[W]{1};
        let t = rvUnramTeich(x);
        assert rvUnramTeich(t) == t;
        assert rvUnramPow[25](t) == t;
        assert rvUnramFrob[1](t) == rvUnramPow[5](t);
        rvUnramVal(rvUnramElem[W]{25,125});
        rvUnramVal(zero[RVUnram(W)]);
        rvunram F[2,1]{1,1,1};
        let a = rvUnramElem[F]{0,1};
        assert rvUnramFrob[1](a) == rvUnramPow[2](a);
        assert rvUnramNorm(a) == rvUnramElem[F]{1};
        card(\r : RVUnram(F) => S.hi);
        rvunram Big[2,63]{1,1};
        assert rvUnramPow[2](rvUnramElem[Big]{-1}) == rvUnramElem[Big]{1};
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvUnramVal = 2\nrvUnramVal = >=5\ncard = 4 / 4\n"
    );
}

#[test]
fn supersingular_period_valuations_are_reported_before_division() {
    let source = r"
        rvunram W[3,3]{1,0,1};
        let lambda = rvUnramElem[W]{2};
        rvLegendreData(lambda);
        rvLegendreTrace(lambda);
        rvLegendreHensel(lambda);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvLegendreData mod 27 = z=[26, 0], B_3=[9, 0], B_2=[24, 0], valuations=(2, 1), supersingular\n",
            "rvLegendreTrace = q=9, trace=-6, points=16\n",
            "rvLegendreHensel mod 27 = q=9, trace=-6, supersingular\n",
        )
    );
}

#[test]
fn unramified_declarations_and_partial_operations_check_their_domains() {
    let setup = "rvunram W[3,3]{1,0,1}; rvunram V[3,3]{4,0,1}; let x = rvUnramElem[W]{0,1};";
    for (code, reason) in [
        ("rvunram B[9,2]{1,1};", "9 is not prime"),
        ("rvunram B[5,0]{1,1};", "precision must be positive"),
        ("rvunram B[5,64]{1,1};", "p^s exceeds u64"),
        ("rvunram B[5,2]{1,0,2};", "must be monic"),
        ("rvunram B[5,2]{1,0,1};", "irreducible modulo p"),
        ("rvunram W[3,3]{1,0,1};", "already defined"),
        ("show rvUnramElem[W]{1,2,3};", "at most 2 coordinates"),
        ("show rvUnramInv(rvUnramElem[W]{3});", "requires a unit"),
        ("show rvUnramMul(x,rvUnramElem[V]{1});", "RVUnram(W)"),
        (
            "rvLegendreTrace(rvUnramElem[W]{4});",
            "parameters 0 and 1 are singular",
        ),
        ("show rvLegendreUnit(rvUnramElem[W]{2});", "supersingular"),
        (
            "rvunram B[2,1]{1,1,1}; rvLegendreTrace(rvUnramElem[B]{0,1});",
            "requires an odd residue characteristic",
        ),
    ] {
        let output = diagnostic(interpret(&format!("{setup}\n{code}")));
        assert!(output.contains(reason), "{code}\n{output}");
    }
}
