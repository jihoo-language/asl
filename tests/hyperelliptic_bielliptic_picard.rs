mod common;

use common::{diagnostic, example, interpret, stdout};

const F9: &str = r"
    rvunram F[3,1]{1,0,1};
    let z=rvUnramElem[F]{0};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let i=rvUnramElem[F]{0,1};
    let a=rvUnramElem[F]{1,1};
    let c=rvUnramElem[F]{2,1};
    rvhyperelliptic C[F]{a,z,neg(a),z,o,o};
    rvbielliptic C[Plus,Minus](z,i,o,two);
";

#[test]
fn split_picard_compositions_hold_for_all_classes_in_the_even_curve_chart() {
    assert_eq!(
        stdout(example("hyperelliptic-bielliptic-picard")),
        "rvJacCount[S] = 64\nrvJacCount[Plus] = 8\nrvJacCount[Minus] = 8\n"
    );
}

#[test]
fn quintic_push_corrects_the_basepoint_and_traces_nonsplit_and_repeated_support() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F9}
            let p=rvJacPoint[C](o,i);
            let plus=rvBiellipticPoint[Plus](o,i);
            let minus=rvBiellipticPoint[Minus](o,i);
            let plusBase=rvJacPoint[Plus](a,z);
            let minusBase=rvJacPoint[Minus](i,z);
            assert rvBiellipticPush[Plus](p)==add(plus,plusBase);
            assert rvBiellipticPush[Minus](p)==add(minus,minusBase);

            let tangent=rvJacScale[2](p);
            assert rvJacU[2](tangent)==o;
            assert rvBiellipticPush[Plus](tangent)==rvJacScale[2](plus);
            assert rvBiellipticPush[Minus](tangent)==rvJacScale[2](minus);

            // u=x²+i*x+i is irreducible over F9. The secant through each
            // conjugate image pair yields the written sum over F9.
            let closed=rvjac[C]{{[0,1],[0,1],[1]}}{{[],[1]}};
            assert rvBiellipticPush[Plus](closed)==rvJacPoint[Plus](z,o);
            assert rvBiellipticPush[Minus](closed)==rvJacPoint[Minus](neg(a),a);

            let pole=rvJacPoint[C](neg(i),c);
            assert rvBiellipticPush[Minus](pole)==minusBase;
            assert rvJacU[2](rvJacScale[2](pole))==o;
            assert rvBiellipticPush[Minus](rvJacScale[2](pole))==zero[RVJac(Minus)];
            assert rvBiellipticPull[C](rvJacPoint[Plus](z,two))==rvJacScale[2](pole);

            // Each basepoint fiber contains infinity and the branch x=i.
            let branch=rvJacPoint[C](i,z);
            assert rvBiellipticPull[C](plusBase)==branch;
            assert rvBiellipticPull[C](minusBase)==branch;
            assert rvBiellipticPull[C](zero[RVJac(Plus)])==zero[RVJac(C)];
            assert rvBiellipticPull[C](zero[RVJac(Minus)])==zero[RVJac(C)];
            "
        ))),
        ""
    );
}

#[test]
fn both_picard_maps_commute_with_changing_the_branch_chart() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F9}
            // The branch x=2 has derivative -i, so the chart rescales y
            // by a genuine extension-field coefficient.
            rvjacchart S[C](two);
            rvJacCount[C];
            rvJacCount[Plus];
            rvJacCount[Minus];
            let pushPlus=\d:RVJac(C)=>rvBiellipticPush[Plus](d);
            let pushMinus=\d:RVJac(C)=>rvBiellipticPush[Minus](d);
            let pullPlus=\p:RVJac(Plus)=>rvBiellipticPull[C](p);
            let pullMinus=\q:RVJac(Minus)=>rvBiellipticPull[C](q);
            assert (\d:RVJac(C)=>rvBiellipticPush[Plus](rvJacRechart[S](d)))==pushPlus;
            assert (\d:RVJac(C)=>rvBiellipticPush[Minus](rvJacRechart[S](d)))==pushMinus;
            assert (\p:RVJac(Plus)=>rvBiellipticPull[S](p))
                ==(\p:RVJac(Plus)=>rvJacRechart[S](pullPlus(p)));
            assert (\q:RVJac(Minus)=>rvBiellipticPull[S](q))
                ==(\q:RVJac(Minus)=>rvJacRechart[S](pullMinus(q)));
            assert (\p:RVJac(Plus)=>pushPlus(pullPlus(p)))
                ==(\p:RVJac(Plus)=>rvJacScale[2](p));
            assert (\q:RVJac(Minus)=>pushMinus(pullMinus(q)))
                ==(\q:RVJac(Minus)=>rvJacScale[2](q));
            assert (\p:RVJac(Plus)=>pushMinus(pullPlus(p)))
                ==(\p:RVJac(Plus)=>zero[RVJac(Minus)]);
            assert (\q:RVJac(Minus)=>pushPlus(pullMinus(q)))
                ==(\q:RVJac(Minus)=>zero[RVJac(Plus)]);
            assert (\d:RVJac(C)=>add(pullPlus(pushPlus(d)),pullMinus(pushMinus(d))))
                ==(\d:RVJac(C)=>rvJacScale[2](d));
            "
        ))),
        "rvJacCount[C] = 96\nrvJacCount[Plus] = 8\nrvJacCount[Minus] = 12\n"
    );
}

#[test]
fn picard_maps_require_the_actual_nominal_factor_and_source() {
    for (input, expected) in [
        (
            "rvhyperelliptic Twin[F]{a,z,neg(a),z,o,o}; show rvBiellipticPush[Plus](zero[RVJac(Twin)]);",
            "must be on the named bielliptic source",
        ),
        (
            "rvhyperelliptic Twin[F]{a,z,neg(a),z,o,o}; rvjacchart S[Twin](i); show rvBiellipticPull[S](zero[RVJac(Plus)]);",
            "must be on the named bielliptic source",
        ),
        (
            "rvhyperelliptic Twin[F]{o,a,a,o}; show rvBiellipticPull[C](zero[RVJac(Twin)]);",
            "has no stored bielliptic quotient",
        ),
        (
            "show rvBiellipticPush[C](zero[RVJac(C)]);",
            "has no stored bielliptic quotient",
        ),
        ("show rvBiellipticPull[C](z);", "an RVJac divisor class"),
        (
            "show add(rvBiellipticPush[Minus](zero[RVJac(C)]),zero[RVJac(Plus)]);",
            "expected type RVJac(Minus)",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F9}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
