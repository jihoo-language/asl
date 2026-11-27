mod common;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zeroW=rvUnramElem[W]{};
    let oneW=rvUnramElem[W]{1};
    let a=rvUnramElem[W]{0,1};
    let r0=rvJacPoint[L](zeroW,zeroW);
    let r1=rvJacPoint[L](oneW,zeroW);
    let r2=rvJacPoint[L](neg(oneW),zeroW);
    let ra=rvJacPoint[L](a,zeroW);
    rvjactorsionbasis B[L,2](neg(oneW)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),ra)};
    let v=rvNakajimaVac[B];
";

#[test]
fn example_covers_the_heisenberg_bracket_and_the_wick_sign() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-nakajima")),
        concat!(
            "rvNakajimaPair = 1\n",
            "rvNakajimaPair = -2\n",
            "rvNakajimaPair = -48\n",
            "rvNakajimaPair = 5\n"
        )
    );
}

#[test]
fn creators_commute_and_different_modes_have_zero_commutator() {
    let output = stdout(interpret(&format!(
        r"{MARKED}
        assert rvNakajimaAnn[1,0](v)==zero[RVNakajima(B,-1,-4)];
        assert rvNakajimaAnn[1,23](v)==zero[RVNakajima(B,-1,0)];
        assert rvNakajimaCreate[1,0](rvNakajimaCreate[2,3](v))
            ==rvNakajimaCreate[2,3](rvNakajimaCreate[1,0](v));
        let x=rvNakajimaCreate[1,2](v);
        assert rvNakajimaAnn[1,1](rvNakajimaCreate[2,2](x))
            ==rvNakajimaCreate[2,2](rvNakajimaAnn[1,1](x));
        let y=rvNakajimaCreate[1,0](v);
        assert add(rvNakajimaAnn[1,23](rvNakajimaCreate[1,0](y)),neg(rvNakajimaCreate[1,0](rvNakajimaAnn[1,23](y))))==y;
        assert (\x:RVNakajima(B,2,1)=>x)==(\x:RVNakajima(B,2,1)=>zero[RVNakajima(B,2,1)]);
        assert rvNakajimaCreate[1,23](zero[RVNakajima(B,1,5)])==zero[RVNakajima(B,2,9)];
        "
    )));
    assert_eq!(output, "");
}

#[test]
fn adjunction_uses_the_marked_gram_and_complementary_grades() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            rvK3LatticePair(rvK3LatticeBasis[B,0],rvK3LatticeBasis[B,1]);
            let x=rvNakajimaCreate[1,0](v);
            let y=rvNakajimaMonomial[B]{{(2,2,1),(1,23,1)}};
            rvNakajimaPair(rvNakajimaCreate[2,1](x),y);
            rvNakajimaPair(x,rvNakajimaAnn[2,1](y));
            rvNakajimaPair(v,v);
            rvNakajimaPair(x,v);
            rvNakajimaPair(x,x);
            rvNakajimaPair(rvNakajimaCreate[2,0](v),rvNakajimaMonomial[B]{{(1,0,1),(1,23,1)}});
            "
        ))),
        concat!(
            "rvK3LatticePair = -1\n",
            "rvNakajimaPair = 2\n",
            "rvNakajimaPair = -2\n",
            "rvNakajimaPair = 1\n",
            "rvNakajimaPair = 0\n",
            "rvNakajimaPair = 0\n",
            "rvNakajimaPair = 0\n",
        )
    );
}

#[test]
fn colored_multiplicities_are_counted_by_the_full_permanent() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let p=rvNakajimaMonomial[B]{{(2,1,2),(2,2,1)}};
            let q=rvNakajimaMonomial[B]{{(2,2,2),(2,1,1)}};
            rvNakajimaPair(p,q);
            rvNakajimaPair(q,p);
            rvNakajimaPair(add(p,rvNakajimaScale[-2](p)),q);
            rvNakajimaPair(rvNakajimaMonomial[B]{{(1,0,4)}},rvNakajimaMonomial[B]{{(1,23,4)}});
            "
        ))),
        "rvNakajimaPair = 144\nrvNakajimaPair = 144\nrvNakajimaPair = -144\nrvNakajimaPair = 24\n"
    );
}

#[test]
fn enormous_exponents_and_wide_cancellation_stay_sparse() {
    let count = 1_000_000_000_000_000_000_000_000u128;
    let previous = count - 1;
    let heavy = (1i128 << 126) + 1;
    let max = i128::MAX;
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let x=rvNakajimaMonomial[B]{{(1,0,{previous}),(1,0,1),(3,23,0)}};
            let y=rvNakajimaMonomial[B]{{(1,0,{count})}};
            assert x==y;
            assert add(x,neg(y))==zero[RVNakajima(B,{count},0)];
            assert rvNakajimaAnn[1,23](x)==rvNakajimaScale[{count}](rvNakajimaMonomial[B]{{(1,0,{previous})}});
            let positive=rvNakajimaScale[{heavy}](rvNakajimaCreate[1,1](v));
            let negative=rvNakajimaScale[-{max}](rvNakajimaCreate[1,2](v));
            assert rvNakajimaAnn[1,1](add(positive,negative))==rvNakajimaScale[-3](v);
            rvNakajimaPair(rvNakajimaMonomial[B]{{(2,0,{count})}},rvNakajimaMonomial[B]{{(1,23,{count}),(1,1,{count})}});
            show x;
            "
        ))),
        format!("rvNakajimaPair = 0\nRVNakajima(B,{count},0){{1*[(1,0,{count})]}}\n")
    );
}

#[test]
fn rejects_wrong_markings_invalid_modes_and_exact_range_overflows() {
    for (source, expected) in [
        ("show rvNakajimaVac[Missing];", "unknown Jacobian torsion basis"),
        ("rvhyperelliptic E[F]{z,neg(o),z,o}; rvjactorsionbasis Wrong[E,2](neg(o)){rvJacPoint[E](z,z),rvJacPoint[E](o,z)}; show rvNakajimaVac[Wrong];", "require a certified genus-two level-2 basis"),
        ("rvjactorsionbasis D[L,2](neg(oneW)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),ra)}; rvNakajimaPair(v,rvNakajimaVac[D]);", "same named level basis"),
        ("show add(v,rvNakajimaCreate[1,0](v));", "expected type RVNakajima(B,0,0)"),
        ("show rvNakajimaCreate[0,0](v);", "mode length must be positive"),
        ("show rvNakajimaAnn[1,24](v);", "color must be between 0 and 23"),
        ("show rvNakajimaCreate[1,0](zero[RVNakajima(B,170141183460469231731687303715884105727,0)]);", "grade exceeds signed 128-bit range"),
        ("show rvNakajimaMonomial[B]{(170141183460469231731687303715884105727,23,1)};", "grade exceeds signed 128-bit range"),
        ("show rvNakajimaScale[2](rvNakajimaScale[170141183460469231731687303715884105727](v));", "coordinate exceeds signed 128-bit range"),
        ("rvNakajimaPair(rvNakajimaMonomial[B]{(1,0,34)},rvNakajimaMonomial[B]{(1,23,34)});", "pairing exceeds signed 128-bit range"),
        ("rvNakajimaPair(rvNakajimaMonomial[B]{(1,0,1000000)},rvNakajimaMonomial[B]{(1,23,1000000)});", "Wick pairing intermediate exceeds 512 bits"),
    ] {
        let problem = diagnostic(interpret(&format!("{MARKED}{source}")));
        assert!(problem.contains(expected), "{source} {problem}");
    }
}
