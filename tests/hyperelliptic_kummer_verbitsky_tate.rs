mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const RANK_ONE: &str = "rvverbitskytate T{3}; rvhodgewick B[21]{0,1}{[2]};";
const RANK_TWO: &str = "rvverbitskytate U{2,2}; rvhodgewick C[20]{0,1}{[0,1],[1,0]};";

#[test]
fn actual_period_maps_distinguish_augmentation_and_stable_center() {
    let output = stdout(example("hyperelliptic-kummer-verbitsky-tate"));
    let rows = output.lines().collect::<Vec<_>>();
    assert_eq!(
        rows[0],
        "rvVerbitskyTatePeriod[T, 0, ordinary, 2, 0, 0] = norm=[(0, 2, 1), (1, 1, 1), (2, 0, 1)]"
    );
    assert_eq!(rows[1], "rvVerbitskyTateCheck[T, 0, ordinary, 0] = map=augmentation rank=3 period_compositions=0 HH=3 Hochschild_cohomology=3");
    assert_eq!(rows[2], "rvVerbitskyTateCheck[T, 0, stable, 0] = map=norm rank=3 period_compositions=0 HH=2 Hochschild_cohomology=2");
    assert_eq!(
        rows[3],
        "rvVerbitskyTatePeriod[T, 0, stable, -2000000, 1, 0] = norm=[(1, 2, 1), (2, 1, 1)]"
    );
    assert_eq!(
        rows[4],
        "rvVerbitskyTateComplexMap[T, 0, stable, -2000001, cochain, 0] = [(2, 3)]"
    );
    assert_eq!(
        rows[5],
        "RVVerbitskyTate(T,stable,-2000000,3000000)[([(0, 0, -1000000)], -2999998)]"
    );
    assert_eq!(rows[6], "RVVerbitskyTate(U,ordinary,1,0)[([(0, 0, 0), (0, 1, 0)], 1), ([(0, 1, 0), (0, 0, 0)], -1)]");
    let output = stdout(interpret(&format!(
        "{RANK_ONE} {RANK_TWO}
         rvVerbitskyTatePeriod[T,0,stable,-2000001,1,0];
         rvVerbitskyTatePeriod[T,0,ordinary,0,1,1];
         rvVerbitskyTatePeriod[T,0,ordinary,0,2,1];
         rvVerbitskyTateComplexMap[T,0,ordinary,2,chain,1];
         rvVerbitskyTateComplexMap[T,0,ordinary,0,cochain,0];
         rvVerbitskyTateComplexMap[T,0,stable,-2,chain,0];
         rvVerbitskyTateCheck[T,0,ordinary,1];
         rvVerbitskyTateTensorPeriod[U]{{1,2}}{{0,0}}{{0,0}};"
    )));
    let lines = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect::<Vec<_>>();
    assert_eq!(lines, [
        "X-Y=[(1, 1, -1), (2, 0, 1)]", "augmentation=[(2, 1)]",
        "augmentation=[]", "[]", "[]", "[(2, 3)]",
        "map=X-Y rank=6 period_compositions=0 HH=2 Hochschild_cohomology=2",
        "[([0, 2], [1, 0], [0, 0], 1), ([0, 2], [0, 0], [1, 0], -1), ([1, 1], [0, 0], [0, 1], -1), ([1, 1], [0, 1], [0, 0], -1)]",
    ]);
}

#[test]
fn rank_one_bar_compares_bv_through_the_third_period_for_all_coefficients() {
    let mut source = RANK_ONE.to_string();
    for k in 0..=3 {
        for a in 0..=1 {
            let id = format!("f{k}a{a}");
            let bar = format!("bar{k}a{a}");
            let scale = 2 - a + 3 * k;
            writeln!(
                source,
                "let {id}=rvVerbitskyTateTerm[T,ordinary]{{{a}}}{{1}}{{{k}}};"
            )
            .unwrap();
            writeln!(source, "let {bar}=rvVerbitskyTateBar[B,1]({id});").unwrap();
            writeln!(source, "assert rvVerbitskyTateBV({id})==rvVerbitskyTateScale[{scale},1](rvVerbitskyTateTerm[T,ordinary]{{{a}}}{{0}}{{{k}}});").unwrap();
            writeln!(
                source,
                "assert rvVerbitskyHHBV({bar})==rvVerbitskyTateBar[B,1](rvVerbitskyTateBV({id}));"
            )
            .unwrap();
            if k <= 2 {
                writeln!(
                    source,
                    "assert rvVerbitskyHodgeProject[{}]({bar})=={bar};",
                    k + 1
                )
                .unwrap();
            }
        }
    }
    source.push_str(
        r"
        let x=rvVerbitskyTateTerm[T,ordinary]{1}{0}{0};
        let e=rvVerbitskyTateTerm[T,ordinary]{0}{1}{0};
        let t=rvVerbitskyTateTerm[T,ordinary]{0}{0}{1};
        assert rvVerbitskyTateBar[B,1](rvVerbitskyTateCup(e,t))==rvVerbitskyHHCup(rvVerbitskyTateBar[B,1](e),rvVerbitskyTateBar[B,1](t));
        assert rvVerbitskyTateBar[B,1](rvVerbitskyTateBracket(e,t))==rvVerbitskyHHBracket(rvVerbitskyTateBar[B,1](e),rvVerbitskyTateBar[B,1](t));
        assert rvVerbitskyTateBar[B,1](x)==rvVerbitskyHHBasis[B,1,0,1,0];
        assert rvVerbitskyTateTerm[T,ordinary]{2}{0}{0}==rvVerbitskyTateCup(x,x);
        assert rvVerbitskyTateTerm[T,ordinary]{2}{0}{1}==zero[RVVerbitskyTate(T,ordinary,2,-1)];
        assert rvVerbitskyTateTerm[T,ordinary]{2}{1}{0}==zero[RVVerbitskyTate(T,ordinary,1,2)];
        ",
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn localized_bv_and_gerstenhaber_signs_hold_for_signed_periods() {
    let output = stdout(interpret(&format!(
        "{RANK_ONE}
         let st=rvVerbitskyTateTerm[T,stable]{{0}}{{0}}{{1000000}};
         let si=rvVerbitskyTateTerm[T,stable]{{0}}{{0}}{{-1000000}};
         let es=rvVerbitskyTateTerm[T,stable]{{0}}{{1}}{{0}};
         assert rvVerbitskyTateCup(st,si)==rvVerbitskyTateTerm[T,stable]{{0}}{{0}}{{0}};
         assert rvVerbitskyTateTerm[T,stable]{{2}}{{0}}{{0}}==zero[RVVerbitskyTate(T,stable,0,2)];
         assert rvVerbitskyTateBracket(es,si)==rvVerbitskyTateScale[3000000,1](si);
         assert rvVerbitskyTateBV(rvVerbitskyTateCup(es,st))==rvVerbitskyTateScale[3000002,1](st);
         assert rvVerbitskyTateBV(rvVerbitskyTateCup(es,si))==rvVerbitskyTateScale[-2999998,1](si);
         assert rvVerbitskyTateBV(rvVerbitskyTateBV(rvVerbitskyTateCup(es,si)))==zero[RVVerbitskyTate(T,stable,-2000001,3000000)];
         let e=rvVerbitskyTateTerm[T,ordinary]{{0}}{{1}}{{0}};
         let xe=rvVerbitskyTateTerm[T,ordinary]{{1}}{{1}}{{0}};
         let t=rvVerbitskyTateTerm[T,ordinary]{{0}}{{0}}{{1}};
         assert rvVerbitskyTateBracket(e,xe)==xe;
         let xet=rvVerbitskyTateBracket(xe,t);
         assert xet==rvVerbitskyTateScale[-3,1](rvVerbitskyTateTerm[T,ordinary]{{1}}{{0}}{{1}});
         assert rvVerbitskyTateBracket(e,xet)==add(rvVerbitskyTateBracket(rvVerbitskyTateBracket(e,xe),t),
              rvVerbitskyTateBracket(xe,rvVerbitskyTateBracket(e,t)));
        "
    )));
    assert!(output.is_empty());
}

#[test]
fn tensor_cup_and_bv_match_the_marked_rank_two_bar_algebra() {
    let output = stdout(interpret(&format!(
        "{RANK_TWO}
         let e1=rvVerbitskyTateTerm[U,ordinary]{{0,0}}{{1,0}}{{0,0}};
         let e2=rvVerbitskyTateTerm[U,ordinary]{{0,0}}{{0,1}}{{0,0}};
         let t1=rvVerbitskyTateTerm[U,ordinary]{{0,0}}{{0,0}}{{1,0}};
         let t2=rvVerbitskyTateTerm[U,ordinary]{{0,0}}{{0,0}}{{0,1}};
         let both=rvVerbitskyTateCup(e1,e2);
         assert rvVerbitskyTateCup(e2,e1)==neg(both);
         assert rvVerbitskyTateBV(both)==add(e2,neg(e1));
         assert rvVerbitskyTateBV(rvVerbitskyTateBV(both))==zero[RVVerbitskyTate(U,ordinary,0,0)];
         assert rvVerbitskyTateCup(e1,e1)==zero[RVVerbitskyTate(U,ordinary,2,0)];
         assert rvVerbitskyTateBracket(e1,t1)==rvVerbitskyTateScale[-2,1](t1);
         assert rvVerbitskyTateBracket(e1,t2)==zero[RVVerbitskyTate(U,ordinary,2,-2)];
         let carried=add(add(rvVerbitskyHHTerm[C,1,2,-2,0]{{1,1}},rvVerbitskyHHTerm[C,1,2,-2,2]{{1,3}}),rvVerbitskyHHTerm[C,1,2,-2,2]{{3,1}});
         assert rvVerbitskyTateBar[C,1](t1)==rvVerbitskyHHClass(carried);
         let cb=rvVerbitskyTateBar[C,1](both);
         assert rvVerbitskyHHCup(rvVerbitskyTateBar[C,1](e1),rvVerbitskyTateBar[C,1](e2))==cb;
         assert rvVerbitskyHHBV(cb)==rvVerbitskyTateBar[C,1](rvVerbitskyTateBV(both));
         let s=rvVerbitskyTateCup(both,t1);
         let sb=rvVerbitskyTateBar[C,1](s);
         assert rvVerbitskyHHBV(sb)==rvVerbitskyTateBar[C,1](rvVerbitskyTateBV(s));
         assert rvVerbitskyTateHodge[3](s)==s;
         assert rvVerbitskyHodgeProject[3](sb)==sb;
        "
    )));
    assert!(output.is_empty());
}

#[test]
fn validates_modes_marks_divisor_axes_and_exponent_bounds() {
    for (query, reason) in [
        ("rvverbitskytate V{1};", "relation exponents at least two"),
        ("rvverbitskytate V{2,2}; let x=rvVerbitskyTateTerm[V,stable]{0,0}{0,0}{0,0};", "supported in rank one"),
        ("let x=rvVerbitskyTateTerm[T,ordinary]{0}{0}{-1};", "nonnegative periodic exponents"),
        ("let x=rvVerbitskyTateTerm[T,ordinary]{0}{2}{0};", "zero or one"),
        ("let x=rvVerbitskyTateTerm[T,ordinary]{0,0}{0}{0};", "one entry per factor"),
        ("let x=rvVerbitskyTateTerm[T,ordinary]{0}{0}{170141183460469231731687303715884105727};", "exceeds the signed integer range"),
        ("let x=rvVerbitskyTateBar[B,1](rvVerbitskyTateTerm[T,stable]{0}{0}{0});", "requires ordinary Tate coordinates"),
        ("let x=rvVerbitskyTateBar[B,1](rvVerbitskyTateTerm[T,ordinary]{0}{0}{1000000});", "work or storage limit"),
        ("rvverbitskytate V{2,2}; rvhodgewick D[20]{0,1}{[2,0],[0,-2]}; let x=rvVerbitskyTateBar[D,1](rvVerbitskyTateTerm[V,ordinary]{0,0}{0,0}{0,0});", "divisor axes do not present"),
        ("rvVerbitskyTatePeriod[T,0,ordinary,-1,0,0];", "nonnegative degree"),
        ("rvVerbitskyTatePeriod[T,0,ordinary,2,3,0];", "outside the truncated ring"),
        ("rvVerbitskyTatePeriod[T,1,ordinary,2,0,0];", "factor index is out of range"),
        ("let x=rvVerbitskyTateCup(rvVerbitskyTateTerm[T,ordinary]{0}{0}{0},rvVerbitskyTateTerm[T,stable]{0}{0}{0});", "same named complete intersection and mode"),
    ] {
        assert!(diagnostic(interpret(&format!("{RANK_ONE} {query}"))).contains(reason), "{query}");
    }
}
