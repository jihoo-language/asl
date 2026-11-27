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
fn vacuum_example_checks_the_central_charge_at_two_and_three() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-virasoro")),
        "rvNakajimaPair = 12\nrvNakajimaPair = 48\n"
    );
}

#[test]
fn every_color_contracts_the_integral_creation_tensor() {
    let mut source = format!("{MARKED} let delta=rvVirasoro[-2](v);");
    for c in 0..24 {
        source.push_str(&format!(
            "assert rvNakajimaAnn[1,{c}](delta)==rvNakajimaCreate[1,{c}](v);"
        ));
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn mode_moves_and_positive_contractions_count_identical_copies_once() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let y=rvNakajimaMonomial[B]{{(2,1,1),(1,0,1)}};
            assert rvVirasoro[1](y)==rvNakajimaScale[2](rvNakajimaMonomial[B]{{(1,0,1),(1,1,1)}});
            assert rvVirasoro[2](rvNakajimaMonomial[B]{{(1,1,3)}})==rvNakajimaScale[-6](rvNakajimaCreate[1,1](v));
            assert rvVirasoro[3](rvNakajimaMonomial[B]{{(1,1,1),(2,2,1)}})==rvNakajimaScale[-2](v);
            assert rvVirasoro[2](rvNakajimaMonomial[B]{{(1,0,1),(1,23,1)}})==v;
            assert rvVirasoro[0](y)==rvNakajimaScale[3](y);
            assert rvVirasoro[-1](zero[RVNakajima(B,2,3)])==zero[RVNakajima(B,3,5)];
            "
        ))),
        ""
    );
}

#[test]
fn virasoro_bracket_on_a_sparse_sum_with_even_and_odd_modes() {
    let mut source = format!(
        r"{MARKED}
        let x=add(rvNakajimaMonomial[B]{{(1,1,2),(3,2,1)}},rvNakajimaScale[-2](rvNakajimaMonomial[B]{{(1,23,1),(2,0,1),(2,1,1)}}));
        "
    );
    for (s, t) in [
        (2, -2),
        (3, -3),
        (1, -1),
        (2, 1),
        (-2, -1),
        (3, -1),
        (-3, 1),
        (0, 2),
        (0, -2),
    ] {
        let central = if s == -t { 2 * (s * s * s - s) } else { 0 };
        let rhs = if central == 0 {
            format!("rvNakajimaScale[{}](rvVirasoro[{}](x))", s - t, s + t)
        } else {
            format!(
                "add(rvNakajimaScale[{}](rvVirasoro[0](x)),rvNakajimaScale[{central}](x))",
                s - t
            )
        };
        source.push_str(&format!("assert add(rvVirasoro[{s}](rvVirasoro[{t}](x)),neg(rvVirasoro[{t}](rvVirasoro[{s}](x))))=={rhs};"));
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn heisenberg_commutators_include_crossing_zero_and_the_zero_mode() {
    let mut source = format!("{MARKED} let x=rvNakajimaMonomial[B]{{(1,1,2),(3,2,1)}};");
    let alpha = |j: i32, argument: &str| {
        if j < 0 {
            format!("rvNakajimaCreate[{},1]({argument})", -j)
        } else {
            format!("rvNakajimaAnn[{j},1]({argument})")
        }
    };
    for (s, j) in [
        (1, -3),
        (3, -2),
        (2, -2),
        (-2, -1),
        (-3, 1),
        (-1, 2),
        (-2, 2),
        (1, 1),
    ] {
        let lhs = format!(
            "add(rvVirasoro[{s}]({}),neg({}))",
            alpha(j, "x"),
            alpha(j, &format!("rvVirasoro[{s}](x)"))
        );
        let rhs = if s + j == 0 {
            "zero[RVNakajima(B,5,10)]".into()
        } else {
            format!("rvNakajimaScale[{}]({})", -j, alpha(s + j, "x"))
        };
        source.push_str(&format!("assert {lhs}=={rhs};"));
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn huge_modes_and_repeats_use_sparse_work_and_preserve_formal_zero() {
    let large = 1_000_000_000_000_000_000_000_000i128;
    let previous = large - 1;
    let far = large - 7;
    let twice = 2 * large;
    let shifted = 2 * large + 1;
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let x=rvNakajimaMonomial[B]{{(1,0,{large})}};
            assert rvVirasoro[-1](x)==rvNakajimaScale[{large}](rvNakajimaMonomial[B]{{(1,0,{previous}),(2,0,1)}});
            let y=rvNakajimaMonomial[B]{{({large},0,1)}};
            assert rvVirasoro[{previous}](y)==rvNakajimaScale[{large}](rvNakajimaCreate[1,0](v));
            assert rvVirasoro[{large}](x)==zero[RVNakajima(B,0,-{twice})];
            assert rvVirasoro[-{large}](zero[RVNakajima(B,-7,1)])==zero[RVNakajima(B,{far},{shifted})];
            "
        ))),
        ""
    );
    for (input, expected) in [
        (
            format!("show rvVirasoro[-{large}](v);"),
            "statement exceeds the computation limit",
        ),
        (
            format!("show rvVirasoro[-{}](v);", i128::MAX),
            "grade exceeds signed 128-bit range",
        ),
        (
            format!(
                "show rvVirasoro[0](rvNakajimaScale[{}](rvNakajimaCreate[2,0](v)));",
                i128::MAX
            ),
            "coordinate exceeds signed 128-bit range",
        ),
        ("show rvVirasoro[1](o);".into(), "a Nakajima state"),
    ] {
        let problem = diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
