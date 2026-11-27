mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[W]{};
    let o=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
    let e1=rvJacPoint[C](z,z);
    let f1=rvJacPoint[C](o,z);
    let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
    let f2=add(add(e1,f1),rvJacPoint[C](i,z));
    rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};
";

#[test]
fn weight_three_three_cycle_excess_is_the_signed_euler_contribution() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-lehn-sorger")),
        "rvNakajimaPair = -216\n"
    );
}

fn hilbert_square_axes() -> Vec<(usize, String)> {
    let degree = |color| match color {
        0 => 0,
        23 => 4,
        _ => 2,
    };
    let mut axes = Vec::new();
    for i in 0..24 {
        for j in i..24 {
            let factors = if i == j {
                format!("(1,{i},2)")
            } else {
                format!("(1,{i},1),(1,{j},1)")
            };
            axes.push((degree(i) + degree(j), factors));
        }
        axes.push((degree(i) + 2, format!("(2,{i},1)")));
    }
    axes
}

#[test]
fn weight_two_cup_matches_geometric_axes_in_each_degree_block() {
    let axes = hilbert_square_axes();
    assert_eq!(axes.len(), 324);
    let mut comparisons = Vec::new();
    for (i, (x_degree, x)) in axes.iter().enumerate() {
        for (j, (y_degree, y)) in axes.iter().enumerate().skip(i) {
            let paired_representative = x == "(1,1,2)"
                || y == "(1,1,2)"
                || x == "(2,1,1)"
                || y == "(2,1,1)"
                || x == "(1,0,1),(1,23,1)"
                || y == "(1,0,1),(1,23,1)";
            if x_degree + y_degree <= 8
                && (x_degree == &0
                    || y_degree == &0
                    || (x_degree == &2 && y_degree == &2)
                    || x == "(2,0,1)"
                    || y == "(2,0,1)"
                    || paired_representative
                    || i == j)
            {
                comparisons.push((x, y));
            }
        }
    }
    for chunk in comparisons.chunks(45) {
        let mut source = MARKED.to_owned();
        for (j, (left, right)) in chunk.iter().enumerate() {
            writeln!(source, r"
                let a{j}=rvNakajimaMonomial[B]{{{left}}};
                let b{j}=rvNakajimaMonomial[B]{{{right}}};
                assert rvNakajimaHilbSquare(rvNakajimaCup(a{j},b{j}))==rvHilbSquareCup(rvNakajimaHilbSquare(a{j}),rvNakajimaHilbSquare(b{j}));
            ").unwrap();
        }
        assert_eq!(stdout(interpret(&source)), "");
    }
}

#[test]
fn scaled_units_and_half_boundaries_match_on_several_partitions() {
    for (n, n_factorial, boundary_factor, axes) in [
        (
            2,
            2,
            2,
            &["(1,0,2)", "(1,1,1),(1,23,1)", "(2,0,1)", "(2,1,1)"][..],
        ),
        (
            3,
            6,
            2,
            &[
                "(1,0,1),(1,1,2)",
                "(1,0,1),(2,1,1)",
                "(1,1,1),(2,23,1)",
                "(3,0,1)",
                "(3,23,1)",
            ],
        ),
        (
            4,
            24,
            4,
            &[
                "(1,0,2),(1,1,2)",
                "(1,1,2),(2,0,1)",
                "(2,0,1),(2,1,1)",
                "(1,0,1),(3,2,1)",
                "(4,0,1)",
                "(4,23,1)",
            ],
        ),
    ] {
        let mut source = format!(
            r"{MARKED}
            let u=rvNakajimaMonomial[B]{{(1,0,{n})}};
            let d=rvNakajimaMonomial[B]{{(1,0,{}),(2,0,1)}};
        ",
            n - 2
        );
        for (index, factors) in axes.iter().enumerate() {
            writeln!(source, r"
                let x{index}=rvNakajimaMonomial[B]{{{factors}}};
                assert rvNakajimaCup(u,x{index})==rvNakajimaScale[{n_factorial}](x{index});
                assert rvNakajimaCup(x{index},u)==rvNakajimaScale[{n_factorial}](x{index});
                assert rvNakajimaCup(d,x{index})==rvNakajimaScale[{boundary_factor}](rvNakajimaBoundary(x{index}));
                assert rvNakajimaCup(x{index},d)==rvNakajimaCup(d,x{index});
            ").unwrap();
        }
        assert_eq!(stdout(interpret(&source)), "");
    }
}

#[test]
fn associativity_commutativity_and_signed_wick_frobenius_hold_for_mixed_terms() {
    let source = format!(
        r"{MARKED}
        let d=rvNakajimaMonomial[B]{{(1,0,1),(2,0,1)}};
        let e=rvNakajimaMonomial[B]{{(1,0,2),(1,1,1)}};
        let p=rvNakajimaMonomial[B]{{(1,0,2),(1,23,1)}};
        let a=rvNakajimaMonomial[B]{{(3,0,1)}};
        let x=add(d,rvNakajimaScale[-2](e));
        let y=add(a,rvNakajimaScale[3](p));
        let xy=rvNakajimaCup(x,y);
        assert xy==rvNakajimaCup(y,x);
        assert rvNakajimaCup(xy,x)==rvNakajimaCup(x,rvNakajimaCup(y,x));
        assert rvNakajimaCup(x,y)==add(rvNakajimaCup(d,y),rvNakajimaScale[-2](rvNakajimaCup(e,y)));
        let z3=rvNakajimaMonomial[B]{{(3,1,1)}};
        rvNakajimaPair(xy,z3);
        rvNakajimaPair(x,rvNakajimaCup(y,z3));
        let b=rvNakajimaMonomial[B]{{(2,0,2)}};
        let c=rvNakajimaMonomial[B]{{(1,0,1),(3,0,1)}};
        let f=rvNakajimaMonomial[B]{{(1,0,1),(1,1,1),(2,0,1)}};
        assert rvNakajimaCup(b,c)==rvNakajimaCup(c,b);
        assert rvNakajimaCup(rvNakajimaCup(b,c),f)==rvNakajimaCup(b,rvNakajimaCup(c,f));
    "
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvNakajimaPair = 72\nrvNakajimaPair = 72\n"
    );
}

#[test]
fn formal_zeros_and_wide_cancellation_do_not_expand_spurious_cycles() {
    let max = i128::MAX;
    let n = 1_000_000_000_000_000_000_000_000i128;
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let v=rvNakajimaVac[B];
        assert rvNakajimaCup(v,v)==v;
        let h=rvNakajimaMonomial[B]{{(1,1,1)}};
        assert rvNakajimaCup(h,h)==rvNakajimaScale[-2](rvNakajimaMonomial[B]{{(1,23,1)}});
        let a=rvNakajimaMonomial[B]{{(1,1,2)}};
        let b=rvNakajimaMonomial[B]{{(1,2,2)}};
        let c=rvNakajimaMonomial[B]{{(1,1,1),(1,2,1)}};
        let wide=add(rvNakajimaScale[{max}](a),rvNakajimaScale[-{max}](b));
        assert rvNakajimaCup(wide,c)==zero[RVNakajima(B,2,8)];
        let minimum=rvNakajimaScale[{}](rvNakajimaMonomial[B]{{(1,0,2)}});
        let minus=neg(rvNakajimaMonomial[B]{{(1,23,2)}});
        assert rvNakajimaCup(minimum,minus)==rvNakajimaScale[{}](rvNakajimaMonomial[B]{{(1,23,2)}});
        let long=rvNakajimaMonomial[B]{{(1,0,{n})}};
        let top=rvNakajimaMonomial[B]{{(1,23,{n})}};
        assert rvNakajimaCup(long,zero[RVNakajima(B,{n},2)])==zero[RVNakajima(B,{n},2)];
        assert rvNakajimaCup(top,top)==zero[RVNakajima(B,{n},{})];
        assert rvNakajimaCup(zero[RVNakajima(B,-8,3)],zero[RVNakajima(B,-8,-3)])==zero[RVNakajima(B,-8,0)];
    ",
            1i128 << 126,
            i128::MIN,
            8 * n
        ))),
        ""
    );
}

#[test]
fn work_limits_and_incompatible_grades_have_diagnostics() {
    for (input, expected) in [
        ("let x=rvNakajimaMonomial[B]{(1,0,9),(1,1,1)}; show rvNakajimaCup(x,x);".to_owned(), "conjugator work limit"),
        ("let x=rvNakajimaMonomial[B]{(2,0,4)}; show rvNakajimaCup(x,x);".to_owned(), "diagonal expansion work limit"),
        ("show rvNakajimaCup(rvNakajimaVac[B],rvNakajimaMonomial[B]{(1,0,1)});".to_owned(), "same named level basis and length"),
        ("rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; show rvNakajimaCup(rvNakajimaVac[B],rvNakajimaVac[D]);".to_owned(), "same named level basis and length"),
        ("show rvNakajimaCup(o,rvNakajimaVac[B]);".to_owned(), "a Nakajima state"),
        (format!("show rvNakajimaCup(zero[RVNakajima(B,1,{})],zero[RVNakajima(B,1,2)]);", i128::MAX), "grade exceeds signed 128-bit range"),
        (format!("let x=rvNakajimaScale[{}](rvNakajimaMonomial[B]{{(1,0,2)}}); show rvNakajimaCup(x,x);", i128::MAX), "coordinate exceeds signed 128-bit range"),
    ] {
        let err = diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(err.contains(expected), "{err}");
    }
}
