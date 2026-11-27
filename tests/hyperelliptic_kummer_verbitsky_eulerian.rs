mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const B: &str = "rvhodgewick B[21]{0,1}{[2]};";
const C: &str = "rvhodgewick C[20]{0,1}{[0,1],[1,0]};";

#[test]
fn exact_adams_identities_and_transpose_use_the_signed_chain_action() {
    let output = stdout(interpret(&format!(
        "{B} {C}
        let c=rvVerbitskyCyclicTerm[B,1,hh,2,4,0]{{1,1,2}};
        assert rvVerbitskyAdams[2](rvVerbitskyAdams[3](c))==rvVerbitskyAdams[6](c);
        assert rvVerbitskyAdams[3](rvVerbitskyCyclicDiff(c))==rvVerbitskyCyclicDiff(rvVerbitskyAdams[3](c));
        assert rvVerbitskyAdams[3](rvVerbitskyChainConnes(c))==rvVerbitskyCyclicScale[3,1](rvVerbitskyChainConnes(rvVerbitskyAdams[3](c)));
        assert rvVerbitskyAdams[1](c)==c;
        assert rvVerbitskyAdams[0](c)==zero[RVVerbitskyCyclicChain(B,1,hh,2,4)];
        let x=rvVerbitskyCyclicTerm[B,1,hh,0,1,0]{{1}};
        assert rvVerbitskyAdams[0](x)==x;
        let f=rvVerbitskyHHTerm[B,1,2,-2,1]{{2,1}};
        rvVerbitskyHHChainPair(rvVerbitskyAdams[2](f),c);
        rvVerbitskyHHChainPair(f,rvVerbitskyAdams[2](c));
        assert rvVerbitskyHHDiff(rvVerbitskyAdams[2](f))==rvVerbitskyAdams[2](rvVerbitskyHHDiff(f));
        let t=rvVerbitskyCyclicTerm[C,1,hh,3,4,0]{{0,1,2,3}};
        let p1=rvVerbitskyHodgeProject[1](t);
        let p2=rvVerbitskyHodgeProject[2](t);
        let p3=rvVerbitskyHodgeProject[3](t);
        assert add(add(p1,p2),p3)==t;
        assert rvVerbitskyHodgeProject[2](p2)==p2;
        assert rvVerbitskyHodgeProject[1](p3)==zero[RVVerbitskyCyclicChain(C,1,hh,3,4)];
        assert rvVerbitskyAdams[3](t)==add(rvVerbitskyCyclicScale[3,1](p1),add(rvVerbitskyCyclicScale[9,1](p2),rvVerbitskyCyclicScale[27,1](p3)));
        "
    )));
    assert_eq!(
        output,
        "rvVerbitskyHHChainPair = -2\nrvVerbitskyHHChainPair = -2\n"
    );
}

#[test]
fn hodge_columns_preserve_mixed_degree_and_shift_under_sbi() {
    let demo = stdout(example("hyperelliptic-kummer-verbitsky-eulerian"));
    let rows = demo.lines().collect::<Vec<_>>();
    assert_eq!(rows[2], "rvVerbitskyHodgeDim[C, 1, hh, 4, 6, 2] = 3");
    assert_eq!(rows[3], "rvMonomialCICotangentDim[4, 2, 6] = 3");
    assert_eq!(rows[5], "RVVerbitskyHH(C,1,2,0)[([2, 1], 3, 1)]");
    assert_eq!(rows[6], "RVVerbitskyCyclic(B,1,hh,1,1)[(0, [0, 1], 1)]");
    let output = stdout(interpret(&format!(
        "{B}
        let c=rvVerbitskyCyclicTerm[B,1,hc,3,3,-1]{{2,1}};
        let dc=rvVerbitskyCyclicDiff(c);
        assert c==rvVerbitskyHodgeProject[2](c);
        assert dc==rvVerbitskyHodgeProject[2](dc);
        assert dc==add(rvVerbitskyCyclicTerm[B,1,hc,2,3,0]{{0,2,1}},neg(rvVerbitskyCyclicTerm[B,1,hc,2,3,0]{{0,1,2}}));
        assert rvVerbitskyHodgeProject[1](dc)==zero[RVVerbitskyCyclicChain(B,1,hc,2,3)];
        let uneg=rvVerbitskyCyclicBasis[B,1,hn,0,0,0];
        let u=rvVerbitskyNegativeU(uneg);
        assert rvVerbitskyHodgeProject[-1](u)==u;
        assert rvVerbitskyHodgeProject[0](u)==zero[RVVerbitskyCyclic(B,1,hn,-2,0)];
        rvVerbitskyHodgeDim[B,1,hc,2,4,1];
        rvVerbitskyHodgeDim[B,1,hc,2,4,2];
        rvVerbitskyHodgeDim[B,1,hc,4,7,2];
        rvVerbitskyHodgeDim[B,1,hn,3,4,2];
        rvVerbitskyHodgeDim[B,1,hn,3,4,1];
        let p=rvVerbitskyCyclicBasis[B,1,hp,-1000000000000,0,0];
        assert rvVerbitskyHodgeProject[-500000000000](p)==p;
        assert rvVerbitskyHodgeProject[0](p)==zero[RVVerbitskyCyclic(B,1,hp,-1000000000000,0)];
        rvVerbitskyHodgeDim[B,1,hp,-1000000000000,0,-500000000000];
        rvVerbitskyHodgeDim[B,1,hp,1000000000000,0,500000000000];
        "
    )));
    assert_eq!(dimensions(&output), [1, 0, 1, 1, 0, 1, 1]);
}

#[test]
fn harrison_carries_and_quantum_alternating_class_are_separate() {
    let output = stdout(interpret(&format!(
        "{C}
        let tu=add(add(rvVerbitskyHHTerm[C,1,2,-2,0]{{1,1}},rvVerbitskyHHTerm[C,1,2,-2,2]{{1,3}}),rvVerbitskyHHTerm[C,1,2,-2,2]{{3,1}});
        let tv=add(add(rvVerbitskyHHTerm[C,1,2,-2,0]{{2,2}},rvVerbitskyHHTerm[C,1,2,-2,1]{{2,3}}),rvVerbitskyHHTerm[C,1,2,-2,1]{{3,2}});
        let u=rvVerbitskyHHClass(tu);
        let v=rvVerbitskyHHClass(tv);
        assert rvVerbitskyHodgeProject[1](u)==u;
        assert rvVerbitskyHodgeProject[1](v)==v;
        assert rvVerbitskyHodgeProject[2](u)==zero[RVVerbitskyHH(C,1,2,-2)];
        assert rvVerbitskyHodgeProject[2](v)==zero[RVVerbitskyHH(C,1,2,-2)];
        rvVerbitskyHHCoeff[0]{{1,1}}(u);
        rvVerbitskyHHCoeff[0]{{2,2}}(v);
        let p=rvVerbitskyHHTerm[C,1,2,0,3]{{2,1}};
        let other=rvVerbitskyHHTerm[C,1,2,0,3]{{1,2}};
        let alternating=rvVerbitskyHHClass(rvVerbitskyHHScale[1,2](add(p,neg(other))));
        let projected=rvVerbitskyHodgeProject[2](rvVerbitskyHHClass(p));
        assert projected==alternating;
        assert rvVerbitskyHodgeProject[1](rvVerbitskyHHClass(p))==zero[RVVerbitskyHH(C,1,2,0)];
        assert rvVerbitskyAdams[3](projected)==rvVerbitskyHHScale[9,1](projected);
        rvVerbitskyHHCoeff[3]{{2,1}}(projected);
        "
    )));
    assert_eq!(dimensions(&output), [1, 1, 1]);
}

fn dimensions(report: &str) -> Vec<usize> {
    report
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn eulerian_dimensions_match_an_independent_derived_exterior_complex() {
    let mut source = format!("{B} {C}");
    let mut cases = Vec::new();
    for (mark, variables, degrees) in [("B", "[3,1]", 4), ("C", "[2,1],[2,1]", 3)] {
        for d in 0..=degrees {
            for index in 0..=d {
                for w in 0..=2 * d + 2 {
                    cases.push((mark, variables, d, index, w));
                }
            }
        }
    }
    for (i, w) in [(2, 5), (2, 6), (3, 4), (2, 4), (3, 6), (4, 4)] {
        cases.push(("C", "[2,1],[2,1]", 4, i, w));
    }
    for (mark, variables, d, i, w) in cases {
        writeln!(source, "rvVerbitskyHodgeDim[{mark},1,hh,{d},{w},{i}];").unwrap();
        writeln!(
            source,
            "rvMonomialCICotangentDim[{d},{i},{w}]{{{variables}}};"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let rows = output.lines().collect::<Vec<_>>();
    for pair in rows.chunks(2) {
        let left = pair[0].split_once(" = ").unwrap().1;
        let right = pair[1].split_once(" = ").unwrap().1;
        assert_eq!(left, right, "{pair:?}");
    }
    assert_eq!(
        dimensions(&rows[rows.len() - 12..].join("\n")),
        [2, 2, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0]
    );
}

#[test]
fn divided_power_differentials_preserve_weights_and_use_wedge_signs() {
    let output = stdout(interpret(
        r"
        rvMonomialCICotangentDiff{[2,1],[2,1]}{0,0}{2,1}{};
        rvMonomialCICotangentDiff{[2,1],[2,1]}{0,0}{0,1}{0};
        rvMonomialCICotangentDiff{[3,2]}{1}{5}{};
        rvMonomialCICotangentDim[4,2,14]{[3,2]};
        rvMonomialCICotangentDim[4,2,12]{[3,2]};
        rvMonomialCICotangentDim[3,2,8]{[3,2]};
        ",
    ));
    let rows = output.lines().collect::<Vec<_>>();
    assert_eq!(
        rows[0],
        "rvMonomialCICotangentDiff = [(2, [1, 0], [1, 1], [0]), (2, [0, 1], [2, 0], [1])]"
    );
    assert_eq!(
        rows[1],
        "rvMonomialCICotangentDiff = [(-2, [0, 1], [0, 0], [0, 1])]"
    );
    assert_eq!(rows[2], "rvMonomialCICotangentDiff = []");
    assert_eq!(dimensions(&rows[3..].join("\n")), [1, 0, 1]);
}

#[test]
fn uniform_words_use_exact_powers_and_general_shuffles_have_checked_work() {
    let word = vec!["1"; 64].join(",");
    let source = format!(
        "{B}
         let z=rvVerbitskyCyclicTerm[B,1,hh,64,64,0]{{0,{word}}};
         assert rvVerbitskyAdams[2](z)==rvVerbitskyCyclicScale[4294967296,1](z);
         assert rvVerbitskyHodgeProject[32](z)==z;
         assert rvVerbitskyHodgeProject[31](z)==zero[RVVerbitskyCyclicChain(B,1,hh,64,64)];"
    );
    assert!(stdout(interpret(&source)).is_empty());
    let word = vec!["1"; 31].join(",");
    for (query, reason) in [
        (format!("let c=rvVerbitskyCyclicTerm[B,1,hh,32,33,0]{{0,{word},2}}; show rvVerbitskyAdams[2](c);"), "computation limit"),
        (String::from("let d=rvVerbitskyAdams[-1](rvVerbitskyCyclicBasis[B,1,hh,0,0,0]);"), "nonnegative"),
        (String::from("let d=rvVerbitskyAdams[2](rvVerbitskyCyclicBasis[B,1,hp,0,0,0]);"), "Hochschild chains or cochains"),
        (String::from("rvMonomialCICotangentDim[0,0,0]{[0,1]};"), "positive relation powers"),
        (String::from("rvMonomialCICotangentDiff{[2,1]}{2}{1}{};"), "exceeds a relation power"),
        (String::from("rvMonomialCICotangentDiff{[2,1]}{0}{1}{0,0};"), "distinct, in range, and increasing"),
    ] {
        assert!(diagnostic(interpret(&format!("{B} {query}"))).contains(reason), "{query}");
    }
}
