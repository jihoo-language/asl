mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARK: &str = "rvhodgewick B[21]{0,1}{[2]};";

#[test]
fn rank_one_full_weighted_cyclic_and_periodic_dimensions_match_oracles() {
    let mut input = MARK.to_owned();
    for degree in 0..=4 {
        for weight in 0..=8 {
            writeln!(input, "rvVerbitskyCyclicDims[B,1,hc,{degree},{weight}];").unwrap();
        }
    }
    for degree in -3..=3 {
        for weight in 0..=4 {
            writeln!(input, "rvVerbitskyCyclicDims[B,1,hp,{degree},{weight}];").unwrap();
        }
    }
    let output = stdout(interpret(&input));
    let dimensions = homology_dimensions(&output);
    let supports = dimensions[..45]
        .chunks(9)
        .map(|degrees| {
            degrees
                .iter()
                .enumerate()
                .filter(|&(_, &dim)| dim == 1)
                .map(|(i, _)| i)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert!(dimensions[..45].iter().all(|&n| n < 2));
    assert_eq!(
        supports,
        [vec![0, 1, 2], vec![], vec![0, 4, 5], vec![], vec![0, 7, 8]]
    );
    for (degree, dimensions) in (-3i32..=3).zip(dimensions[45..].chunks(5)) {
        let expected = if degree % 2 == 0 { 1 } else { 0 };
        assert_eq!(
            dimensions,
            [expected, 0, 0, 0, 0],
            "periodic degree {degree}"
        );
    }
    let far = stdout(interpret(&format!(
        "{MARK}
         let x=rvVerbitskyCyclicTerm[B,1,hp,-1000000000000,0,500000000000]{{0}};
         assert rvVerbitskyCyclicClass(x)==rvVerbitskyCyclicBasis[B,1,hp,-1000000000000,0,0];
         rvVerbitskyCyclicDims[B,1,hp,1000000000000,0];"
    )));
    assert!(far.ends_with("(C=1, Z=1, B=0, H=1)\n"));
}

fn homology_dimensions(report: &str) -> Vec<usize> {
    report
        .lines()
        .map(|row| {
            row.split_once(", H=")
                .unwrap()
                .1
                .trim_end_matches(')')
                .parse()
                .unwrap()
        })
        .collect()
}

#[test]
fn sbi_maps_and_negative_lift_match_actual_representatives() {
    let output = stdout(example("hyperelliptic-kummer-verbitsky-cyclic"));
    let rows = output.lines().collect::<Vec<_>>();
    assert_eq!(
        rows[3],
        "rvVerbitskySBI[B, 1, 2, 0] = B_into_HH=0 I=0 S=1 B=0 I_next=0 exact=true"
    );
    assert_eq!(
        rows[4],
        "rvVerbitskySBI[B, 1, 2, 1] = B_into_HH=0 I=0 S=0 B=1 I_next=0 exact=true"
    );
    assert_eq!(rows[5], "RVVerbitskyCyclic(B,1,hh,1,1)[(0, [0, 1], 1)]");
    assert_eq!(rows[6], "RVVerbitskyCyclic(B,1,hn,1,1)[(0, [0, 1], 1)]");
    assert_eq!(rows[7], "RVVerbitskyCyclicChain(B,1,hp,2,1)[(-1, [1], 1)]");
    let mut input = MARK.to_owned();
    for degree in 0..=4 {
        for weight in 0..=8 {
            writeln!(input, "rvVerbitskySBI[B,1,{degree},{weight}];").unwrap();
        }
    }
    let output = stdout(interpret(&input));
    for (index, row) in output.lines().enumerate() {
        assert!(row.ends_with("exact=true"));
        assert_eq!(row.contains("S=1"), index == 18 || index == 36, "{row}");
    }
}

#[test]
fn mixed_chain_identities_and_fujiki_adjunction_include_cyclic_signs() {
    let output = stdout(interpret(&format!(
        "{MARK}
        let x2=rvVerbitskyCyclicTerm[B,1,hh,0,2,0]{{2}};
        let bx2=rvVerbitskyChainConnes(x2);
        assert bx2==rvVerbitskyCyclicTerm[B,1,hh,1,2,0]{{0,2}};
        let h=rvVerbitskyCyclicClass(rvVerbitskyCyclicTerm[B,1,hh,1,2,0]{{1,1}});
        assert rvVerbitskyCyclicClass(bx2)==rvVerbitskyCyclicScale[2,1](h);
        assert h==rvVerbitskyCyclicBasis[B,1,hh,1,2,0];
        assert rvVerbitskyCyclicB(rvVerbitskyCyclicI(rvVerbitskyCyclicClass(x2)))==rvVerbitskyCyclicClass(bx2);
        let tower=rvVerbitskyCyclicBasis[B,1,hc,2,0,0];
        assert rvVerbitskyCyclicS(tower)==rvVerbitskyCyclicBasis[B,1,hc,0,0,0];
        assert rvVerbitskyCyclicS(rvVerbitskyCyclicBasis[B,1,hc,2,4,0])==zero[RVVerbitskyCyclic(B,1,hc,0,4)];
        let uneg=rvVerbitskyCyclicBasis[B,1,hn,0,0,0];
        assert rvVerbitskyNegativeU(uneg)==rvVerbitskyCyclicBasis[B,1,hn,-2,0,0];
        let d=add(rvVerbitskyHHTerm[B,1,1,0,1]{{1}},rvVerbitskyHHTerm[B,1,1,0,2,2]{{2}});
        rvVerbitskyHHChainPair(d,bx2);
        rvVerbitskyHHChainPair(rvVerbitskyHHBV(d),x2);
        rvVerbitskyHHChainPair(rvVerbitskyHHClass(d),rvVerbitskyCyclicClass(bx2));
        let f=rvVerbitskyHHTerm[B,1,1,-1,0]{{1}};
        let word=rvVerbitskyCyclicTerm[B,1,hh,2,3,0]{{1,1,1}};
        rvVerbitskyHHChainPair(rvVerbitskyHHDiff(f),word);
        rvVerbitskyHHChainPair(f,rvVerbitskyCyclicDiff(word));
        assert rvVerbitskyCyclicDiff(rvVerbitskyCyclicDiff(word))==zero[RVVerbitskyCyclicChain(B,1,hh,0,3)];
        let bw=rvVerbitskyCyclicDiff(rvVerbitskyChainConnes(word));
        assert add(bw,rvVerbitskyChainConnes(rvVerbitskyCyclicDiff(word)))==zero[RVVerbitskyCyclicChain(B,1,hh,2,3)];
        assert bw==rvVerbitskyCyclicScale[3,1](add(rvVerbitskyCyclicTerm[B,1,hh,2,3,0]{{0,1,2}},
               neg(rvVerbitskyCyclicTerm[B,1,hh,2,3,0]{{0,2,1}})));
        assert rvVerbitskyChainConnes(rvVerbitskyChainConnes(word))==zero[RVVerbitskyCyclicChain(B,1,hh,4,3)];
        let f2=rvVerbitskyHHTerm[B,1,2,-1,2]{{2,1}};
        let r=rvVerbitskyCyclicTerm[B,1,hh,1,3,0]{{2,1}};
        rvVerbitskyHHChainPair(f2,rvVerbitskyChainConnes(r));
        rvVerbitskyHHChainPair(rvVerbitskyHHBV(f2),r);
        "
    )));
    assert_eq!(output, "rvVerbitskyHHChainPair = 4\nrvVerbitskyHHChainPair = 4\nrvVerbitskyHHChainPair = 4\nrvVerbitskyHHChainPair = 4\nrvVerbitskyHHChainPair = 4\nrvVerbitskyHHChainPair = 2\nrvVerbitskyHHChainPair = 2\n");
}

#[test]
fn rank_two_pbw_oracle_and_a_nonzero_signed_homology_lift() {
    let mut input = "rvhodgewick C[20]{0,1}{[0,1],[1,0]};".to_owned();
    for degree in 0..=3 {
        for flavor in ["hh", "hc", "hn", "hp"] {
            for weight in 0..=4 {
                writeln!(
                    input,
                    "rvVerbitskyCyclicDims[C,1,{flavor},{degree},{weight}];"
                )
                .unwrap();
            }
        }
    }
    let output = stdout(interpret(&input));
    let dims = homology_dimensions(&output);
    assert_eq!(
        dims.chunks(5).collect::<Vec<_>>(),
        [
            [1, 2, 1, 0, 0],
            [1, 2, 1, 0, 0],
            [1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 2, 2, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 2, 1, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 1, 2, 2],
            [1, 0, 0, 2, 2],
            [0, 0, 1, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 0, 0, 2, 4],
            [0, 0, 0, 0, 2],
            [0, 0, 0, 2, 2],
            [0, 0, 0, 0, 0],
        ]
    );
    let lifted = stdout(interpret(
        r"
        rvhodgewick C[20]{0,1}{[0,1],[1,0]};
        rvVerbitskySBI[C,1,3,2];
        let c=rvVerbitskyCyclicTerm[C,1,hh,1,2,0]{1,2};
        let z=rvVerbitskyChainConnes(c);
        assert z==add(rvVerbitskyCyclicTerm[C,1,hh,2,2,0]{0,1,2},neg(rvVerbitskyCyclicTerm[C,1,hh,2,2,0]{0,2,1}));
        let hh=rvVerbitskyCyclicClass(z);
        assert hh==rvVerbitskyCyclicBasis[C,1,hh,2,2,0];
        let lift=rvVerbitskyNegativeLift(hh);
        assert rvVerbitskyNegativeEval(lift)==hh;
        let raw=rvVerbitskyNegativePeriodic(rvVerbitskyCyclicRepresentative(lift));
        let primitive=rvVerbitskyCyclicTerm[C,1,hp,3,2,-1]{1,2};
        assert rvVerbitskyCyclicDiff(primitive)==raw;
        assert rvVerbitskyNegativePeriodic(lift)==zero[RVVerbitskyCyclic(C,1,hp,2,2)];
        show lift;
        ",
    ));
    assert_eq!(lifted, "rvVerbitskySBI[C, 1, 3, 2] = B_into_HH=0 I=0 S=0 B=1 I_next=0 exact=true\nRVVerbitskyCyclic(C,1,hn,2,2)[(0, [0, 1, 2], 1), (0, [0, 2, 1], -1)]\n");
}

#[test]
fn detects_noncycles_nonlifts_bad_terms_and_work_limits() {
    for (query, reason) in [
        (
            "let x=rvVerbitskyCyclicTerm[B,1,hh,1,0,0]{0,0};",
            "bars must have positive weight",
        ),
        (
            "let x=rvVerbitskyCyclicTerm[B,1,hh,0,1,0]{};",
            "leading ring factor",
        ),
        (
            "let x=rvVerbitskyCyclicTerm[B,1,hh,0,2,0]{1};",
            "wrong internal weight",
        ),
        (
            "let x=rvVerbitskyCyclicTerm[B,1,hc,0,1,1]{1};",
            "wrong homological degree or u exponent",
        ),
        (
            "let x=rvVerbitskyCyclicClass(rvVerbitskyCyclicTerm[B,1,hp,2,1,-1]{1});",
            "not a cycle",
        ),
        (
            "let x=rvVerbitskyNegativeLift(rvVerbitskyCyclicBasis[B,1,hh,0,1,0]);",
            "no full negative cyclic lift",
        ),
        (
            "let x=rvVerbitskyCyclicPrimitive(rvVerbitskyCyclicTerm[B,1,hp,0,0,0]{0});",
            "no boundary primitive",
        ),
        (
            "let x=rvVerbitskyCyclicB(rvVerbitskyCyclicTerm[B,1,hc,0,1,0]{1});",
            "does not accept this hc chain",
        ),
        (
            "let x=rvVerbitskyCyclicBasis[B,1,hp,0,1,0];",
            "basis index out of range",
        ),
        (
            "rvVerbitskyCyclicDims[B,1,hp,0,1000000000000];",
            "work or storage limit",
        ),
        (
            "rvVerbitskyCyclicDims[B,1,unknown,0,0];",
            "expected hh, hc, hn, or hp",
        ),
    ] {
        assert!(
            diagnostic(interpret(&format!("{MARK} {query}"))).contains(reason),
            "{query}"
        );
    }
}
