mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn a_nonformal_boundary_identifies_two_representatives_of_an_induced_connes_class() {
    assert_eq!(
        stdout(example("rank-hochhom")),
        concat!(
            "rvHHdim[Flag, 4, [1, 1, 1, 1, 1, 1]] = 1\n",
            "rvHHdim[Flag, 5, [1, 1, 1, 1, 1, 1]] = 3\n",
            "rvHHdim[Flag, 6, [1, 1, 1, 1, 1, 1]] = 2\n",
            "RVHH(Flag, 5, [1, 1, 1, 1, 1, 1]){g0}\n",
        )
    );
}

#[test]
fn exterior_homology_and_its_norm_have_the_expected_weight_parity() {
    let mut program = String::from(
        "rvideal E[[0]] { rvmon[[0]]{1} };
         let unit = rvmtorgen[E,0,[0],0];
         let e = rvmtorgen[E,1,[1],0];",
    );
    for n in 1..=8 {
        let factors = vec!["e"; n].join(",");
        let low = 2 * n - 1;
        let high = 2 * n;
        let expected = if n % 2 == 0 {
            format!("zero[RVHH(E,{high},[{n}])]")
        } else {
            format!("rvHHgen[E,{high},[{n}],0]")
        };
        program.push_str(&format!(
            "assert rvHHClass(rvHHWord({factors})) == rvHHgen[E,{low},[{n}],0];
             assert rvHHClass(rvHHWord(unit,{factors})) == rvHHgen[E,{high},[{n}],0];
             assert rvConnesB(rvHHgen[E,{low},[{n}],0]) == {expected};
             assert rvConnesB(rvHHgen[E,{high},[{n}],0]) == zero[RVHH(E,{above},[{n}])];",
            above = high + 1
        ));
    }
    program.push_str("assert rvHHClass(rvHHWord(unit)) == rvHHgen[E,0,[0],0]; rvHHdim[E,-1,[0]];");
    assert_eq!(stdout(interpret(&program)), "rvHHdim[E, -1, [0]] = 0\n");
}

#[test]
fn two_generators_return_independent_connes_images_and_support_addition() {
    let program = r"
        rvideal J[[0,1,2]] { rvmon[[0,1,2]]{1,1,0}, rvmon[[0,1,2]]{1,0,1}, rvmon[[0,1,2]]{0,1,1} };
        let a = rvHHgen[J,2,[1,1,1],0];
        let b = rvHHgen[J,2,[1,1,1],1];
        let A = rvHHgen[J,3,[1,1,1],0];
        let B = rvHHgen[J,3,[1,1,1],1];
        assert rvConnesB(a) == A;
        assert rvConnesB(b) == B;
        assert rvConnesB(add(a,b)) == add(A,B);
        assert rvHHClass(rvHHRep(add(A,B))) == add(A,B);
        assert rvConnesB(add(A,A)) == zero[RVHH(J,4,[1,1,1])];
        card \x: RVHH(J,2,[1,1,1]) => S.hi;
        card \x: RVHH(J,3,[1,1,1]) => S.hi;
        ";
    assert_eq!(stdout(interpret(program)), "card = 4 / 4\ncard = 4 / 4\n");
}

#[test]
fn projections_check_cycles_and_report_wrong_kinds_and_indices() {
    let source = include_str!("../examples/rank-hochhom.asl");
    for (extra, expected) in [
        ("rvHHClass(w)", "requires a cycle"),
        ("rvHHRep(w)", "requires an RVHH class"),
        ("rvHHClass(coefficient)", "requires an RVHHChain"),
        ("rvHHgen[Flag,2,[1,1,1,1,1,1],0]", "generator index"),
    ] {
        let error = diagnostic(interpret(&format!("{source} show {extra};")));
        assert!(error.contains(expected), "{extra}: {error}");
    }
}
