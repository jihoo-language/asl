mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn odd_exterior_weight_has_a_connes_isomorphism_and_no_periodic_homology() {
    assert_eq!(
        stdout(example("rank-cyclic")),
        concat!(
            "rvHPdim[Exterior, 5, [3]] = 0\n",
            "rvHPdim[Exterior, 6, [3]] = 0\n",
            "RVHC(Exterior, 5, [3]){g0}\n",
            "RVHN(Exterior, 6, [3]){g0}\n",
        )
    );
    let mut src = String::from("rvideal J[[0]] { rvmon[[0]]{1} };");
    let mut expected = String::new();
    for n in 0..=6i64 {
        for b in -6..=2 * n + 6 {
            let (hc, hn, hp) = if n == 0 {
                (b % 2 == 0 && b >= 0, b % 2 == 0 && b <= 0, b % 2 == 0)
            } else if n % 2 == 0 {
                (b >= 2 * n - 1, b <= 2 * n, true)
            } else {
                (b == 2 * n - 1, b == 2 * n, false)
            };
            for (mode, rank) in [("HC", hc), ("HN", hn), ("HP", hp)] {
                src.push_str(&format!("rv{mode}dim[J,{b},[{n}]];"));
                expected.push_str(&format!(
                    "rv{mode}dim[J, {b}, [{n}]] = {}\n",
                    usize::from(rank)
                ));
            }
        }
    }
    assert_eq!(stdout(interpret(&src)), expected);
}

#[test]
fn periodic_u_is_invertible_and_ordinary_s_forgets_the_top_coefficient() {
    let source = r"
        rvideal E[[0]] { rvmon[[0]]{1} };
        let e = rvmtorgen[E,1,[1],0];
        let x = rvHHWord(e,e);
        let periodic = rvCyclicClass(rvHPTerm[0](x));
        assert periodic == rvHPgen[E,3,[2],0];
        assert rvCyclicU[-1](rvCyclicU[1](periodic)) == periodic;
        assert rvCyclicU[2](rvCyclicU[-3](periodic)) == rvCyclicU[-1](periodic);
        assert rvCyclicD(rvCyclicU[-1](rvHPTerm[0](x))) == rvCyclicU[-1](rvCyclicD(rvHPTerm[0](x)));
        assert rvCyclicPart[-1](rvCyclicRep(rvCyclicU[-1](periodic))) == x;
        let h = rvCyclicClass(rvHCTerm[-1](x));
        assert h == rvHCgen[E,5,[2],0];
        assert rvCyclicS(h) == rvCyclicI(rvHHClass(x));
        assert rvCyclicS(rvCyclicS(h)) == zero[RVHC(E,1,[2])];
        assert rvCyclicB(h) == zero[RVHH(E,6,[2])];
        assert rvCyclicI(add(rvHHClass(x),rvHHClass(x))) == zero[RVHC(E,3,[2])];
        card \x: RVHP(E,3,[2]) => S.hi;
        card \x: RVHPChain(E,3,[2]) => S.hi;
        ";
    assert_eq!(stdout(interpret(source)), "card = 2 / 2\ncard = 2 / 2\n");
}

#[test]
fn flag_sbi_uses_the_actual_nonformal_homology_and_boundary() {
    let program = format!(
        r"{}
        let ordinary = rvCyclicI(coefficient);
        assert ordinary == rvHCgen[Flag,4,[1,1,1,1,1,1],0];
        assert rvCyclicB(ordinary) == rvConnesB(coefficient);
        assert rvCyclicI(rvCyclicB(ordinary)) == zero[RVHC(Flag,5,[1,1,1,1,1,1])];
        let coeff = rvHHWord(top);
        assert rvCyclicD(rvHCTerm[-1](coeff)) == rvHCTerm[0](left);
        assert rvCyclicB(rvCyclicI(coeff)) == left;
        assert rvCyclicD(rvCyclicI(w)) == rvCyclicI(rvHHd(w));
        assert rvCyclicClass(rvCyclicI(left)) == rvCyclicClass(rvCyclicI(ends));
        assert rvCyclicS(ordinary) == zero[RVHC(Flag,2,[1,1,1,1,1,1])];
        rvHCdim[Flag,5,[1,1,1,1,1,1]];
        rvHNdim[Flag,5,[1,1,1,1,1,1]];
        rvHNdim[Flag,6,[1,1,1,1,1,1]];
        rvHPdim[Flag,6,[1,1,1,1,1,1]];
        ",
        include_str!("../examples/rank-hochhom.asl"),
    );
    assert!(stdout(interpret(&program)).ends_with(concat!(
        "rvHCdim[Flag, 5, [1, 1, 1, 1, 1, 1]] = 2\n",
        "rvHNdim[Flag, 5, [1, 1, 1, 1, 1, 1]] = 1\n",
        "rvHNdim[Flag, 6, [1, 1, 1, 1, 1, 1]] = 2\n",
        "rvHPdim[Flag, 6, [1, 1, 1, 1, 1, 1]] = 0\n",
    )));
}

#[test]
fn signed_degree_windows_and_power_diagnostics_are_exact() {
    let program = r"
        rvideal E[[0]] { rvmon[[0]]{1} };
        let unit = rvHHWord(rvmtorgen[E,0,[0],0]);
        let lower = rvHPTerm[4611686018427387904](unit);
        let upper = rvHCTerm[-4611686018427387903](unit);
        assert rvCyclicPart[4611686018427387904](lower) == unit;
        assert rvCyclicClass(lower) == rvHPgen[E,-9223372036854775808,[0],0];
        assert rvCyclicClass(upper) == rvHCgen[E,9223372036854775806,[0],0];
        assert rvCyclicU[-4611686018427387904](rvCyclicClass(lower)) == rvHPgen[E,0,[0],0];
        assert rvCyclicS(rvCyclicClass(upper)) == rvHCgen[E,9223372036854775804,[0],0];
        rvHPdim[E,9223372036854775807,[0]];
        rvHNdim[E,-9223372036854775808,[0]];
        ";
    assert_eq!(
        stdout(interpret(program)),
        "rvHPdim[E, 9223372036854775807, [0]] = 0\nrvHNdim[E, -9223372036854775808, [0]] = 1\n"
    );
    let source = include_str!("../examples/rank-cyclic.asl");
    for (expr, expected) in [
        ("rvHNTerm[-1](low)", "nonnegative"),
        ("rvHCTerm[1](low)", "nonpositive"),
        ("rvCyclicU[1](hcLow)", "requires a periodic"),
        ("rvCyclicS(hnHigh)", "requires an ordinary"),
        ("rvCyclicB(hnHigh)", "requires an ordinary"),
        ("rvCyclicClass(rvHPTerm[0](low))", "requires a cycle"),
        (
            "rvCyclicD(zero[RVHPChain(Exterior,-9223372036854775808,[0])])",
            "degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{source} show {expr};")));
        assert!(error.contains(expected), "{expr}: {error}");
    }
}
