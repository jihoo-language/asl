mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn odd_exterior_weight_has_exactly_the_two_extreme_survivors() {
    assert_eq!(
        stdout(example("rank-negative")),
        concat!(
            "RVNCChain(Exterior, 4, 6, [3]){u^0 * (0, [0]):g0[(1, [1]):g0 | (1, [1]):g0 | (1, [1]):g0]}\n",
            "RVNCChain(Exterior, 4, -1, [3]){u^3 * (1, [1]):g0[(1, [1]):g0 | (1, [1]):g0]}\n",
        )
    );
    let mut program = String::from("rvideal E[[0]] { rvmon[[0]]{1} };");
    let mut expected = String::new();
    for n in 1..=5i64 {
        for length in 1..=6i64 {
            let low = 2 * n - 1 - 2 * (length - 1);
            let high = 2 * n;
            for b in low - 1..=high + 1 {
                program.push_str(&format!("rvNCdim[E,{length},{b},[{n}]];"));
                let rank =
                    usize::from(b == low || b == high || (n % 2 == 0 && low < b && b < high));
                expected.push_str(&format!("rvNCdim[E, {length}, {b}, [{n}]] = {rank}\n"));
            }
        }
    }
    assert_eq!(stdout(interpret(&program)), expected);
}

#[test]
fn even_weight_retains_the_whole_u_module_and_maps_to_shorter_truncations() {
    let source = r"
        rvideal E[[0]] { rvmon[[0]]{1} };
        let one = rvmtorgen[E,0,[0],0];
        let e = rvmtorgen[E,1,[1],0];
        let top = rvNCClass(rvNCTerm[3,0](rvHHWord(one,e,e)));
        let t1 = rvNCU[1](top);
        let t2 = rvNCU[2](top);
        assert t1 == rvNCgen[E,3,2,[2],0];
        assert t2 == rvNCgen[E,3,0,[2],0];
        assert rvNCU[1](t1) == t2;
        assert rvNCU[1](t2) == zero[RVNC(E,3,-2,[2])];
        assert rvNCTruncate[2](t1) == rvNCgen[E,2,2,[2],0];
        assert rvNCTruncate[2](t2) == zero[RVNC(E,2,0,[2])];
        assert rvNCU[1](rvNCTruncate[2](top)) == rvNCTruncate[2](rvNCU[1](top));
        assert rvNCU[1](rvNCClass(rvNCTerm[3,0](rvHHWord(e,e)))) == rvNCgen[E,3,1,[2],0];
        assert rvNCTruncate[1](rvNCTruncate[2](top)) == rvNCTruncate[1](top);
        assert rvHHClass(rvNCPart[0](rvNCRep(rvNCTruncate[1](top)))) == rvHHgen[E,4,[2],0];
        assert rvNCU[1](add(top,top)) == zero[RVNC(E,3,2,[2])];
        card \x : RVNC(E,3,2,[2]) => S.hi;
        card \x : RVNCChain(E,3,2,[2]) => S.hi;
        ";
    assert_eq!(stdout(interpret(source)), "card = 2 / 2\ncard = 2 / 2\n");
}

#[test]
fn the_total_differential_cancels_two_nonzero_mixed_cross_terms() {
    let source = r"
        rvideal J[[0,1,2]] { rvmon[[0,1,2]]{1,0,0}, rvmon[[0,1,2]]{0,1,0}, rvmon[[0,1,2]]{0,0,1} };
        let a = rvmtorgen[J,1,[1,0,0],0];
        let b = rvmtorgen[J,1,[0,1,0],0];
        let c = rvmtorgen[J,1,[0,0,1],0];
        let x = rvHHWord(a,b,c);
        let z = rvNCTerm[3,0](x);
        let lower = rvNCTerm[3,0](rvHHd(x));
        let upper = rvNCTerm[3,1](rvConnesB(x));
        assert rvNCd(z) == add(lower,upper);
        assert rvNCd(lower) == rvNCd(upper);
        assert rvNCd(lower) == rvNCTerm[3,1](rvConnesB(rvHHd(x)));
        assert rvNCd(rvNCd(z)) == zero[RVNCChain(J,3,3,[1,1,1])];
        assert rvNCClass(rvNCd(z)) == zero[RVNC(J,3,4,[1,1,1])];
        assert rvNCd(rvNCU[1](z)) == rvNCU[1](rvNCd(z));
        assert rvNCd(rvNCTruncate[1](z)) == rvNCTruncate[1](rvNCd(z));
        show rvNCd(lower);
        ";
    let output = stdout(interpret(source));
    assert!(output.contains(" + "), "{output}");
}

#[test]
fn the_nonformal_flag_slice_uses_the_actual_hochschild_differential() {
    let source = format!(
        r"{}
        let w0 = rvNCTerm[2,0](w);
        assert rvNCd(w0) == rvNCTerm[2,0](add(left,ends));
        assert rvNCClass(rvNCd(w0)) == zero[RVNC(Flag,2,5,[1,1,1,1,1,1])];
        assert rvNCClass(rvNCTerm[2,0](left)) == rvNCClass(rvNCTerm[2,0](ends));
        rvNCdim[Flag,2,2,[1,1,1,1,1,1]];
        rvNCdim[Flag,2,3,[1,1,1,1,1,1]];
        rvNCdim[Flag,2,4,[1,1,1,1,1,1]];
        rvNCdim[Flag,2,5,[1,1,1,1,1,1]];
        rvNCdim[Flag,2,6,[1,1,1,1,1,1]];
        ",
        include_str!("../examples/rank-hochhom.asl"),
    );
    assert!(stdout(interpret(&source)).ends_with(concat!(
        "rvNCdim[Flag, 2, 2, [1, 1, 1, 1, 1, 1]] = 1\n",
        "rvNCdim[Flag, 2, 3, [1, 1, 1, 1, 1, 1]] = 2\n",
        "rvNCdim[Flag, 2, 4, [1, 1, 1, 1, 1, 1]] = 0\n",
        "rvNCdim[Flag, 2, 5, [1, 1, 1, 1, 1, 1]] = 1\n",
        "rvNCdim[Flag, 2, 6, [1, 1, 1, 1, 1, 1]] = 2\n",
    )));
}

#[test]
fn powers_and_signed_grades_are_checked_without_enumerating_a_huge_truncation() {
    let source = include_str!("../examples/rank-negative.asl");
    for (extra, expected) in [
        ("rvNCClass(start)", "requires a cycle"),
        ("rvNCTerm[4,4](low)", "below the truncation length"),
        ("rvNCPart[4](start)", "below the truncation length"),
        ("rvNCTruncate[5](top)", "cannot increase"),
        ("rvNCTruncate[0](top)", "must be positive"),
        ("zero[RVNC(Exterior,0,0,[1])]", "must be positive"),
        ("rvNCRep(start)", "requires an RVNC class"),
        (
            "rvNCd(zero[RVNCChain(Exterior,1,-9223372036854775808,[1])])",
            "degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{source} show {extra};")));
        assert!(error.contains(expected), "{extra}: {error}");
    }
    if usize::BITS == 64 {
        let program = r"
            rvideal E[[0]] { rvmon[[0]]{1} };
            let unit = rvHHWord(rvmtorgen[E,0,[0],0]);
            let x = rvNCTerm[9223372036854775808,4611686018427387904](unit);
            assert rvNCPart[4611686018427387904](x) == unit;
            assert rvNCClass(x) == rvNCgen[E,9223372036854775808,-9223372036854775808,[0],0];
            rvNCdim[E,9223372036854775808,-9223372036854775807,[0]];
            ";
        assert_eq!(
            stdout(interpret(program)),
            "rvNCdim[E, 9223372036854775808, -9223372036854775807, [0]] = 0\n",
        );
    }
}
