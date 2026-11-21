mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exterior_circle_norm_counts_all_rotations_with_multiplicity() {
    assert_eq!(
        stdout(example("rank-hochconnes")),
        concat!(
            "card = 2 / 2\n",
            "RVHHChain(Exterior, 6, [3]){(0, [0]):g0[(1, [1]):g0 | (1, [1]):g0 | (1, [1]):g0]}\n",
        )
    );
    let source = include_str!("../examples/rank-hochconnes.asl");
    let mut checks = String::new();
    for n in 1..=7 {
        let coefficients = vec!["e"; n].join(",");
        let lower = 2 * n - 1;
        let upper = 2 * n;
        checks.push_str(&format!(
            "assert rvHHd(rvHHWord(e{tail})) == zero[RVHHChain(Exterior,{below},[{n}])];
             assert rvHHd(rvHHWord(one,{coefficients})) == zero[RVHHChain(Exterior,{lower},[{n}])];
             card \\x: RVHHChain(Exterior,{lower},[{n}]) => S.hi;
             card \\x: RVHHChain(Exterior,{upper},[{n}]) => S.hi;
             card \\x: RVHHChain(Exterior,{above},[{n}]) => S.hi;",
            tail = if n == 1 {
                String::new()
            } else {
                format!(",{}", vec!["e"; n - 1].join(","))
            },
            below = lower - 1,
            above = upper + 1
        ));
        let expected = if n % 2 == 0 {
            format!("zero[RVHHChain(Exterior,{upper},[{n}])]")
        } else {
            format!("rvHHWord(one,{coefficients})")
        };
        checks.push_str(&format!(
            "assert rvConnesB(rvHHWord({coefficients})) == {expected};"
        ));
    }
    let output = stdout(interpret(&format!("{source}{checks}")));
    assert_eq!(output.matches("card = 2 / 2").count(), 15);
    assert_eq!(output.matches("card = 1 / 1").count(), 7);
}

#[test]
fn the_mixed_relation_cancels_two_nonzero_composites() {
    let source = r"
        rvideal J[[0,1,2]] { rvmon[[0,1,2]]{1,0,0}, rvmon[[0,1,2]]{0,1,0}, rvmon[[0,1,2]]{0,0,1} };
        let unit = rvmtorgen[J,0,[0,0,0],0];
        let a = rvmtorgen[J,1,[1,0,0],0];
        let b = rvmtorgen[J,1,[0,1,0],0];
        let c = rvmtorgen[J,1,[0,0,1],0];
        let ab = rvTorProd(a,b);
        let bc = rvTorProd(b,c);
        let ca = rvTorProd(c,a);
        let word = rvHHWord(a,b,c);
        let expected = add(add(rvHHWord(ab,c),rvHHWord(a,bc)),rvHHWord(ca,b));
        assert rvHHd(word) == expected;
        let norm = rvConnesB(rvHHd(word));
        show norm;
        assert rvHHd(rvConnesB(word)) == norm;
        assert rvHHd(rvHHd(word)) == zero[RVHHChain(J,3,[1,1,1])];
        assert rvConnesB(rvConnesB(word)) == zero[RVHHChain(J,7,[1,1,1])];
        assert rvConnesB(word) == add(add(rvHHWord(unit,a,b,c),rvHHWord(unit,b,c,a)),rvHHWord(unit,c,a,b));
        assert rvHHWord(add(a,a),b,c) == zero[RVHHChain(J,5,[1,1,1])];
        ";
    let output = stdout(interpret(source));
    assert!(output.contains(" + "), "{output}");
}

#[test]
fn normalized_hochschild_differential_contains_the_nonformal_flag_operation() {
    let source = format!(
        r"{}
        let unit = rvmtorgen[Flag,0,[0,0,0,0,0,0],0];
        let word = rvHHWord(unit,a,b,c);
        let d = rvHHd(word);
        assert rvHHLength[1](d) == rvHHWord(unit,observed);
        assert d == add(add(rvHHWord(a,b,c),rvHHWord(c,a,b)),rvHHWord(unit,observed));
        assert rvHHd(d) == zero[RVHHChain(Flag,4,[1,1,1,1,1,1])];
        assert rvConnesB(d) == rvHHd(rvConnesB(word));
        assert rvConnesB(rvConnesB(word)) == zero[RVHHChain(Flag,8,[1,1,1,1,1,1])];
        let cyclic = rvHHWord(a,b,c);
        let ca = rvTorProd(c,a);
        let cyclic_boundary = rvHHWord(ca,b);
        assert rvHHd(cyclic) == cyclic_boundary;
        let norm = add(rvHHWord(unit,ca,b),rvHHWord(unit,b,ca));
        assert rvHHd(rvConnesB(cyclic)) == norm;
        assert rvConnesB(rvHHd(cyclic)) == norm;
        ",
        include_str!("../examples/rank-hochster.asl")
    );
    stdout(interpret(&source));
}

#[test]
fn degrees_normalization_and_ideal_agreement_are_checked() {
    let source = include_str!("../examples/rank-hochconnes.asl");
    for (expr, expected) in [
        ("rvHHWord()", "requires a coefficient"),
        ("rvHHWord(e,one)", "positive multidegree"),
        (
            "rvHHWord(e,zero[RVMTor(Exterior,0,[0])])",
            "positive multidegree",
        ),
        ("rvHHd(e)", "requires an RVHHChain"),
        ("rvHHLength[1](e)", "requires an RVHHChain"),
        ("rvHHWord(rvtorgen[Exterior,1,1,0])", "requires RVMTor"),
        (
            "rvConnesB(zero[RVHHChain(Exterior,9223372036854775807,[0])])",
            "degree is out of range",
        ),
        (
            "rvHHd(zero[RVHHChain(Exterior,-9223372036854775808,[0])])",
            "degree is out of range",
        ),
        ("zero[RVHHChain(Exterior,0,[0,0])]", "one exponent"),
    ] {
        let error = diagnostic(interpret(&format!("{source} show {expr};")));
        assert!(error.contains(expected), "{expr}: {error}");
    }
    let error = diagnostic(interpret(&format!(
        "{source} rvideal A[[0]] {{ rvmon[[0]]{{1}} }}; show rvHHWord(e,rvmtorgen[A,1,[1],0]);"
    )));
    assert!(error.contains("same rank ideal"), "{error}");
}

#[test]
fn finite_basis_and_multilinearity_include_both_generators_of_one_weight() {
    let source = r"
        rvideal J[[0,1,2]] { rvmon[[0,1,2]]{1,1,0}, rvmon[[0,1,2]]{1,0,1}, rvmon[[0,1,2]]{0,1,1} };
        let unit = rvmtorgen[J,0,[0,0,0],0];
        let a = rvmtorgen[J,2,[1,1,1],0];
        let b = rvmtorgen[J,2,[1,1,1],1];
        let A = rvHHWord(unit,a);
        let B = rvHHWord(unit,b);
        assert rvHHWord(unit,add(a,b)) == add(A,B);
        assert rvConnesB(rvHHWord(add(a,b))) == add(A,B);
        assert rvHHd(add(A,B)) == zero[RVHHChain(J,2,[1,1,1])];
        card \x: RVHHChain(J,2,[1,1,1]) => S.hi;
        card \x: RVHHChain(J,3,[1,1,1]) => S.hi;
        card \x: RVHHChain(J,4,[1,1,1]) => S.hi;
        ";
    assert_eq!(
        stdout(interpret(source)),
        "card = 4 / 4\ncard = 4 / 4\ncard = 1 / 1\n"
    );
}
