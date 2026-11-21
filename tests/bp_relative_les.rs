mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn height_les_detects_the_weight_six_and_eight_extensions() {
    assert_eq!(
        stdout(example("bp-height-les")),
        concat!(
            "bporder[6] = 2\n",
            "BPExt(2, 1, 2, 6){g0}\n",
            "bporder[8] = 2\n",
            "BPRelExt(2, 1, 2, 1, 8){g0}\n",
        )
    );
}

#[test]
fn relative_classes_are_first_class_and_connecting_is_linear() {
    let source = r"
        let x = bprelgen[2,1,2,1,6,0];
        let d = bpHeightDelta(x);
        let z = zero[BPRelExt(2,1,2,1,6)];
        assert bpRelRep(z) == zero[BPRel(2,1,2,1)];
        assert bpRelClass[6](bpRelRep(z)) == z;
        assert bpHeightDelta(bpScale[-1](x)) == bpNeg(d);
        assert bpHeightDelta(add(x,x)) == add(d,d);
        assert bpD(bpNull[6](bpRelRep(x))) == bpScale[2](bpRelRep(x));
        assert bpRelClass[6](bpScale[2](bpRelRep(x))) == z;
        let lift = \q : BPRelExt(2,1,2,1,6) => bpHeightDelta(q);
        assert lift(x) == d;
        card(\q : BPRelExt(2,1,2,1,6) => S.hi);

        let y = bprelgen[3,1,2,1,16,0];
        assert bpRelClass[16](bpRelRep(y)) == y;
        assert bpNeg(y) == bpScale[2](y);
        assert add(y,bpNeg(y)) == zero[BPRelExt(3,1,2,1,16)];
        bporder[16](bpRelRep(y));
        assert bpD(bpNull[16](bpRelRep(y))) == bpScale[3](bpRelRep(y));
        card(\q : BPRelExt(3,1,2,1,16) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(source)),
        "card = 2 / 2\nbporder[16] = 3\ncard = 3 / 3\n"
    );
}

#[test]
fn height_zero_quotient_preserves_positive_degree_classes_only() {
    let source = r"
        let unit = bpextgen[2,0,0,0,0];
        assert bpRaise[2](unit) == bpextgen[2,2,0,0,0];
        assert bpRel[0](bpRaise[2](unit)) == zero[BPRelExt(2,0,2,0,0)];
        let t = bpextgen[2,2,1,2,0];
        let r = bpRel[0](t);
        assert r == bprelgen[2,0,2,1,2,0];
        assert bpHeightDelta(r) == zero[BPExt(2,0,2,2)];
        assert bpHeightDelta(bpRelRep(r)) == zero[BP(2,0,2)];
        bporder[2](bpRelRep(r));
    ";
    assert_eq!(stdout(interpret(source)), "bporder[2] = 2\n");
}

#[test]
fn the_height_connecting_map_checks_cycle_and_types() {
    for (source, reason) in [
        (
            "show bpHeightDelta(bpRel[1](bpv[2,2,0,2]));",
            "bpHeightDelta requires a closed relative cochain",
        ),
        (
            "show bpRelClass[6](bpRel[1](bpv[2,2,0,2]));",
            "BP Ext requires a closed cochain",
        ),
        (
            "show bprelgen[2,1,2,2,6,0];",
            "bprelgen index 0 is out of range for rank 0",
        ),
        (
            "show bpRelRep(bpextgen[2,2,1,8,0]);",
            "bpRelRep requires a BPRelExt class",
        ),
        (
            "show bpRel[2](bpextgen[2,2,1,8,0]);",
            "relative heights must satisfy n < m",
        ),
        (
            "let x=bprelgen[2,1,2,1,6,0]; show bpCup(x,x);",
            "bpCup requires a BP polynomial",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
