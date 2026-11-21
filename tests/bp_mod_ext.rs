mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn divided_greek_classes_survive_both_first_class_connecting_maps() {
    assert_eq!(
        stdout(example("bp-divided-classes")),
        concat!(
            "BPThickExt(2, 2, 0, 2, 2, 12){g0}\n",
            "BPIExt(2, 2, 1, 1, 8){g0 + g2}\n",
            "BPExt(2, 2, 2, 8){g0 + g1}\n",
            "BPThickExt(3, 2, 0, 2, 3, 48){g0}\n",
            "BPIExt(3, 2, 1, 1, 36){2*g0 + g1}\n",
            "BPExt(3, 2, 2, 36){g0}\n",
        )
    );
}

#[test]
fn multiplication_reduction_and_connecting_are_exact_in_the_divided_examples() {
    let program = r"
        let i2 = bpIExtgen[2,2,0,1,4,0];
        let u2 = bpIExtgen[2,2,0,1,12,0];
        assert bpCup(i2,bpIExtgen[2,2,0,1,8,0]) == u2;
        assert bpThick[2](u2) == zero[BPThickExt(2,2,0,2,2,12)];
        let x2 = bpThickExtgen[2,2,0,2,2,12,0];
        let y2 = bpThickDelta(x2);
        let g0 = bpIExtgen[2,2,1,1,8,0];
        let g1 = bpIExtgen[2,2,1,1,8,1];
        let g2 = bpIExtgen[2,2,1,1,8,2];
        assert bpCup(i2,g0) == bpIExtgen[2,2,1,1,12,0];
        assert bpCup(i2,g1) == bpIExtgen[2,2,1,1,12,1];
        assert bpCup(i2,g2) == bpIExtgen[2,2,1,1,12,0];
        assert y2 == add(g0,g2);
        assert bpCup(i2,y2) == zero[BPIExt(2,2,1,1,12)];
        bpIExtdim[2,2,0,1,8];
        bpIExtdim[2,2,0,1,12];
        bpThickExtdim[2,2,0,2,2,12];
        bpIExtdim[2,2,1,1,8];
        bpIExtdim[2,2,1,1,12];
        let i3 = bpIExtgen[3,2,0,1,12,0];
        let u3 = bpIExtgen[3,2,0,1,48,0];
        assert bpCup(i3,bpIExtgen[3,2,0,1,36,0]) == u3;
        assert bpThick[3](u3) == zero[BPThickExt(3,2,0,2,3,48)];
        let x3 = bpThickExtgen[3,2,0,2,3,48,0];
        let y3 = bpThickDelta(x3);
        let h0 = bpIExtgen[3,2,1,1,36,0];
        let h1 = bpIExtgen[3,2,1,1,36,1];
        assert bpCup(i3,h0) == bpIExtgen[3,2,1,1,48,0];
        assert bpCup(i3,h1) == bpIExtgen[3,2,1,1,48,0];
        assert y3 == add(pscale[2](h0),h1);
        assert bpScale[3](bpDelta(y3)) == zero[BPExt(3,2,2,36)];
        bpIExtdim[3,2,0,1,36];
        bpIExtdim[3,2,0,1,48];
        bpThickExtdim[3,2,0,2,3,48];
        bpIExtdim[3,2,1,1,36];
        bpIExtdim[3,2,1,1,48];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "bpIExtdim[2, 2, 0, 1, 8] = 1\n",
            "bpIExtdim[2, 2, 0, 1, 12] = 1\n",
            "bpThickExtdim[2, 2, 0, 2, 2, 12] = 1\n",
            "bpIExtdim[2, 2, 1, 1, 8] = 3\n",
            "bpIExtdim[2, 2, 1, 1, 12] = 2\n",
            "bpIExtdim[3, 2, 0, 1, 36] = 1\n",
            "bpIExtdim[3, 2, 0, 1, 48] = 1\n",
            "bpThickExtdim[3, 2, 0, 2, 3, 48] = 1\n",
            "bpIExtdim[3, 2, 1, 1, 36] = 2\n",
            "bpIExtdim[3, 2, 1, 1, 48] = 1\n",
        )
    );
}

#[test]
fn quotient_products_signs_field_operations_and_representatives_are_exact() {
    let program = r"
        let t = bpReduce[1](bpt[3,2,1,1,1]);
        let x = bpIExtClass[4](t);
        assert bpCup(x,x) == zero[BPIExt(3,2,2,1,8)];
        assert bpIExtClass[4](bpIExtRep(pscale[2](x))) == pscale[2](x);
        assert add(x,pscale[2](x)) == zero[BPIExt(3,2,1,1,4)];
        let alpha = bpIExtClass[2](bpReduce[1](bpt[2,2,1,1,1]));
        let square = bpCup(alpha,alpha);
        assert square == bpIExtgen[2,2,2,1,4,0];
        let reduced = bpThick[2](square);
        assert bpCup(bpThick[2](alpha),bpThick[2](alpha)) == reduced;
        assert bpThickExtClass[4](bpThickExtRep(reduced)) == reduced;
        let v2 = bpReduce[1](bpv[3,2,0,2]);
        let thick_one = bpThickExtClass[16](bpThick[1](v2));
        let preceding = bpReduce[2](thick_one);
        assert bpDelta(preceding) == bpThickDelta(thick_one);
        assert bpDelta(bpDelta(preceding)) == bpExtClass[12](bpDelta(bpDelta(bpReduce[2](v2))));
        let raw3 = bpThick[2](bpPow[2](bpReduce[2](bpv[2,3,0,3])));
        let height_three = bpThickExtClass[28](raw3);
        let divided = bpThickDelta(height_three);
        assert divided == bpIExtClass[16](bpThickDelta(raw3));
        assert divided == add(bpIExtgen[2,3,1,2,16,0],bpIExtgen[2,3,1,2,16,1]);
        assert bpDelta(divided) == zero[BPIExt(2,3,2,1,14)];
        betti BPIExt(3,2,1,1,4);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "betti(BPIExt(3, 2, 1, 1, 4)) = [3]\n"
    );
}

#[test]
fn projection_and_connecting_respect_witness_weights_and_cocycles() {
    for (term, reason) in [
        (
            "bpIExtClass[2](bpReduce[1](bpt[3,2,1,1,1]))",
            "outside its normalized cobar basis",
        ),
        (
            "bpIExtClass[0](bpReduce[1](bpint[2,2,1,1]))",
            "outside its normalized cobar basis",
        ),
        (
            "bpIExtClass[16](bpReduce[1](bpv[3,2,0,2]))",
            "requires a closed normalized cochain",
        ),
        (
            "bpThickExtClass[16](bpThick[2](bpReduce[1](bpv[3,2,0,2])))",
            "requires a closed normalized cochain",
        ),
        ("bpIExtClass[0](zero[BPThick(2,2,0,2,2)])", "a BPI cochain"),
        (
            "bpThickExtRep(zero[BPIExt(2,2,0,1,0)])",
            "a BPThickExt class",
        ),
        (
            "bpThickDelta(zero[BPThickExt(2,2,0,2,2,2)])",
            "negative output internal degree",
        ),
        (
            "bpDelta(zero[BPIExt(3,2,0,2,0)])",
            "negative output internal degree",
        ),
        (
            "bpThickLift(zero[BPThickExt(2,2,0,2,2,0)])",
            "requires a BP cochain",
        ),
        (
            "bpThick[0](zero[BPIExt(2,2,0,1,0)])",
            "exponent must be positive",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
}
