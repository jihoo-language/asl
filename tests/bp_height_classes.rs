mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn typed_height_pages_lift_the_extension_and_correct_a_representative() {
    assert_eq!(
        stdout(example("bp-height-classes")),
        concat!(
            "BPHeight(2, 2, 1, 1, 2, 6){g0}\n",
            "BP(2, 2, 2){v1*t1_1*t2_1^2 - 3*v1*t1_1^2*t2_1 - 2*v1^2*t1_1*t2_1 - 2*t1_1*t2_1^3 - 2*t1_1^2*t2_1^2 - 4*t1_1^3*t2_1}\n",
            "BPHeight(2, 2, 3, 1, 1, 8){4*g0}\n",
            "BPHeight(2, 2, 3, 1, 1, 8){g0}\n",
            "BPHeight(2, 2, 3, 2, 1, 8){g0}\n",
            "bporder[8] = 8\n",
            "bporder[8] = 16\n",
        )
    );
}

#[test]
fn height_page_zero_has_integral_coordinates_and_differential() {
    let source = r"
        let v2 = bpv[2,2,0,2];
        let t2 = bpt[2,2,1,1,2];
        let x = bpHeightgen[2,2,0,2,0,6,0];
        let y = bpHeightgen[2,2,0,2,1,6,0];
        assert bpHeightClass[0,2,6](v2) == x;
        assert bpHeightClass[0,2,6](t2) == y;
        assert bpHeightD(x) == bpScale[2](y);
        assert bpHeightD(bpHeightD(x)) == zero[BPHeight(2,2,0,2,2,6)];
        assert bpHeightRep(bpNeg(bpScale[5](y))) == bpMul(bpint[2,2,1,-5],t2);
        assert bpHeightNext(y) == bpHeightgen[2,2,1,2,1,6,0];

        let k = bpHeightgen[2,2,0,1,1,6,0];
        let l = bpHeightgen[2,2,0,1,1,6,1];
        let mixed = add(bpScale[-12](k),bpScale[7](l));
        assert bpHeightClass[0,1,6](bpHeightRep(mixed)) == mixed;
        let u = bpHeightgen[3,0,7,0,0,0,0];
        assert bpHeightRep(bpScale[-9](u)) == bpint[3,0,0,-9];
        assert bpHeightClass[7,0,0](bpHeightRep(u)) == u;
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn signed_differentials_and_mixed_finite_coordinates_are_first_class() {
    let source = r"
        let a = bpHeightgen[2,2,3,1,2,8,0];
        let b = bpHeightgen[2,2,3,1,2,8,1];
        let z = zero[BPHeight(2,2,3,1,2,8)];
        show add(a,b);
        assert bpScale[2](a) == z;
        assert bpScale[2](b) == z;
        assert bpHeightClass[3,1,8](bpHeightRep(add(a,b))) == add(a,b);
        assert bpHeightD(a) == zero[BPHeight(2,2,3,-2,3,8)];
        assert bpHeightNext(bpHeightPrev(add(a,b))) == add(a,b);
        card(\x : BPHeight(2,2,3,1,2,8) => S.hi);

        let c = bpHeightgen[3,2,1,2,1,16,0];
        let d = bpHeightD(c);
        show d;
        assert bpScale[3](c) == zero[BPHeight(3,2,1,2,1,16)];
        assert bpHeightD(bpNeg(c)) == bpNeg(d);
        assert bpHeightClass[1,1,16](bpHeightRep(d)) == d;
        assert bpHeightD(d) == zero[BPHeight(3,2,1,0,3,16)];
        card(\x : BPHeight(3,2,1,1,2,16) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "BPHeight(2, 2, 3, 1, 2, 8){g0 + g1}\n",
            "card = 4 / 4\n",
            "BPHeight(3, 2, 1, 1, 2, 16){2*g0}\n",
            "card = 3 / 3\n",
        ),
    );
}

#[test]
fn huge_pages_use_stabilized_lattices() {
    let r = usize::MAX;
    let target = 1 - r as i128;
    let source = format!(
        "let a = bpHeightgen[2,2,{r},1,1,8,0];\n\
         assert bpHeightClass[{r},1,8](bpHeightRep(a)) == a;\n\
         assert bpHeightD(a) == zero[BPHeight(2,2,{r},{target},2,8)];\n\
         assert bpHeightNext(bpHeightPrev(a)) == a;"
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn typed_height_pages_reject_invalid_terms_and_transitions() {
    for (source, reason) in [
        (
            "show bpHeightClass[0,1,6](bpt[2,2,1,1,2]);",
            "cochain exceeds the requested height filtration",
        ),
        (
            "show bpHeightClass[2,2,6](bpt[2,2,1,1,2]);",
            "requires a cycle on the requested height page",
        ),
        (
            "show bpHeightClass[1,2,8](bpt[2,2,1,1,2]);",
            "requested homogeneous internal degree",
        ),
        (
            "show bpHeightNext(bpHeightgen[2,2,1,2,1,6,0]);",
            "bpHeightNext requires a cycle for the page differential",
        ),
        (
            "show bpHeightPrev(zero[BPHeight(2,2,0,1,1,6)]);",
            "bpHeightPrev requires a positive page",
        ),
        (
            "show bpHeightgen[2,2,2,-1,1,6,0];",
            "index 0 is out of range for rank 0",
        ),
        (
            r"card(\x : BPHeight(2,2,0,2,1,6) => S.hi);",
            "cannot enumerate the infinite free BP height page",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
