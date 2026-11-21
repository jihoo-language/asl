mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn relative_raw_quotient_closes_the_weight_six_class_and_gives_the_exact_order() {
    assert_eq!(
        stdout(example("bp-relative")),
        concat!(
            "bprelext[2, 1, 2, 1, 6] = free 0, torsion [2]\n",
            "bporder[6] = 2\n",
            "BPRel(2, 1, 2, 1){2*t1_2}\n",
        )
    );
}

#[test]
fn quotient_splits_high_terms_and_commutes_with_the_actual_differential() {
    let source = r"
        let v1 = bpv[2,2,1,1];
        let v2 = bpv[2,2,1,2];
        let t1 = bpt[2,2,1,1,1];
        let t2 = bpt[2,2,1,1,2];
        let mixed = add(bpMul(v2,t1),bpMul(v1,t1));
        let q = bpRel[1](mixed);
        assert bpRelLift(q) == bpMul(v2,t1);
        assert bpRel[1](bpMul(v1,t1)) == zero[BPRel(2,1,2,1)];
        assert bpRel[1](bpRelLift(q)) == q;
        assert bpD(q) == bpRel[1](bpD(mixed));
        assert bpD(bpD(q)) == zero[BPRel(2,1,2,3)];
        assert add(q,bpNeg(q)) == zero[BPRel(2,1,2,1)];
        assert bpScale[-3](q) == bpNeg(add(q,add(q,q)));
        let s = bpRel[1](t2);
        let four = bpScale[4](s);
        assert bpD(bpNull[6](four)) == four;
        bporder[6](four);
    ";
    assert_eq!(stdout(interpret(source)), "bporder[6] = 1\n");
}

#[test]
fn higher_relative_weights_and_height_zero_have_the_expected_invariants() {
    let source = format!(
        "bprelext[2,0,1,0,0]; bprelext[2,0,1,1,2]; bprelext[2,1,2,1,8];\n\
         bprelext[2,1,2,2,8]; bprelext[3,1,2,1,16]; bprelext[2,7,8,16,16];\n\
         bprelext[2,1,2,16,{}];",
        usize::MAX,
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "bprelext[2, 0, 1, 0, 0] = free 0, torsion []\n\
             bprelext[2, 0, 1, 1, 2] = free 0, torsion [2]\n\
             bprelext[2, 1, 2, 1, 8] = free 0, torsion [2]\n\
             bprelext[2, 1, 2, 2, 8] = free 0, torsion [2]\n\
             bprelext[3, 1, 2, 1, 16] = free 0, torsion [3]\n\
             bprelext[2, 7, 8, 16, 16] = free 0, torsion []\n\
             bprelext[2, 1, 2, 16, {}] = free 0, torsion []\n",
            usize::MAX,
        ),
    );
}

#[test]
fn relative_shapes_and_cochains_validate_bounds_and_normalization() {
    for (source, reason) in [
        (
            "bprelext[2,2,2,1,6];",
            "relative heights must satisfy n < m",
        ),
        (
            "show bpRel[1](bpReduce[1](bpt[2,2,1,1,2]));",
            "requires an integral BP cochain",
        ),
        (
            "show bpRel[1](bpv[2,2,1,2]);",
            "normalized cochain using every arrow slot",
        ),
        (
            "show bpRelLift(bpt[2,2,1,1,2]);",
            "requires a BPRel cochain",
        ),
        (
            r"card(\q : BPRel(2,1,2,1) => S.hi);",
            "cannot enumerate an integral BP relative cochain space",
        ),
        (
            "show bpCup(bpRel[1](bpt[2,2,1,1,2]),bpRel[1](bpt[2,2,1,1,2]));",
            "bpCup requires a BP polynomial",
        ),
        (
            "bporder[6](bpRel[1](bpv[2,2,0,2]));",
            "BP Ext requires a closed cochain",
        ),
        (
            "show bpNull[8](bpRel[1](bpt[2,2,1,1,2]));",
            "requested homogeneous internal degree",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
