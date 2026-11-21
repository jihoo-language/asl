mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn full_integral_bracket_detects_a_singleton_and_a_proper_indeterminacy() {
    assert_eq!(
        stdout(example("bp-triple")),
        concat!(
            "bpTripleDefined = yes\n",
            "bpTripleContains = yes\n",
            "bpTripleContains = no\n",
            "BPTriple(3, 2, 2, 12){base = BPExt(3, 2, 2, 12){2*g0}, size = 1}\n",
            "bpTripleSize = 1\n",
            "bpTripleContains = yes\n",
            "bpTripleContains = no\n",
            "bpTripleContains = no\n",
            "bpTripleSize = 2^1\n",
            "bpTripleDefined = no\n",
            "BPTriple(2, 2, 2, 8){base = BPExt(2, 2, 2, 8){g1}, size = 2^1}\n",
            "bpTripleContains = no\n",
            "bpTripleContains = yes\n",
            "bpTripleContains = yes\n",
        )
    );
}

#[test]
fn equivalent_direction_generators_give_the_same_first_class_coset() {
    let program = r"
        let eta = bpextgen[2,2,1,2,0];
        let start = zero[BPExt(2,2,1,2)];
        let mid = zero[BPExt(2,2,1,4)];
        let left = bpTriple(eta,mid,start);
        let right = bpTriple(start,mid,eta);
        assert left == right;
        assert bpTripleBase(left) == zero[BPExt(2,2,2,8)];
        bpTripleContains(left,bpTripleBase(right));
        let all_zero = bpTriple(start,mid,start);
        bpTripleSize(all_zero);
        bpTripleContains(all_zero,bpextgen[2,2,2,8,0]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "bpTripleContains = yes\nbpTripleSize = 1\nbpTripleContains = no\n"
    );
}

#[test]
fn undefined_brackets_and_wrong_grades_are_diagnosed() {
    let setup = r"
        let eta = bpextgen[2,2,1,2,0];
        let zero_eta = zero[BPExt(2,2,1,2)];
        let mid = zero[BPExt(2,2,1,4)];
        let bracket = bpTriple(eta,mid,eta);
    ";
    for (statement, expected) in [
        ("show bpTriple(eta,eta);", "requires three BPExt classes"),
        ("show bpTriple(eta,eta,zero_eta);", "A cup B is nonzero"),
        ("show bpTriple(zero_eta,eta,eta);", "B cup C is nonzero"),
        (
            "bpTripleDefined(zero[BPExt(2,2,0,0)],eta,eta);",
            "positive cobar degrees",
        ),
        ("bpTripleContains(bracket,eta);", "BPExt(2, 2, 2, 8)"),
        ("bpTripleSize(eta);", "a BPTriple coset"),
        (
            "show bpTriple(eta,mid,zero[BPExt(3,2,1,4)]);",
            "same prime, height",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(expected), "{statement}\n{error}");
    }
}

#[test]
fn supplied_integral_witnesses_agree_with_the_all_nonzero_affine_example() {
    let source = r"
        let eta = bpextgen[2,2,1,2,0];
        let alpha = bpextgen[2,2,1,4,0];
        let t = bpt[2,2,1,1,1];
        let t2 = bpt[2,2,1,1,2];
        let v = bpv[2,2,1,1];
        let w = add(bpMul(v,t),bpPow[2](t));
        assert alpha == bpExtClass[4](w);
        let vvt = bpMul(bpPow[2](v),t);
        let u = add(add(bpNeg(vvt),bpMul(bpint[2,2,1,-2],bpPow[3](t))),bpMul(bpint[2,2,1,5],t2));
        let v_witness = add(add(vvt,bpPow[3](t)),bpMul(bpint[2,2,1,-3],t2));
        let checked = bpExtClass[8](bpMassey(t,w,t,u,v_witness));
        let coset = bpTriple(eta,alpha,eta);
        assert checked == bpTripleBase(coset);
        bpTripleContains(coset,checked);
        bpTripleContains(coset,zero[BPExt(2,2,2,8)]);
    ";
    assert_eq!(
        stdout(interpret(source)),
        "bpTripleContains = yes\nbpTripleContains = no\n"
    );
}
