mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn normalized_bp_yoneda_and_signed_massey_comparison_run_in_the_cli() {
    assert_eq!(
        stdout(example("bp-ext-ring")),
        concat!(
            "BPExt(2, 2, 2, 4){g0}\n",
            "bpext[2, 2, 2, 8] = free 0, torsion [2, 2]\n",
            "BPExt(2, 2, 2, 8){g0 + g1}\n",
            "card = 4 / 4\n",
            "bporder[12] = 3\n",
            "BPExt(3, 2, 2, 12){0}\n",
        )
    );
}

#[test]
fn invariant_factors_control_exact_coordinates_and_concrete_representatives() {
    let source = r"
        let x = bpextgen[2,2,2,8,0];
        let y = bpextgen[2,2,2,8,1];
        let mixed = add(x,y);
        let z = zero[BPExt(2,2,2,8)];
        let rx = bpExtRep(x);
        let ry = bpExtRep(y);
        assert bpD(rx) == zero[BP(2,2,3)];
        assert bpD(ry) == zero[BP(2,2,3)];
        assert bpExtClass[8](rx) == x;
        assert bpExtClass[8](ry) == y;
        assert bpExtClass[8](add(rx,ry)) == mixed;
        assert add(add(x,mixed),y) == z;
        assert add(x,x) == z;
        assert add(y,y) == z;
        bporder[8](rx);
        bporder[8](ry);
        bporder[8](bpExtRep(mixed));
        let sigma = bpextgen[2,2,1,8,0];
        let twice = add(sigma,sigma);
        assert twice == bpScale[2](sigma);
        assert bpExtClass[8](bpExtRep(twice)) == twice;
        assert bpScale[-1](sigma) == bpScale[15](sigma);
        assert bpScale[16](sigma) == zero[BPExt(2,2,1,8)];
        assert bpExtClass[8](bpD(bpPow[4](bpv[2,2,0,1]))) == zero[BPExt(2,2,1,8)];
        bporder[8](bpExtRep(sigma));
        bporder[8](bpExtRep(twice));
        let longer = bpextgen[2,2,2,20,2];
        assert bpExtClass[20](bpExtRep(longer)) == longer;
        bporder[20](bpExtRep(longer));
        card(\q : BPExt(2,2,2,20) => S.hi);

    ";
    assert_eq!(
        stdout(interpret(source)),
        "bporder[8] = 2\nbporder[8] = 2\nbporder[8] = 2\nbporder[8] = 16\nbporder[8] = 8\nbporder[20] = 4\ncard = 16 / 16\n"
    );
}

#[test]
fn cup_product_respects_bilinearity_and_the_raw_cobar_representatives() {
    let source = r"
        let a = bpExtClass[2](bpt[2,2,1,1,1]);
        let b = bpExtClass[4](bpBockP[2](bpPow[2](bpv[2,2,0,1])));
        let c = bpextgen[2,2,2,8,0];
        let d = bpextgen[2,2,2,8,1];
        assert bpCup(a,b) == bpNeg(bpCup(b,a));
        assert bpCup(b,c) == bpExtClass[12](bpCup(bpExtRep(b),bpExtRep(c)));
        assert bpCup(b,add(c,d)) == add(bpCup(b,c),bpCup(b,d));
        assert bpCup(bpCup(a,a),b) == bpCup(a,bpCup(a,b));
        let t = bpt[3,2,1,1,1];
        let v = bpv[3,2,1,1];
        let w = add(bpMul(v,t),bpPow[2](t));
        let m = bpExtClass[12](bpMassey(t,t,t,w,w));
        let beta = bpExtClass[12](bpDelta(bpDelta(bpReduce[2](bpv[3,2,0,2]))));
        let shift = add(m,bpNeg(beta));
        assert shift == beta;
        assert bpExtClass[12](bpExtRep(m)) == m;
        bporder[12](bpExtRep(shift));
    ";
    assert_eq!(stdout(interpret(source)), "bporder[12] = 3\n");
}

#[test]
fn free_scalars_are_signed_and_act_modulo_torsion_without_raw_overflow() {
    let min = i128::MIN;
    let max = i128::MAX;
    let largest_degree = usize::MAX;
    let source = format!(
        r"
        let e = bpextgen[2,2,0,0,0];
        let a = bpScale[{min}](e);
        let b = bpScale[{max}](e);
        assert add(a,b) == bpNeg(e);
        assert bpExtRep(a) == bpint[2,2,0,{min}];
        assert bpExtClass[0](bpExtRep(a)) == a;
        let x = bpextgen[2,2,1,8,0];
        assert bpScale[{max}](x) == bpScale[15](x);
        assert bpCup(b,bpScale[15](x)) == x;
        assert bpCup(x,b) == bpScale[15](x);
        show bpScale[-7](e);
        card(\h : BPExt(2,2,2,{largest_degree}) => S.hi);
        ",
    );
    assert_eq!(
        stdout(interpret(&source)),
        "BPExt(2, 2, 0, 0){-7*g0}\ncard = 1 / 1\n"
    );
}

#[test]
fn bp_ext_checks_integrality_closure_types_and_enumeration_limits() {
    for (source, reason) in [
        (
            r"card(\x : BPExt(2,2,0,0) => S.hi);",
            "cannot enumerate the infinite free BP Ext group",
        ),
        (
            r"card(\x : BPExt(65537,1,1,131072) => S.hi);",
            "type exceeds the enumeration limit",
        ),
        (
            "show bpextgen[2,2,2,8,2];",
            "bpextgen index 2 is out of range for rank 2",
        ),
        (
            "show bpExtClass[2](bpReduce[1](bpt[2,2,1,1,1]));",
            "bpExtClass requires an integral BP cochain",
        ),
        (
            "show bpExtClass[6](bpt[2,2,1,1,2]);",
            "BP Ext requires a closed cochain",
        ),
        (
            "show bpExtClass[2](bpv[2,2,1,1]);",
            "normalized cochain using every arrow slot",
        ),
        (
            "show bpExtClass[4](bpt[2,2,1,1,1]);",
            "does not have the requested homogeneous internal degree",
        ),
        (
            "show bpCup(bpextgen[2,2,1,2,0],bpextgen[2,1,1,2,0]);",
            "bpCup requires the same prime, height, and ideal depth",
        ),
        (
            "show bpCup(bpextgen[2,2,1,2,0],bpt[2,2,1,1,1]);",
            "a BPExt class",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
    let overflow = diagnostic(interpret(&format!(
        "show bpCup(zero[BPExt(2,2,0,{})],zero[BPExt(2,2,0,1)]);",
        usize::MAX,
    )));
    assert!(
        overflow.contains("bpCup internal degree is out of range"),
        "{overflow}"
    );
    let error = diagnostic(interpret(&format!(
        "show bpNeg(bpScale[{}](bpextgen[2,1,0,0,0]));",
        i128::MIN
    )));
    assert!(
        error.contains("BP Ext coefficient exceeds the signed i128 range"),
        "{error}"
    );
}
