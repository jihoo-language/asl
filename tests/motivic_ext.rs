mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn motivic_cubic_relation_and_tau_torsion_hold_in_ext() {
    assert_eq!(
        stdout(example("motivic-ext")),
        concat!(
            "MExt(3, 6, 2){g0}\n",
            "MExt(4, 8, 4){g0}\n",
            "mextdim[4, 8, 4] = 1\n",
            "mextdim[4, 8, 3] = 0\n",
        ),
    );
}

#[test]
fn motivic_representatives_and_realization_are_multiplicative() {
    let source = r"
        let one = mextgen[0,0,0,0];
        let tau = tauShift[1](one);
        let h0 = mextgen[1,1,0,0];
        let h1 = mextgen[1,2,1,0];
        let h2 = mextgen[1,4,2,0];
        assert extRep(tau) == mcword[1]{};
        assert yoneda(one,h0) == h0;
        assert yoneda(h1,one) == h1;
        assert tauShift[3](h0) == yoneda(tauShift[2](tau),h0);
        assert extClass(extRep(tauShift[5](h1))) == tauShift[5](h1);
        assert tauShift[1](add(h1,h1)) == zero[MExt(1,2,0)];
        assert extClass(concat(extRep(h1),extRep(h2))) == yoneda(h1,h2);
        assert realize(yoneda(h1,h2)) == yoneda(realize(h1),realize(h2));
        assert realize(tauShift[8](h0)) == realize(h0);
        assert realize(tau) == extgen[0,0,0];
        assert realize(h0) == extgen[1,1,0];
        assert realize(h1) == extgen[1,2,0];
        assert realize(h2) == extgen[1,4,0];
        let g = \x : MExt(1,2,1) => tauShift[2](x);
        assert g(h1) == tauShift[2](h1);
        card(\x : MExt(1,2,1) => S.hi);
    ";
    assert_eq!(stdout(interpret(source)), "card = 2 / 2\n");
}

#[test]
fn motivic_ext_bidegrees_slice_at_signed_weights() {
    let source = r"
        mextdim[3,6,3];
        mextdim[3,6,2];
        mextdim[3,6,1];
        mextdim[3,6,-2];
        mextdim[0,0,-10];
        mextdim[0,0,1];
        mextdim[3,2,-2];
        mextdim[1,2,2];
        mextdim[12,24,12];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "mextdim[3, 6, 3] = 1\n",
            "mextdim[3, 6, 2] = 1\n",
            "mextdim[3, 6, 1] = 1\n",
            "mextdim[3, 6, -2] = 1\n",
            "mextdim[0, 0, -10] = 1\n",
            "mextdim[0, 0, 1] = 0\n",
            "mextdim[3, 2, -2] = 0\n",
            "mextdim[1, 2, 2] = 0\n",
            "mextdim[12, 24, 12] = 1\n",
        ),
    );
}

#[test]
fn motivic_ext_rejects_wrong_types_noncocycles_and_nonboundaries() {
    for (source, reason) in [
        (
            "show extClass(mcword[0]{([],[0,1])});",
            "extClass requires a closed cobar cochain",
        ),
        (
            "show extPrimitive(mcword[0]{([1],[]),([1],[]),([1],[]),([1],[])});",
            "extPrimitive requires a cobar boundary",
        ),
        (
            "show yoneda(mextgen[1,2,1,0],extgen[1,2,0]);",
            "yoneda requires an MExt class",
        ),
        ("show mextgen[1,2,2,0];", "index 0"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
