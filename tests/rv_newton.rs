mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn newton_slopes_track_frobenius_cycles_cancellation_and_tate_twists() {
    assert_eq!(
        stdout(example("rv-newton")),
        concat!(
            "rvNewton[3, 1] = {0: 1, 1: 1}\n",
            "rvNewton[3, 1] = {0: 1, 1: 1}\n",
            "rvNewton[3, 1] = {1/2: 2}\n",
            "rvNewton[3, 1] = {1/2: 2}\n",
            "rvNewton[3, 1] = {1/3: 3}\n",
            "rvNewton[3, 1] = {-5/3: 3}\n",
        ),
    );
}

#[test]
fn laurent_maps_form_a_category_and_frobenius_transports_the_point_and_coefficients() {
    let program = r"
        qmodule X[0] { x0:0,x1:0 } {};
        let f = rvU[-3](rvgen[One,X,[0],7,[2],1]);
        let g = rvU[2](rvScale[2](rvgen[X,One,[0],7,[2],1]));
        let c = rvU[-1](rvScale[2](rvId[One,[0],7,[2]]));
        let zero = zero[RVLaurent(One,One,[0],7,[2])];
        assert ecomp(g,f) == c;
        assert add(c,c) == zero;
        assert rvScale[0](c) == zero;
        assert rvScale[3](c) == rvU[-1](rvId[One,[0],7,[2]]);
        assert rvFrob[1](c) == rvU[-1](rvScale[3](rvId[One,[0],7,[3]]));
        assert rvFrob[1](ecomp(g,f)) == ecomp(rvFrob[1](g),rvFrob[1](f));
        assert ecomp(rvU[0](rvId[X,[0],7,[2]]),f) == f;
        show c;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "RVLaurent(One, One, [0], 7, [2]){2*u^-1*g0}\n",
    );
}

#[test]
fn direct_sums_and_non_monomial_determinants_have_exact_signed_slopes() {
    let program = r"
        qmodule X[0] { x0:0,x1:0,x2:0 } {};
        let off = add(rvU[0](rvgen[X,X,[0],3,[1],1]),rvU[1](rvgen[X,X,[0],3,[1],3]));
        let direct = add(off,rvU[-1](rvgen[X,X,[0],3,[1],8]));
        rvNewton[3,1](direct);
        let unit = rvId[One,[0],3,[1]];
        let nonunit = add(rvU[0](unit),rvU[1](unit));
        rvNewton[3,1](nonunit);
        assert rvUNorm[3,1](nonunit) == nonunit;
        assert rvNewtonCoeff[3,1,0](nonunit) == nonunit;
        assert rvNewtonCoeff[3,1,99](nonunit) == zero[RVLaurent(One,One,[0],3,[1])];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvNewton[3, 1] = {-1: 1, 1/2: 2}\nrvNewton[3, 1] = {0: 1}\n",
    );
}

#[test]
fn newton_handles_relative_frobenius_and_empty_margolis_homology() {
    let program = r"
        qmodule A[0] { a:1,b:0 } { 0(a)=>[b] };
        let z = rvU[0](rvId[A,[0],19,[6]]);
        rvNewton[7,6](z);
        assert rvUNorm[7,6](z) == z;
        assert rvNewtonCoeff[7,6,0](z) == rvU[0](rvId[One,[0],7,[2]]);
        let a = rvU[-5](rvScale[2](rvId[One,[0],19,[6]]));
        rvNewton[7,6](a);
        assert rvUNorm[7,6](a) == rvU[-10](rvScale[6](rvId[One,[0],19,[6]]));
        assert rvNewtonCoeff[7,6,0](a) == rvU[-10](rvScale[2](rvId[One,[0],7,[2]]));
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvNewton[7, 6] = {}\nrvNewton[7, 6] = {-5: 1}\n",
    );
}

#[test]
fn newton_rejects_singular_matrices_incompatible_points_and_exponent_overflow() {
    for (program, expected) in [
        ("rvNewton[3,1](zero[RVLaurent(One,One,[0],3,[1])]);", "requires an invertible Laurent endomorphism"),
        ("show rvUNorm[3,1](zero[RVLaurent(One,One,[0],3,[1])]);", "requires an invertible Laurent endomorphism"),
        ("rvNewton[3,1](rvU[0](rvId[One,[0],7,[2]]));", "affine point in the specified subfield"),
        ("show ecomp(rvU[0](rvId[One,[0],7,[1]]),rvU[0](rvId[One,[0],7,[2]]));", "same affine Margolis field point"),
        ("show ecomp(rvU[9223372036854775807](rvId[One,[0],3,[1]]),rvU[1](rvId[One,[0],3,[1]]));", "Laurent exponent is out of the signed 64-bit range"),
        ("rvNewton[3,1](rvU[9223372036854775807](rvId[One,[0],7,[1]]));", "Laurent exponent is out of the signed 64-bit range"),
        ("show rvU[9223372036854775808](rvId[One,[0],3,[1]]);", "Laurent exponent is out of range"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
    assert_eq!(
        stdout(interpret(
            "rvNewton[3,1](rvU[-9223372036854775808](rvId[One,[0],3,[1]]));"
        )),
        "rvNewton[3, 1] = {-9223372036854775808: 1}\n",
    );
}
