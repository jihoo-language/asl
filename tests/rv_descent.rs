mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn constructive_frames_and_rectangular_twisted_descent() {
    assert_eq!(
        stdout(example("rv-hilbert90")),
        concat!(
            "RVHom(X, X, [0], 7, [1]){g0 + 3*g2 + g3}\n",
            "RVHom(X, Y, [0], 3, [1]){g1}\n",
            "RVHom(X, X, [0], 11, [1]){3*g0 + g1 + g2 + 5*g3}\n",
        ),
    );
}

#[test]
fn inverses_accept_isomorphisms_between_distinct_endpoints_and_zero_dimension() {
    let program = r"
        qmodule X[0] { x0:0, x1:0 } {};
        qmodule Y[0] { y0:3, y1:5 } {};
        let f = add(rvgen[X,Y,[0],7,[0],1], rvScale[3](rvgen[X,Y,[0],7,[0],2]));
        let g = rvInverse(f);
        let expected = add(rvScale[2](rvgen[Y,X,[0],7,[0],1]), rvgen[Y,X,[0],7,[0],2]);
        assert g == expected;
        assert ecomp(g,f) == rvId[X,[0],7,[0]];
        assert ecomp(f,g) == rvId[Y,[0],7,[0]];

        qmodule A[0] { a:1, b:0 } {0(a)=>[b]};
        let u = rvId[A,[0],7,[1]];
        assert rvInverse(u) == u;
        assert rvHilbert90[3,1](u) == u;
        assert rvTwistedDescend[3,1](u,u,u) == rvId[A,[0],3,[1]];
        assert rvTwistedExtend[3,1](u,rvId[A,[0],3,[1]],u) == u;
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn descent_checks_norm_intertwining_and_points() {
    let setup = "qmodule X[0] { x0:0,x1:0 } {};";
    for (program, expected) in [
        ("show rvInverse(zero[RVHom(One,One,[0],7,[1])]);", "requires an invertible Margolis map"),
        ("show rvInverse(zero[RVHom(X,One,[0],7,[1])]);", "requires equal Margolis dimensions"),
        ("show rvHilbert90[7,6](rvScale[2](rvId[One,[0],19,[1]]));", "cocycle norm to equal the identity"),
        ("show rvHilbert90[3,1](rvId[One,[0],7,[2]]);", "affine point in the specified subfield"),
        ("show rvHilbert90[3,1](rvgen[X,One,[0],7,[1],0]);", "requires RVHom endomorphism data"),
        ("let u = rvId[One,[0],7,[1]]; show rvTwistedDescend[3,1](u,rvScale[2](u),u);", "intertwines the two descent data"),
        ("let u = rvId[One,[0],7,[1]]; show rvTwistedExtend[3,1](u,rvId[One,[0],3,[0]],u);", "matching endpoints and affine field point"),
        ("let u = rvId[One,[0],7,[1]]; show rvTwistedDescend[3,1](u,u,rvId[One,[0],7,[0]]);", "the same affine Margolis field point"),
        ("let u = rvId[One,[0],7,[1]]; show rvTwistedExtend[3,1](u,rvId[One,[0],3,[1]],zero[RVHom(One,One,[0],7,[1])]);", "cocycle norm to equal the identity"),
    ] {
        let message = diagnostic(interpret(&format!("{setup}\n{program}")));
        assert!(message.contains(expected), "{program}\n{message}");
    }
}

#[test]
fn twisted_data_over_a_nonbinary_rational_point_use_the_relative_frobenius() {
    let program = r"
        let u = rvId[One,[0],19,[6]];
        let change = rvScale[2](u);
        let a = ecomp(change, rvFrob[2](rvInverse(change)));
        assert ecomp(a, rvFrob[2](a)) == u;
        let b = rvHilbert90[7,6](a);
        assert ecomp(a, rvFrob[2](b)) == b;

        let f = rvScale[3](rvId[One,[0],7,[2]]);
        let large = rvTwistedExtend[7,6](u, f, a);
        assert rvFrob[2](large) == ecomp(large, a);
        assert rvTwistedDescend[7,6](u,large,a) == f;
        assert ecomp(large,b) == rvExtend[19,6](f);
    ";
    assert_eq!(stdout(interpret(program)), "");
}
