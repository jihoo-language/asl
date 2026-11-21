mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn affine_flags_refine_hodge_positions_and_see_relative_frobenius() {
    assert_eq!(
        stdout(example("rv-iwahori")),
        concat!(
            "rvHodge = [0, 0]\n",
            "rvHodge = [0, 0]\n",
            "rvIw = [0, 1]\n",
            "rvIw = [1, 0]\n",
            "rvIwLength = 1\n",
            "rvIw = [-1, 2]\n",
            "rvIwLength = 1\n",
            "rvIw = [0, 1]\n",
            "rvIw = [1, 0]\n",
            "rvIwPosition = [1, 0]\n",
            "rvIwAt[3, 1] = [1, 0]\n",
        ),
    );
}

#[test]
fn affine_rank_dominance_handles_large_central_translates_and_components() {
    let program = r"
        qmodule X[0] { x0:0,x1:0 } {};
        let i = rvU[0](rvId[X,[0],3,[1]]);
        let translate = rvU[5000000000000000000](rvId[X,[0],3,[1]]);
        rvIw(translate);
        rvIwLength(translate);
        assert rvInIwADLV[3,1,[10000000000000000000,10000000000000000001]](translate,i) == S.hi;
        assert rvInIwADLVClosure[3,1,[10000000000000000001,10000000000000000000]](translate,i) == S.hi;
        assert rvInIwADLVClosure[3,1,[10000000000000000001,10000000000000000002]](translate,i) == S.lo;
        let rotation = add(rvU[0](rvgen[X,X,[0],3,[1],1]),rvU[1](rvgen[X,X,[0],3,[1],2]));
        rvIw(rotation);
        rvIwLength(rotation);
        assert rvInIwADLVClosure[3,1,[1,2]](i,i) == S.lo;
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvIw = [10000000000000000000, 10000000000000000001]\n",
            "rvIwLength = 0\n",
            "rvIw = [1, 2]\n",
            "rvIwLength = 0\n",
        ),
    );
}

#[test]
fn negative_and_positive_translation_powers_have_their_affine_length() {
    let program = r"
        qmodule X[0] { x0:0,x1:0 } {};
        let a = add(rvU[-3](rvgen[X,X,[0],3,[1],0]),rvU[2](rvgen[X,X,[0],3,[1],3]));
        rvIw(a);
        rvIwLength(a);
        assert rvInIwADLV[3,1,[-6,5]](a,rvU[0](rvId[X,[0],3,[1]])) == S.hi;
        let extreme = add(
            rvU[-9223372036854775808](rvgen[X,X,[0],3,[1],0]),
            rvU[9223372036854775807](rvgen[X,X,[0],3,[1],3])
        );
        rvIw(extreme);
        rvIwLength(extreme);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvIw = [-6, 5]\nrvIwLength = 5\n",
            "rvIw = [-18446744073709551616, 18446744073709551615]\n",
            "rvIwLength = 18446744073709551615\n",
        )
    );
}

#[test]
fn rank_zero_has_the_empty_affine_window_and_rank_one_checks_its_component() {
    let program = r"
        qmodule A[0] { a:1,b:0 } {0(a)=>[b]};
        let z = rvU[0](rvId[A,[0],7,[1]]);
        rvIw(z);
        rvIwLength(z);
        rvIwPosition(z,z);
        rvIwAt[3,1](z,z);
        assert rvInIwADLV[3,1,[]](z,z) == S.hi;
        assert rvInIwADLVClosure[3,1,[]](z,z) == S.hi;
        let unit = rvU[0](rvId[One,[0],7,[1]]);
        let t = rvU[-4](rvId[One,[0],7,[1]]);
        rvIwAt[3,1](t,unit);
        assert rvInIwADLVClosure[3,1,[-4]](t,unit) == S.hi;
        assert rvInIwADLVClosure[3,1,[-3]](t,unit) == S.lo;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvIw = []\nrvIwLength = 0\nrvIwPosition = []\nrvIwAt[3, 1] = []\nrvIwAt[3, 1] = [-4]\n",
    );
}

#[test]
fn malformed_affine_windows_singular_frames_and_unrational_points_fail() {
    let prefix = "qmodule X[0] { x0:0,x1:0 } {}; let i=rvU[0](rvId[X,[0],7,[1]]);";
    for (source, expected) in [
        (
            "rvIw(rvU[0](rvgen[One,X,[0],7,[1],0]));",
            "requires a Laurent Margolis endomorphism",
        ),
        (
            "rvIw(zero[RVLaurent(X,X,[0],7,[1])]);",
            "requires invertible Laurent endomorphisms",
        ),
        (
            "rvIwPosition(zero[RVLaurent(X,X,[0],7,[1])],i);",
            "requires invertible Laurent endomorphisms",
        ),
        (
            "rvIwAt[3,1](i,rvU[0](rvId[X,[0],7,[0]]));",
            "same endpoint and affine field point",
        ),
        (
            "let x=rvU[0](rvId[X,[0],7,[2]]); rvIwAt[3,1](x,x);",
            "affine point in the specified subfield",
        ),
        (
            "show rvInIwADLV[3,1,[0]](i,i);",
            "one entry per Margolis dimension",
        ),
        (
            "show rvInIwADLVClosure[3,1,[0,2]](i,i);",
            "window residues must form a permutation",
        ),
        (
            "show rvInIwADLV[3,1,[170141183460469231731687303715884105728,1]](i,i);",
            "affine window entry exceeds the signed i128 range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{prefix}\n{source}")));
        assert!(error.contains(expected), "{source}\n{error}");
    }
}
