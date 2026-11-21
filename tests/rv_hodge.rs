mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn supersingular_lattices_distinguish_hodge_newton_and_adlv_closure() {
    assert_eq!(
        stdout(example("rv-adlv")),
        concat!(
            "rvNewton[3, 1] = {1/2: 2}\n",
            "rvHodge = [0, 1]\n",
            "rvPosition = [0, 2]\n",
            "rvMuAt[3, 1] = [0, 1]\n",
            "rvMuAt[3, 1] = [-1, 2]\n",
            "rvPosition = [0, 0]\n",
            "rvMuAt[3, 1] = [-1, 2]\n",
            "rvHodge = [0, 1]\n",
        ),
    );
}

#[test]
fn relative_frobenius_changes_the_answer_and_fixes_subfield_coefficients() {
    let program = r"
        qmodule X[0] { x0:0,x1:0 } {};
        let a = rvU[0](rvId[X,[0],19,[6]]);
        let upper = rvgen[X,X,[0],19,[6],2];
        let fixed = add(a,rvU[-1](rvScale[6](upper)));
        let changed = add(a,rvU[-1](rvScale[2](upper)));
        rvMuAt[7,6](a,fixed);
        rvMuAt[7,6](a,changed);
        assert rvInADLV[7,6,[-1,1]](a,changed) == S.hi;
        assert rvInADLVClosure[7,6,[-2,2]](a,changed) == S.hi;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvMuAt[7, 6] = [0, 0]\nrvMuAt[7, 6] = [-1, 1]\n",
    );
}

#[test]
fn hodge_exponents_include_multiplicity_and_allow_rank_zero() {
    let program = r"
        qmodule X[0] { x0:0,x1:0,x2:0 } {};
        let a = add(
            rvU[2](add(rvgen[X,X,[0],3,[1],0],rvgen[X,X,[0],3,[1],4])),
            rvU[-1](rvgen[X,X,[0],3,[1],8])
        );
        rvHodge(a);
        qmodule A[0] { a:1,b:0 } { 0(a)=>[b] };
        let z = rvU[0](rvId[A,[0],7,[1]]);
        rvHodge(z);
        rvPosition(z,z);
        rvMuAt[3,1](z,z);
        assert rvInADLV[3,1,[]](z,z) == S.hi;
        assert rvInADLVClosure[3,1,[]](z,z) == S.hi;
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvHodge = [-1, 2, 2]\n",
            "rvHodge = []\n",
            "rvPosition = []\n",
            "rvMuAt[3, 1] = []\n",
        ),
    );
}

#[test]
fn positions_are_exact_when_exponent_differences_exceed_i64() {
    let program = r"
        let unit = rvId[One,[0],3,[1]];
        rvHodge(rvU[-9223372036854775808](unit));
        rvPosition(rvU[-9223372036854775808](unit),rvU[9223372036854775807](unit));
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvHodge = [-9223372036854775808]\nrvPosition = [18446744073709551615]\n",
    );
}

#[test]
fn hodge_and_adlv_reject_invalid_matrices_points_and_coweights() {
    let setup = "qmodule X[0] { x0:0,x1:0 } {};";
    for (program, expected) in [
        (
            "rvHodge(zero[RVLaurent(X,X,[0],3,[1])]);",
            "requires invertible Laurent endomorphisms",
        ),
        (
            "rvHodge(rvU[0](rvgen[One,X,[0],3,[1],0]));",
            "requires a Laurent Margolis endomorphism",
        ),
        (
            "rvPosition(rvU[0](rvId[X,[0],7,[1]]),rvU[0](rvId[X,[0],7,[2]]));",
            "same endpoint and affine field point",
        ),
        (
            "rvPosition(rvU[0](rvId[X,[0],3,[1]]),zero[RVLaurent(X,X,[0],3,[1])]);",
            "requires invertible Laurent endomorphisms",
        ),
        (
            "rvPosition(zero[RVLaurent(X,X,[0],3,[1])],rvU[0](rvId[X,[0],3,[1]]));",
            "requires invertible Laurent endomorphisms",
        ),
        (
            "let x=rvU[0](rvId[X,[0],7,[2]]); rvMuAt[3,1](x,x);",
            "affine point in the specified subfield",
        ),
        (
            "let x=rvU[0](rvId[X,[0],7,[1]]); show rvInADLV[3,1,[1,0]](x,x);",
            "coweight must be weakly increasing",
        ),
        (
            "let x=rvU[0](rvId[X,[0],7,[1]]); show rvInADLVClosure[3,1,[0]](x,x);",
            "one entry per Margolis dimension",
        ),
    ] {
        let message = diagnostic(interpret(&format!("{setup}\n{program}")));
        assert!(message.contains(expected), "{program}\n{message}");
    }
}

#[test]
fn polynomial_matrix_storage_is_bounded() {
    let generators = (0..550)
        .map(|index| format!("x{index}:0"))
        .collect::<Vec<_>>()
        .join(",");
    let program =
        format!("qmodule Large[0] {{{generators}}} {{}};\nrvHodge(rvU[0](rvId[Large,[0],3,[1]]));");
    let message = diagnostic(interpret(&program));
    assert!(message.contains("storage limit"), "{message}");
}
