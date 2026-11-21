mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn shintani_descends_conjugate_eigenvalues_and_distinguishes_invariant_factors() {
    assert_eq!(
        stdout(example("rv-shintani")),
        concat!(
            "rvInvariantFactors[3, 1] = [[1, 1, 1]]\n",
            "RVHom(X, X, [0], 3, [1]){g1 + g2 + g3}\n",
            "rvInvariantFactors[3, 1] = [[1, 1], [1, 1]]\n",
            "rvInvariantFactors[3, 1] = [[1, 0, 1]]\n",
        ),
    );
}

#[test]
fn shintani_handles_nonbinary_subfields_and_exact_affine_points() {
    let program = r"
        let u = rvId[One,[0],19,[6]];
        let a = rvScale[2](u);
        let n = rvFrobNorm[7,6](a);
        assert n == ecomp(a, rvFrob[2](a));
        assert n == rvScale[6](u);
        let c = rvShintani[7,6](a);
        assert c == rvScale[2](rvId[One,[0],7,[2]]);
        assert n == rvExtend[19,6](c);
        assert rvShintaniFrame[7,6](a) == u;
        rvInvariantFactors[7,6](a);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvInvariantFactors[7, 6] = [[2, 1]]\n"
    );
}

#[test]
fn shintani_zero_dimensional_norm_frame_and_factors_are_well_defined() {
    let program = r"
        qmodule A[0] { a:1,b:0 } { 0(a)=>[b] };
        let u = rvId[A,[0],7,[1]];
        assert rvFrobNorm[3,1](u) == u;
        assert rvShintaniFrame[3,1](u) == u;
        assert rvShintani[3,1](u) == rvId[A,[0],3,[1]];
        rvInvariantFactors[3,1](u);
        let v = rvId[One,[],7,[]];
        assert rvFrobNorm[7,3](v) == v;
        assert rvShintani[7,3](v) == v;
        assert rvShintaniFrame[7,3](v) == v;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvInvariantFactors[3, 1] = []\n"
    );
}

#[test]
fn shintani_requires_invertible_endomorphisms_at_rational_points() {
    let setup = "qmodule X[0] { x0:0,x1:0 } {};";
    for (program, expected) in [
        (
            "show rvFrobNorm[3,1](zero[RVHom(One,One,[0],7,[1])]);",
            "requires an invertible Margolis endomorphism",
        ),
        (
            "show rvShintani[3,1](rvgen[X,X,[0],7,[1],1]);",
            "requires an invertible Margolis endomorphism",
        ),
        (
            "show rvShintaniFrame[3,1](rvgen[X,One,[0],7,[1],0]);",
            "requires RVHom endomorphism data",
        ),
        (
            "show rvShintani[3,1](rvId[One,[0],7,[2]]);",
            "affine point in the specified subfield",
        ),
        (
            "rvInvariantFactors[3,1](rvId[One,[0],7,[2]]);",
            "affine point in the specified subfield",
        ),
        (
            "rvInvariantFactors[7,2](rvId[One,[0],19,[1]]);",
            "does not satisfy the smaller modulus",
        ),
        (
            "show rvFrobNorm[11,2](rvId[One,[0],19,[1]]);",
            "smaller degree to divide",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{program}")));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}

#[test]
fn shintani_bounds_the_polynomial_matrix_before_allocating_it() {
    let generators = (0..550)
        .map(|index| format!("x{index}:0"))
        .collect::<Vec<_>>()
        .join(",");
    let program = format!(
        "qmodule Large[0] {{{generators}}} {{}};\nshow rvShintani[3,1](rvId[Large,[0],7,[1]]);"
    );
    let error = diagnostic(interpret(&program));
    assert!(error.contains("storage limit"), "{error}");
}

#[test]
fn invariant_factors_combine_distinct_primary_blocks() {
    let matrix = [0, 7, 8, 16, 25, 31, 32, 41, 47, 48]
        .iter()
        .fold("zero[RVHom(X,X,[0],3,[1])]".to_owned(), |sum, index| {
            format!("add({sum},rvgen[X,X,[0],3,[1],{index}])")
        });
    let program = format!(
        "qmodule X[0] {{ x0:0,x1:0,x2:0,x3:0,x4:0,x5:0,x6:0 }} {{}};\nlet a = {matrix};\n{}",
        r"
        let c = rvShintani[3,1](a);
        let b = rvShintaniFrame[3,1](a);
        assert ecomp(a,b) == ecomp(b,c);
        assert ecomp(rvInverse(b),b) == rvId[X,[0],3,[1]];
        rvInvariantFactors[3,1](a);
        "
    );
    assert_eq!(
        stdout(interpret(&program)),
        "rvInvariantFactors[3, 1] = [[1, 0, 0, 1], [1, 1, 0, 1, 1]]\n",
    );
}
