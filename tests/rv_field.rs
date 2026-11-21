mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn frobenius_extension_descent_and_trace_agree_with_specialization() {
    assert_eq!(
        stdout(example("rv-descent")),
        concat!(
            "RVHom(One, Cusp, [0, 1, 2], 7, [1, 2, 1]){2*g0}\n",
            "RVHom(One, Cusp, [0, 1, 2], 19, [1, 7, 1]){7*g0}\n",
            "RVHom(One, One, [0], 7, [1]){g0}\n",
        ),
    );
}

#[test]
fn transport_preserves_composition_and_trace_satisfies_both_projection_formulas() {
    let setup = include_str!("../examples/rv-descent.asl")
        .split("\nlet rawB1")
        .next()
        .unwrap();
    let program = format!(
        "{setup}\n{}",
        r"
        let left = rvScale[2](rvgen[Cusp,Cusp,[0,1,2],7,[1,3,1],1]);
        let right = rvScale[3](rvgen[Cusp,Cusp,[0,1,2],7,[1,3,1],2]);
        let l = rvExtend[19,6](left);
        let r = rvExtend[19,6](right);
        assert rvExtend[19,6](ecomp(left,right)) == ecomp(l,r);
        assert rvDescend[7,6](ecomp(l,r)) == ecomp(left,right);
        assert rvFrob[1](ecomp(left,right)) == ecomp(rvFrob[1](left),rvFrob[1](right));
        assert rvExtend[19,6](rvId[Cusp,[0,1,2],7,[1,3,1]]) == rvId[Cusp,[0,1,2],19,[1,7,1]];
        assert rvFrob[1](rvId[Cusp,[0,1,2],7,[1,3,1]]) == rvId[Cusp,[0,1,2],7,[1,2,1]];

        let w = rvScale[2](r);
        assert rvTrace[7,6](w) == right;
        assert rvTrace[7,6](ecomp(l,w)) == ecomp(left,rvTrace[7,6](w));
        assert rvTrace[7,6](ecomp(w,l)) == ecomp(rvTrace[7,6](w),left);
        assert rvTrace[7,6](add(w,l)) == rvTrace[7,6](w);
        assert rvTrace[7,6](ecomp(l,w)) == rvgen[Cusp,Cusp,[0,1,2],7,[1,3,1],3];
        assert rvTrace[7,6](ecomp(w,l)) == rvgen[Cusp,Cusp,[0,1,2],7,[1,3,1],0];
    "
    );
    assert_eq!(stdout(interpret(&program)), "");
}

#[test]
fn field_transport_rejects_invalid_extensions_points_and_coefficients() {
    for (program, expected) in [
        (
            "show rvExtend[19,2](rvId[One,[],7,[]]);",
            "does not satisfy the smaller modulus",
        ),
        (
            "show rvExtend[19,16](rvId[One,[],7,[]]);",
            "field coefficient 16 is out of range",
        ),
        (
            "show rvExtend[11,2](rvId[One,[],7,[]]);",
            "smaller degree to divide",
        ),
        (
            "show rvExtend[5,0](rvId[One,[],3,[]]);",
            "reducible over F2",
        ),
        (
            "show rvDescend[7,6](rvId[One,[0],19,[2]]);",
            "affine point to lie in the specified subfield",
        ),
        (
            "show rvTrace[7,6](rvId[One,[0],19,[2]]);",
            "affine point to lie in the specified subfield",
        ),
        (
            "show rvDescend[7,6](rvScale[2](rvId[One,[0],19,[6]]));",
            "every map coefficient to lie in the specified subfield",
        ),
        ("show rvFrob[1](One.pt);", "rvFrob requires an RVHom map"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}

#[test]
fn frobenius_handles_large_powers_and_degree_one_embeddings() {
    let program = format!(
        "{}\n{}\n{}",
        r"
        let x = rvScale[2](rvId[One,[0],7,[3]]);
        assert rvExtend[7,3](x) == rvFrob[1](x);
        assert rvDescend[7,3](rvFrob[1](x)) == x;
        assert rvTrace[7,3](rvFrob[1](x)) == x;
        ",
        format_args!("assert rvFrob[{}](x) == rvFrob[1](x);", usize::MAX),
        r"
        let unit = rvId[One,[0],2,[0]];
        assert rvExtend[7,0](unit) == rvId[One,[0],7,[0]];
        assert rvDescend[2,0](rvId[One,[0],7,[0]]) == unit;
        assert rvTrace[2,0](rvScale[2](rvId[One,[0],7,[0]])) == unit;
        ",
    );
    assert_eq!(stdout(interpret(&program)), "");
}
