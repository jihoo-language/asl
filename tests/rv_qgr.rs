mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn quantum_schubert_and_gromov_witten_examples() {
    assert_eq!(
        stdout(example("rv-quantum-grassmannian")),
        concat!(
            "RVQGr(3, 6){(0, [3, 3, 3]): 1, (1, [1, 1, 1]): 1, (1, [2, 1]): 2, (1, [3]): 1}\n",
            "rvQGrCoeff = 2\nrvQGrGW = 2\nrvQGrGW = 1\n",
            "RVQGr(36, 37){(1, []): 1}\n",
        ),
    );
}

#[test]
fn signed_quantum_ring_is_associative_and_transposes() {
    let program = r"
        let u = rvQGrOne[2,5];
        let a = rvQGrBasis[2,5,[2,1]];
        let b = rvQGrBasis[2,5,[3,2]];
        let c = rvQGrBasis[2,5,[1,1]];
        assert rvQGrBasis[2,5,[0,0]] == u;
        assert rvQGrBasis[2,5,[1,0]] == rvQGrBasis[2,5,[1]];
        let x = add(rvQGrScale[3](a),rvQGrShift[2](b));
        let y = add(c,rvQGrShift[1](neg(a)));
        assert rvQGrMul(rvQGrMul(x,y),b) == rvQGrMul(x,rvQGrMul(y,b));
        assert rvQGrMul(x,y) == rvQGrMul(y,x);
        assert add(x,neg(x)) == zero[RVQGr(2,5)];
        assert rvQGrTranspose(rvQGrMul(x,y)) == rvQGrMul(rvQGrTranspose(x),rvQGrTranspose(y));
        assert rvQGrTranspose(rvQGrTranspose(x)) == x;
        assert rvQGrClassical(rvQGrMul(x,y)) == rvQGrClassical(rvQGrMul(rvQGrClassical(x),rvQGrClassical(y)));
        rvQGrCoeff[2,[3,2]](x);
        rvQGrCoeff[0,[3,2]](x);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvQGrCoeff = 1\nrvQGrCoeff = 0\n"
    );
}

#[test]
fn rim_hooks_reduce_multiple_powers_and_projective_space() {
    let program = r"
        let point = rvQGrBasis[3,6,[3,3,3]];
        assert rvQGrMul(point,point) == rvQGrShift[3](rvQGrOne[3,6]);
        rvQGrGW[3,6,3,[3,3,3],[3,3,3],[3,3,3]];
        rvQGrGW[3,6,2,[3,3,3],[3,3,3],[3,3,3]];
        let h = rvQGrBasis[1,37,[1]];
        let h2 = rvQGrMul(h,h);
        let h4 = rvQGrMul(h2,h2);
        let h8 = rvQGrMul(h4,h4);
        let h16 = rvQGrMul(h8,h8);
        let h32 = rvQGrMul(h16,h16);
        assert rvQGrMul(h32,rvQGrMul(h4,h)) == rvQGrShift[1](rvQGrOne[1,37]);
        rvQGrGW[1,37,1,[36],[1],[36]];
        rvQGrGW[36,37,1,[1],[1],[1]];
        let x = rvQGrBasis[1,1000000000,[999999999]];
        let y = rvQGrBasis[1,1000000000,[1]];
        assert rvQGrMul(x,y) == rvQGrShift[1](rvQGrOne[1,1000000000]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvQGrGW = 1\nrvQGrGW = 0\nrvQGrGW = 1\nrvQGrGW = 0\n",
    );
}

#[test]
fn qgr_reports_bad_partitions_and_checked_bounds() {
    let exponent_overflow = format!(
        "show rvQGrShift[1](rvQGrShift[{}](rvQGrOne[2,4]));",
        usize::MAX,
    );
    for (program, expected) in [
        ("show rvQGrOne[0,9];", "requires 0 < k < n"),
        ("show rvQGrOne[9,9];", "requires 0 < k < n"),
        (
            "show rvQGrOne[25,50];",
            "requires min(k,n-k) to be at most 24",
        ),
        ("show rvQGrBasis[2,5,[1,2]];", "must be weakly decreasing"),
        (
            "show rvQGrBasis[2,5,[4]];",
            "outside the Grassmannian rectangle",
        ),
        ("show rvQGrMul(rvQGrOne[1,4],rvQGrOne[2,4]);", "rvQGrMul"),
        (
            "rvQGrCoeff[2,[1]](S.lo);",
            "requires a quantum Grassmannian element",
        ),
        (
            r"show \a: RVQGr(2,4) => a;",
            "cannot enumerate a quantum Grassmannian ring",
        ),
        (
            "show rvQGrTranspose(rvQGrBasis[1,1000000,[300001]]);",
            "transposed Schubert partition exceeds",
        ),
        (
            exponent_overflow.as_str(),
            "integer arithmetic is out of range",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
