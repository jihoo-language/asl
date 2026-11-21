mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn adjoint_fusion_changes_at_each_affine_level() {
    assert_eq!(
        stdout(example("rv-verlinde")),
        concat!(
            "RVFusion(3, 2){[0, 0, 0]: 1, [0, 1, 2]: 1}\n",
            "RVFusion(3, 3){[0, 0, 0]: 1, [0, 0, 3]: 1, [0, 1, 2]: 2, [0, 3, 3]: 1}\n",
            "RVFusion(3, 4){[0, 0, 0]: 1, [0, 0, 3]: 1, [0, 1, 2]: 2, [0, 2, 4]: 1, [0, 3, 3]: 1}\n",
            "rvFusionCoeff = 1\n",
        ),
    );
}

#[test]
fn signed_fusion_is_distributive_associative_and_compatible_with_duality() {
    let program = r"
        let u = rvFusionOne[3,2];
        let a = rvFusionBasis[2,[-10,-10,-9]];
        let b = rvFusionDual(a);
        assert b == rvFusionBasis[2,[0,1,1]];
        let c = rvFusionBasis[2,[0,0,2]];
        let v = add(rvFusionScale[2](a), neg(b));
        let w = add(b, rvFusionScale[-3](c));
        let adj = rvFusionBasis[2,[0,1,2]];
        let z = add(adj,rvFusionScale[-2](u));
        assert rvFusionMul(v,add(w,z)) == add(rvFusionMul(v,w),rvFusionMul(v,z));
        assert rvFusionMul(rvFusionMul(v,w),z) == rvFusionMul(v,rvFusionMul(w,z));
        assert rvFusionMul(v,w) == rvFusionMul(w,v);
        assert rvFusionDual(rvFusionMul(v,w)) == rvFusionMul(rvFusionDual(v),rvFusionDual(w));
        assert add(v,neg(v)) == zero[RVFusion(3,2)];
        rvFusionCoeff[[-15,-15,-14]](v);
        rvFusionCoeff[[0,1,1]](v);
        rvFusionCoeff[[0,1,2]](v);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvFusionCoeff = 2\nrvFusionCoeff = -1\nrvFusionCoeff = 0\n",
    );
}

#[test]
fn rank_two_interval_and_level_zero_unit_are_exact() {
    let program = r"
        let x = rvFusionBasis[5,[-5,-1]];
        assert rvFusionMul(x,x) == add(rvFusionOne[2,5], rvFusionBasis[5,[0,2]]);
        let y = rvFusionBasis[3,[0,3]];
        assert rvFusionMul(y,y) == rvFusionOne[2,3];
        assert y == rvFusionDual(y);
        let u = rvFusionBasis[0,[-170141183460469231731687303715884105728,-170141183460469231731687303715884105728]];
        assert u == rvFusionOne[2,0];
        let c = rvFusionScale[-170141183460469231731687303715884105728](u);
        assert rvFusionMul(c,u) == c;
        rvFusionCoeff[[0,0]](c);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvFusionCoeff = -170141183460469231731687303715884105728\n",
    );
}

#[test]
fn fusion_validates_alcoves_types_and_ring_arithmetic() {
    for (program, expected) in [
        (
            "show rvFusionOne[1,2];",
            "fusion rank must be between 2 and 24",
        ),
        (
            "show rvFusionBasis[2,[0,2,1]];",
            "must be weakly increasing",
        ),
        ("show rvFusionBasis[2,[5,5,8]];", "outside its level alcove"),
        ("show rvFusionBasis[0,[0,0,1]];", "outside its level alcove"),
        (
            "rvFusionCoeff[[0]](rvFusionOne[2,1]);",
            "must match the ring rank",
        ),
        (
            "show rvFusionMul(rvFusionOne[2,1],rvFusionOne[2,2]);",
            "rvFusionMul",
        ),
        (
            "show rvFusionDual(S.lo);",
            "rvFusionDual requires a fusion element",
        ),
        (
            r"show \x: RVFusion(3,1) => x;",
            "cannot enumerate a Verlinde fusion ring",
        ),
        (
            "show neg(rvFusionScale[-170141183460469231731687303715884105728](rvFusionOne[2,1]));",
            "integer arithmetic is out of range",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
