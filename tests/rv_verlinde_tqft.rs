mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_weyl_values_and_conformal_blocks() {
    assert_eq!(
        stdout(example("rv-verlinde-tqft")),
        concat!(
            "RVCyclo(15){1, 0, -1, 1, 0, 0, 0, -1}\n",
            "RVCyclo(15){1, -1, -1, 1, -1, 1, 0, -1}\n",
            "RVCyclo(15){-1, 1, 0, 0, 1, -1}\n",
            "RVFusion(3, 2){[0, 0, 0]: 6, [0, 1, 2]: 3}\n",
            "rvConformalDim = 1\nrvConformalDim = 6\nrvConformalDim = 45\n",
            "rvConformalDim = 405\nrvConformalDim = 4050\nrvConformalDim = 42525\n",
            "rvConformalDim = 1\nrvConformalDim = 0\n",
        ),
    );
}

#[test]
fn cyclotomic_arithmetic_uses_the_exact_monic_quotient() {
    let program = r"
        let one = rvCyclo[12,0];
        let z = rvCyclo[12,1];
        assert rvCyclo[12,4] == add(rvCyclo[12,2],neg(one));
        assert rvCyclo[12,6] == neg(one);
        assert rvCyclo[12,-1] == rvCycloConj(z);
        assert rvCycloMul(z,rvCyclo[12,-1]) == one;
        assert rvCyclo[12,-170141183460469231731687303715884105728] == rvCyclo[12,4];
        assert rvCyclo[1,-10] == rvCyclo[1,0];
        let f = add(rvCycloScale[-3](rvCyclo[12,4]), rvCycloScale[2](z));
        assert add(f,neg(f)) == zero[RVCyclo(12)];
        assert rvCycloConj(rvCycloConj(f)) == f;
        rvCycloCoeff[0](f);
        rvCycloCoeff[1](f);
        rvCycloCoeff[2](f);
        rvCycloCoeff[7](f);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvCycloCoeff = 3\nrvCycloCoeff = 2\nrvCycloCoeff = -3\nrvCycloCoeff = 0\n",
    );
}

#[test]
fn handle_sewing_with_punctures_and_zero_level() {
    let program = r"
        let h = rvFusionHandle[3,2];
        let a = rvFusionBasis[2,[0,1,2]];
        let b = rvFusionBasis[2,[0,0,1]];
        assert rvQDim(rvFusionMul(h,a)) == rvCycloMul(rvQDim(h),rvQDim(a));
        assert rvVerlinde[[-7,-7,-6]](b) == rvVerlinde[[0,0,1]](b);
        rvFusionCoeff[[0,0,0]](rvFusionMul(rvFusionMul(h,h),a));
        rvConformalDim[3,2,2,[[0,1,2]]];
        rvConformalDim[3,2,2,[[0,0,1]]];
        rvConformalDim[4,0,100000000000,[[2,2,2,2],[-1,-1,-1,-1]]];
        assert rvQDim(rvFusionOne[4,0]) == rvCyclo[16,0];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvFusionCoeff = 45\nrvConformalDim = 45\nrvConformalDim = 0\nrvConformalDim = 1\n",
    );
}

#[test]
fn cyclotomic_and_sewing_inputs_have_checked_bounds() {
    for (program, expected) in [
        (
            "show rvCyclo[0,1];",
            "cyclotomic order must be between 1 and 4096",
        ),
        (
            "show rvCyclo[4097,1];",
            "cyclotomic order must be between 1 and 4096",
        ),
        ("show rvCycloMul(rvCyclo[3,1],rvCyclo[6,1]);", "rvCycloMul"),
        ("show rvQDim(S.lo);", "requires a fusion element"),
        ("show rvQDim(rvFusionOne[3,2000]);", "cyclotomic order"),
        (
            "show rvVerlinde[[0,0,3]](rvFusionOne[3,2]);",
            "outside its level alcove",
        ),
        ("rvConformalDim[3,2,1,[[0,1]]];", "must match the ring rank"),
        (
            "show rvFusionHandle[3,10000];",
            "limit of 4096 integrable weights",
        ),
        (
            r"show \x: RVCyclo(3) => x;",
            "cannot enumerate an integral cyclotomic ring",
        ),
        (
            "show neg(rvCycloScale[-170141183460469231731687303715884105728](rvCyclo[3,0]));",
            "integer arithmetic is out of range",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
