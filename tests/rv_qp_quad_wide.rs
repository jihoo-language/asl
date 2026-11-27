mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn wide_cancellation_controls_valuation_inverse_precision_angle_and_flag_ranks() {
    assert_eq!(
        stdout(example("rv-quadratic-wide")),
        concat!(
            "rvQpQuadVal = 31\n",
            "rvQpQuadVal = 0\n",
            "RVResSpace(7, 2){[1, 6]}\n",
            "rvQpQuadResidue = 946652292209074839690256458\n",
            "rvQpQuadResidue = 19655363547886327328095587\n",
            "rvQpQuadResidue = 650756891687269553788266559\n",
            "rvQpQuadVal = 29\n",
            "RVResSpace(7, 2){[1, 2]}\n",
            "rvQpQuadVal = -2\n",
            "RVResSpace(7, 2){[1, 2]}\n",
            "rvPHodgeData[Flag, [0]] = tH 30, tN 31, degree -1, Hodge [30], Newton [31]\n",
            "rvPHodgeWeak[Flag] = weakly admissible\n",
        )
    );
}

#[test]
fn normalization_handles_a_wide_gcd_and_true_379_bit_norms() {
    let source = r"
        rvqpadic K[5,-85070591730234615865843651857942052864,1];
        let x = rvQpQuad[K,0,85070591730234615865843651857942052863,85070591730234615865843651857942052864];
        let y = rvQpQuad[K,0,-1,85070591730234615865843651857942052863];
        let one = rvQpQuad[K,1,0,1];
        assert rvQpQuadInv(x) == y;
        assert rvQpQuadMul(x,y) == one;
        assert rvQpQuadInv(y) == x;
        let u = rvQpQuad[K,85070591730234615865843651857942052864,85070591730234615865843651857942052866,75];
        rvQpQuadVal(u);
        show rvResSpecializeQuad[K,2,[0,-1]]{[one,u]};
        rvqpadic J[5,-85070591730234615865843651857942052864,4];
        rvQpQuadVal(rvQpQuad[J,85070591730234615865843651857942052864,85070591730234615865843651857942052866,75]);
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvQpQuadVal = -1\nRVResSpace(5, 2){[1, 4]}\nrvQpQuadVal = -2\n"
    );
}

#[test]
fn endpoint_coefficients_and_wide_affine_denominators_are_checked_after_normalization() {
    let source = r"
        rvqpadic M[3,-170141183460469231731687303715884105728,2];
        let u = rvQpQuad[M,-170141183460469231731687303715884105728,1,5];
        assert rvQpQuadMul(u,rvQpQuad[M,1,0,1]) == u;
        rvQpQuadVal(u);
        show rvResSpecializeQuad[M,2,[0,1]]{[rvQpQuad[M,1,0,1],u]};
        let a = rvQpQuad[M,1,0,63802943797675961899382738893456539648];
        let b = rvQpQuad[M,3,0,106338239662793269832304564822427566080];
        let target = rvQpQuad[M,7,0,159507359494189904748456847233641349120];
        assert add(a,b) == target;
        assert add(target,neg(a)) == b;
        let half = rvQpQuad[M,-1,0,2];
        assert rvQpQuadMul(half,rvQpQuad[M,-170141183460469231731687303715884105728,0,1])
            == rvQpQuad[M,85070591730234615865843651857942052864,0,1];
        rvqpadic K[7,2,3];
        let pell = rvQpQuadPow[90](rvQpQuad[K,1,1,1]);
        assert rvQpQuadInv(pell) == rvQpQuadConj(pell);
        assert rvQpQuadMul(pell,rvQpQuadInv(pell)) == rvQpQuad[K,1,0,1];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvQpQuadVal = 1\nRVResSpace(3, 2){[1, 2]}\n"
    );
}

#[test]
fn arithmetic_still_rejects_genuinely_unrepresentable_outputs() {
    let setup = "rvqpadic M[3,-170141183460469231731687303715884105728,2];";
    for statement in [
        "show rvQpQuadMul(rvQpQuad[M,170141183460469231731687303715884105727,1,1],rvQpQuad[M,2,0,1]);",
        "show rvQpQuadInv(rvQpQuad[M,0,1,1]);",
        "show neg(rvQpQuad[M,-170141183460469231731687303715884105728,0,1]);",
    ] {
        let output = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(output.contains("quadratic exact integer arithmetic is out of range"), "{output}");
    }
}
