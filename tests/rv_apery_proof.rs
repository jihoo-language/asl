#[allow(dead_code)]
mod common;

use common::{example, stdout};

#[test]
fn computed_quantum_operator_certifies_the_laurent_apery_period() {
    assert_eq!(
        stdout(example("rv-grassmannian-apery-proof")),
        concat!(
            "RVOre{(0, 2): 1, (1, 0): -3, (1, 1): -11, (1, 2): -11, ",
            "(2, 0): -1, (2, 1): -2, (2, 2): -1}\n",
            "rvQGrRegularizedRhs = {}\n",
            "rvQGrRegularized = [1, 3, 19, 147, 1251, 11253]\n",
            "rvCTWZCheck = certified Laurent identity and recurrence (2 initial coefficients checked)\n",
            "rvCTExact[Mirror, 0] = 1\n",
            "rvCTExact[Mirror, 1] = 3\n",
            "rvCTExact[Mirror, 2] = 19\n",
            "rvCTExact[Mirror, 3] = 147\n",
            "RVCTResidue(3, 2){6}\n",
            "RVCTResidue(3, 3){19}\n",
        )
    );
}
