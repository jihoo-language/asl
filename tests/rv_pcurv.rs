mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn characteristic_and_nontrivial_nilpotent_vector() {
    assert_eq!(
        stdout(example("rv-quantum-pcurvature")),
        concat!(
            "RVQGrP(2, 2, 4){(0, [1]): 1, (0, [1, 1]): 1, (0, [2]): 1}\n",
            "RVQGrP(2, 2, 4){(0, [1, 1]): 1, (0, [2]): 1}\n",
            "rvQGrCurvatureChar = [[], [], [], [], [], [], [(0, 1)]]\n",
            "RVQGrP(3, 2, 4){(0, [1]): 2, (0, [2, 1]): 2}\n",
            "rvQGrCurvatureCoeff = [(3, 2)]\n",
            "rvQGrCurvatureHorizontal = [[], [], [(1, 2)], [], [], [], [(0, 1)]]\n",
            "rvQGrCurvatureCoeff = [(5, 4)]\nrvQGrCurvatureCoeff = [(10, 4)]\n",
            "rvQGrCurvatureCoeff = [(15, 1)]\nrvQGrCurvatureCoeff = [(5, 4)]\n",
        ),
    );
}

#[test]
fn exact_reduction_is_compatible_with_connection_and_action() {
    let program = r"
        let h = rvQGrBasis[3,6,[3,1,1]];
        let u = rvQGrOne[3,6];
        let integral = add(rvQGrScale[-3](h),rvQGrShift[1](rvQGrScale[7](u)));
        let v = rvQGrMod[5](integral);
        assert v == add(rvQGrScale[2](rvQGrPBasis[5,3,6,[3,1,1]]),rvQGrShift[1](rvQGrScale[2](rvQGrPBasis[5,3,6,[]])));
        assert rvQGrMod[5](rvQGrNabla(integral)) == rvQGrNabla(v);
        let c = rvQGrCurvature(v);
        let n1 = rvQGrNabla(v);
        let n2 = rvQGrNabla(n1);
        let n3 = rvQGrNabla(n2);
        let n4 = rvQGrNabla(n3);
        let n5 = rvQGrNabla(n4);
        assert c == add(n5,neg(n1));
        assert c == rvQGrCurvatureSatake(v);
        assert rvQGrCurvature(rvQGrShift[7](v)) == rvQGrShift[7](c);
        assert rvQGrNabla(c) == rvQGrCurvature(n1);
        assert rvQGrTranspose(c) == rvQGrCurvature(rvQGrTranspose(v));
        assert add(c,neg(c)) == zero[RVQGrP(5,3,6)];
        let minimum = rvQGrScale[-170141183460469231731687303715884105728](u);
        assert rvQGrMod[3](minimum) == rvQGrPBasis[3,3,6,[]];
        rvQGrCoeff[1,[]](v);
        rvQGrCoeff[0,[3,1,1]](rvQGrClassical(v));
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvQGrCoeff = 2\nrvQGrCoeff = 2\n"
    );
}

#[test]
fn projective_power_and_prime_dividing_n() {
    let program = r"
        let u = rvQGrPBasis[5,1,5,[]];
        let c1 = rvQGrCurvature(u);
        let c2 = rvQGrCurvature(c1);
        let c3 = rvQGrCurvature(c2);
        let c4 = rvQGrCurvature(c3);
        let c5 = rvQGrCurvature(c4);
        assert c5 == rvQGrShift[5](u);
        let v = rvQGrPBasis[3,4,6,[2,1,1,1]];
        assert rvQGrCurvature(v) == rvQGrCurvatureSatake(v);
        assert rvQGrTranspose(rvQGrCurvature(v)) == rvQGrCurvature(rvQGrTranspose(v));
        rvQGrCurvatureCoeff[2,3,6,8];
        rvQGrCurvatureCoeff[3,3,6,2];
        rvQGrCurvatureCoeff[7,1,5,0];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvQGrCurvatureCoeff = [(4, 1)]\nrvQGrCurvatureCoeff = [(9, 2)]\nrvQGrCurvatureCoeff = [(7, 6)]\n",
    );
}

#[test]
fn invalid_prime_types_and_budgets_are_diagnosed() {
    for (program, expected) in [
        ("show rvQGrPBasis[1,1,5,[]];", "is not prime"),
        ("rvQGrCurvatureChar[4,2,4];", "is not prime"),
        ("show rvQGrCurvatureCol[9,2,4,[]];", "is not prime"),
        (
            "show rvQGrCurvature(rvQGrOne[2,4]);",
            "requires a finite-field Schubert section",
        ),
        (
            "show rvQGrMod[3](rvQGrPBasis[5,2,4,[]]);",
            "requires an integral quantum Grassmannian element",
        ),
        (
            "show add(rvQGrPBasis[3,2,4,[]],rvQGrPBasis[5,2,4,[]]);",
            "RVQGrP(3, 2, 4)",
        ),
        (
            "show rvQGrCurvatureCol[4099,2,4,[]];",
            "characteristic exceeds the limit of 4096",
        ),
        (
            "rvQGrCurvatureChar[7,4,9];",
            "Schubert rank exceeds the limit of 64",
        ),
        (
            "show rvQGrCurvatureSatake(rvQGrPBasis[3,1,257,[]]);",
            "ambient dimension exceeds the limit of 256",
        ),
        (
            r"show \v: RVQGrP(3,2,4) => v;",
            "cannot enumerate a polynomial Schubert module",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
