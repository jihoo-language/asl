mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn primitive_annihilators_and_resonant_periods() {
    assert_eq!(
        stdout(example("rv-quantum-connection")),
        concat!(
            "RVOre{(0, 5): 1, (1, 0): -2, (1, 1): -4}\n",
            "rvQGrPeriod = [1, 2, 3/8, 5/324, 35/165888]\n",
            "rvOreCoeff = -1\nrvOreCoeff = 1\n",
            "rvQGrPeriod = [1, 3, 19/32, 49/2592]\n",
            "rvOreCoeff = 1\n",
            "rvQGrPeriod = [1, 1, 1/32, 1/7776, 1/7962624]\n",
        ),
    );
}

#[test]
fn ore_products_and_action_are_noncommutative() {
    let program = r"
        let a = rvOre[2,3];
        let b = rvOre[4,1];
        let c = add(rvOre[0,5,-2],rvOre[1,0]);
        let ab = rvOreMul(a,b);
        assert ab == add(add(rvOre[6,4],rvOre[6,3,12]),add(rvOre[6,2,48],rvOre[6,1,64]));
        assert rvOreMul(b,a) == add(rvOre[6,4],rvOre[6,3,2]);
        assert add(ab,neg(ab)) == zero[RVOre];
        assert rvOreMul(ab,c) == rvOreMul(a,rvOreMul(b,c));
        let u = rvQGrBasis[2,5,[2,1]];
        let v = add(rvQGrShift[1](u),rvQGrScale[-3](rvQGrOne[2,5]));
        assert rvOreApply(ab,v) == rvOreApply(a,rvOreApply(b,v));
        assert rvOreApply(c,v) == add(rvOreApply(rvOre[0,5,-2],v),rvQGrShift[1](v));
        assert rvQGrNabla(v) == rvOreApply(rvOre[0,1],v);
        assert rvQGrTranspose(rvQGrNabla(v)) == rvQGrNabla(rvQGrTranspose(v));
        rvOreCoeff[6,3](ab);
        rvOreCoeff[6,3](rvOreMul(b,a));
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvOreCoeff = 12\nrvOreCoeff = 2\n"
    );
}

#[test]
fn projective_and_complementary_connections_agree() {
    let program = r"
        let l = rvQGrAnnihilator[1,37];
        assert l == add(rvOre[0,37],rvOre[1,0,-1]);
        assert l == rvQGrAnnihilator[36,37];
        assert rvOreApply(l,rvQGrOne[1,37]) == zero[RVQGr(1,37)];
        assert rvOreApply(l,rvQGrOne[36,37]) == zero[RVQGr(36,37)];
        assert rvQGrAnnihilator[2,5] == rvQGrAnnihilator[3,5];
        rvQGrPeriod[3,5,3];
        rvQGrPeriod[1,1000000000,1];
        rvOreCoeff[0,1000000000](rvQGrAnnihilator[1,1000000000]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvQGrPeriod = [1, 3, 19/32, 49/2592]\nrvQGrPeriod = [1, 1]\nrvOreCoeff = 1\n",
    );
}

#[test]
fn quantum_connection_diagnoses_type_and_resource_bounds() {
    let overflowing_order = format!("show rvOreMul(rvOre[0,{}],rvOre[0,1]);", usize::MAX);
    for (program, expected) in [
        ("show rvOreMul(rvOre[0,1],S.hi);", "rvOreMul"),
        (
            "show rvOreApply(rvOre[0,1],S.hi);",
            "requires a quantum Grassmannian element",
        ),
        (
            "show rvQGrNabla(rvOre[0,1]);",
            "requires a quantum Grassmannian element",
        ),
        (
            "show rvQGrAnnihilator[4,9];",
            "Schubert rank exceeds the limit of 64",
        ),
        ("rvQGrPeriod[3,6,5000];", "degree exceeds the limit of 4096"),
        (
            "rvQGrPeriod[1,1000,2];",
            "integer arithmetic is out of range",
        ),
        (
            overflowing_order.as_str(),
            "integer arithmetic is out of range",
        ),
        (
            r"show \x: RVOre => x;",
            "cannot enumerate an Ore operator ring",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
