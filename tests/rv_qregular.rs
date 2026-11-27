mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn apery_order_two_and_projective_forcing_are_certified() {
    assert_eq!(
        stdout(example("rv-quantum-regularization")),
        concat!(
            "RVOre{(0, 2): 1, (1, 0): -3, (1, 1): -11, (1, 2): -11, ",
            "(2, 0): -1, (2, 1): -2, (2, 2): -1}\n",
            "rvQGrRegularizedRhs = {}\n",
            "rvQGrRegularized = [1, 3, 19, 147, 1251, 11253]\n",
            "RVOre{(0, 0): 1, (1, 0): -1}\n",
            "rvQGrRegularizedRhs = {0: 1}\n",
            "rvOreCoeff = 1\n",
            "rvQGrRegularized = [1, 6, 126, 3948, 149310]\n",
        )
    );
}

#[test]
fn factorial_respects_ore_products_and_left_cancellation() {
    let program = r"
        let l = add(rvOre[2,3,-7],rvOre[0,2,3]);
        let t = rvOreFactorial[2](l);
        assert t == add(rvOre[0,2,3],add(rvOre[2,7,-7],add(rvOre[2,6,-42],add(rvOre[2,5,-91],add(rvOre[2,4,-84],rvOre[2,3,-28])))));
        let m = add(rvOre[1,2,3],rvOre[3,0,-2]);
        assert rvOreFactorial[0](m) == m;
        assert rvOreFactorial[2](rvOreMul(l,m)) == rvOreMul(t,rvOreFactorial[2](m));
        let one = rvOre[0,0];
        let zero = zero[RVOre];
        assert rvOreFactorial[3](zero) == zero;
        assert rvOreEulerFactor(zero) == one;
        assert rvOreEulerReduce(zero) == zero;
        let d = add(rvOre[0,1,2],one);
        let b = add(rvOre[0,1,3],rvOre[2,0,-6]);
        let f = rvOreMul(d,b);
        assert rvOreEulerFactor(f) == d;
        assert rvOreEulerReduce(f) == b;
        assert rvOreEulerFactor(b) == one;
        assert rvOreMul(d,rvOreEulerReduce(f)) == f;
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn central_binomial_and_complementary_periods_use_the_scaled_recurrence() {
    let program = r"
        let f = rvOreFactorial[4](rvQGrAnnihilator[2,4]);
        let d = rvOreEulerFactor(f);
        let b = rvOreEulerReduce(f);
        assert d == rvOre[0,4];
        assert b == add(rvOre[0,1],add(rvOre[1,1,-4],rvOre[1,0,-2]));
        assert f == rvOreMul(d,b);
        rvQGrRegularizedRhs[2,4];
        rvQGrRegularized[2,4,8];
        rvQGrRegularized[3,5,6];
        rvQGrRegularizedRhs[3,5];
        rvQGrRegularized[1,1000000000,5];
        rvQGrRegularizedRhs[1000000000,1000000001];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvQGrRegularizedRhs = {}\n",
            "rvQGrRegularized = [1, 2, 6, 20, 70, 252, 924, 3432, 12870]\n",
            "rvQGrRegularized = [1, 3, 19, 147, 1251, 11253, 104959]\n",
            "rvQGrRegularizedRhs = {}\n",
            "rvQGrRegularized = [1, 1, 1, 1, 1, 1]\n",
            "rvQGrRegularizedRhs = {0: 1}\n",
        )
    );
    let rows = stdout(interpret("rvQGrRegularized[2,4,40];"));
    assert!(rows.ends_with(", 107507208733336176461620]\n"), "{rows}");
}

#[test]
fn regularization_checks_type_degree_storage_and_arithmetic() {
    let q_overflow = format!("show rvOreFactorial[3](rvOre[{},1]);", usize::MAX);
    let j_overflow = format!("show rvOreFactorial[1](rvOre[1,{}]);", usize::MAX);
    let j_sparse = format!("show rvOreFactorial[1](rvOre[0,{}]);", usize::MAX);
    assert_eq!(
        stdout(interpret(&j_sparse)),
        format!("RVOre{{(0, {}): 1}}\n", usize::MAX)
    );
    for (program, expected) in [
        ("show rvOreFactorial[2](S.hi);", "requires an Ore operator"),
        ("show rvOreEulerFactor(S.lo);", "requires an Ore operator"),
        (
            "show rvOreEulerReduce(rvOre[1,300001]);",
            "polynomial storage limit",
        ),
        (
            "show rvOreFactorial[1000000](rvOre[1,0]);",
            "polynomial storage limit",
        ),
        (
            "show rvOreFactorial[3](rvOre[200,0]);",
            "integer arithmetic is out of range",
        ),
        ("rvQGrRegularized[2,4,4097];", "period degree exceeds"),
        (
            "rvQGrRegularizedRhs[4,9];",
            "Schubert rank exceeds the limit",
        ),
        (q_overflow.as_str(), "integer arithmetic is out of range"),
        (j_overflow.as_str(), "integer arithmetic is out of range"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
