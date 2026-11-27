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

#[test]
fn euler_gcd_certifies_the_small_factor_when_pseudo_remainders_overflow() {
    let program = r"
        let first = add(rvOre[1,4],rvOre[1,0]);
        let transformed = rvOreFactorial[4](add(first,rvOre[6,0]));
        let factor = rvOre[0,4];
        let reduced = add(first,rvOreMul(rvOre[1,0],rvOreFactorial[4](rvOre[5,0])));
        assert rvOreMul(factor,reduced) == transformed;
        assert rvOreEulerFactor(transformed) == factor;
        assert rvOreEulerReduce(transformed) == reduced;

        let left = add(add(rvOre[0,3,-3],rvOre[0,2,-2]),add(rvOre[0,0,-2],add(rvOre[2,2,3],rvOre[2,1,-1])));
        let right = add(rvOre[2,3,-2],add(rvOre[2,2,1],rvOre[2,1,-1]));
        let second = rvOreFactorial[3](rvOreMul(left,right));
        let expected = add(rvOre[0,7],add(rvOre[0,6,-5],add(rvOre[0,5,9],add(rvOre[0,4,-7],rvOre[0,3,2]))));
        assert rvOreEulerFactor(second) == expected;
        assert rvOreMul(expected,rvOreEulerReduce(second)) == second;
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn euler_reduction_checks_actual_quotient_range_after_exact_cancellation() {
    let high = i128::MAX;
    let cancellation = format!(
        "let divisor = add(rvOre[0,0],add(rvOre[0,1],rvOre[0,2]));
         let a = add(rvOre[0,0,{high}],add(rvOre[0,2,{}],add(rvOre[0,3,{}],rvOre[0,5,{high}])));
         let shifted_divisor = add(rvOre[1,0,3],add(rvOre[1,1,3],rvOre[1,2]));
         let input = add(a,shifted_divisor);
         let q = add(rvOre[0,0,{high}],add(rvOre[0,1,{}],add(rvOre[0,2,{}],rvOre[0,3,{high}])));
         assert rvOreEulerFactor(input) == divisor;
         assert rvOreEulerReduce(input) == add(q,rvOre[1,0]);",
        -high, -high, -high, -high,
    );
    assert_eq!(stdout(interpret(&cancellation)), "");

    let out_of_range = format!(
        "let a = add(rvOre[0,0,-1],add(rvOre[0,1,{high}],add(rvOre[0,2,{high}],rvOre[0,3,-1])));
         let input = add(a,add(rvOre[1,0,2],rvOre[1,1]));
         assert rvOreEulerFactor(input) == add(rvOre[0,0],rvOre[0,1]);
         show rvOreEulerReduce(input);",
    );
    let error = diagnostic(interpret(&out_of_range));
    assert!(
        error.contains("integer arithmetic is out of range"),
        "{error}"
    );
}

#[test]
fn oversized_shifted_blocks_can_have_a_small_euler_factor() {
    let program = r"
        let one = rvOre[0,0];
        let t = rvOre[0,1];
        let large = rvOre[8,50];
        let coprime = add(t,large);
        assert rvOreEulerFactor(coprime) == one;
        assert rvOreEulerReduce(coprime) == coprime;

        let linear = add(t,rvOre[0,0,-8]);
        let input = add(linear,large);
        assert rvOreEulerFactor(input) == linear;
        assert rvOreEulerReduce(input) == add(one,rvOre[8,49]);

        let input2 = add(rvOre[0,2],add(rvOre[8,51],rvOre[8,50,8]));
        assert rvOreEulerFactor(input2) == t;
        assert rvOreEulerReduce(input2) == add(t,large);

        let constant = add(rvOre[8,0,-7],rvOre[19,70,-3]);
        assert rvOreEulerFactor(constant) == one;
        assert rvOreEulerReduce(constant) == constant;
        assert rvOreEulerFactor(zero[RVOre]) == one;
    ";
    assert_eq!(stdout(interpret(program)), "");
    let huge = format!(
        "let x = add(rvOre[0,1],rvOre[{},20,-2]);
         assert rvOreEulerFactor(x) == rvOre[0,0];
         assert rvOreEulerReduce(x) == x;",
        usize::MAX
    );
    assert_eq!(stdout(interpret(&huge)), "");
}

#[test]
fn shifted_cofactors_are_separately_certified_in_original_block_coordinates() {
    let m = 1i128 << 126;
    let program = format!(
        "let f = add(rvOre[0,1,{m}],rvOre[0,0,{}]);
         let input = add(f,add(rvOre[2,2,{}],add(rvOre[2,1,-1],rvOre[2,0,{}])));
         assert rvOreEulerFactor(input) == f;
         assert rvOreEulerReduce(input) == add(rvOre[0,0],add(rvOre[2,0],rvOre[2,1,-1]));
         let single = add(rvOre[8,50,-6],rvOre[8,0,12]);
         assert rvOreEulerReduce(single) == rvOre[8,0,-6];
         assert rvOreEulerReduce(rvOre[8,50,{}]) == rvOre[8,0,{}];
         let shifted = rvOreMul(rvOre[8,0],rvOreFactorial[50](rvOre[1,0]));
         assert rvOreEulerReduce(add(rvOre[8,50],shifted)) == add(rvOre[8,0],rvOre[9,0]);",
        1 - m,
        -m,
        m + 1,
        i128::MIN,
        i128::MIN,
    );
    assert_eq!(stdout(interpret(&program)), "");

    let unrepresentable = diagnostic(interpret("show rvOreEulerFactor(rvOre[8,50]);"));
    assert!(
        unrepresentable.contains("integer arithmetic is out of range"),
        "{unrepresentable}"
    );
}
