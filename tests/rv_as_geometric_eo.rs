mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_same_zeta_polynomial_and_a_number_can_have_different_geometric_eo_types() {
    let output = stdout(example("rank-as-geometric-eo"));
    assert_eq!(
        output,
        concat!(
            "rvASLPolynomial[A] = [1, -1, 2, -2, 4, -4, 8, -8, 16, -16, 32]\n",
            "rvASLPolynomial[B] = [1, -1, 2, -2, 4, -4, 8, -8, 16, -16, 32]\n",
            "rvASPrank[A] = 1\n",
            "rvASPrank[B] = 1\n",
            "rvASaNumber[A] = 2\n",
            "rvASaNumber[B] = 2\n",
            "rvASPoleOrders[A] = [(1, 1), (9, 1)]\n",
            "rvASPoleOrders[B] = [(5, 2)]\n",
            "rvASEO[A] = [1, 1, 2, 2, 3]\n",
            "rvASEO[B] = [1, 1, 1, 2, 3]\n",
        )
    );
    let mut source = include_str!("../examples/rank-as-geometric-eo.asl").to_owned();
    let points = [2, 8, 8, 24, 32, 80, 128, 288, 512, 1088];
    let mut expected = output;
    for (i, points) in points.into_iter().enumerate() {
        let m = i + 1;
        for curve in ["A", "B"] {
            writeln!(source, "rvASPoints[{curve},{m}];").unwrap();
            writeln!(expected, "rvASPoints[{curve}, {m}] = {points}").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn the_multiplicity_digits_count_all_geometric_roots_in_nonrational_orbits() {
    let input = r"
        let s = rvmon[[0,1]]{1,0};
        let t = rvmon[[0,1]]{0,1};
        let q2 = add(rvMul(s,s),add(rvMul(s,t),rvMul(t,t)));
        let q3 = add(rvmon[[0,1]]{3,0},add(rvmon[[0,1]]{2,1},rvmon[[0,1]]{0,3}));
        let cubes = rvMul(rvMul(q2,q2),q2);
        let mixed = rvMul(rvMul(cubes,q3),rvmon[[0,1]]{0,5});
        rvascover Orbit[[0,1]] {rvmon[[0,1]]{6,0},cubes};
        rvascover Mixed[[0,1]] {rvmon[[0,1]]{14,0},mixed};
        rvascover FiniteOnly[[0,1]] {rvmon[[0,1]]{5,0},rvmon[[0,1]]{0,5}};
        rvascover InfinityOnly[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
        rvASPoleOrders[Orbit];
        rvASGenus[Orbit];
        rvASPrank[Orbit];
        rvASEO[Orbit];
        rvASPoleOrders[Mixed];
        rvASGenus[Mixed];
        rvASPrank[Mixed];
        rvASEO[Mixed];
        rvASPoleOrders[FiniteOnly]; rvASEO[FiniteOnly];
        rvASPoleOrders[InfinityOnly]; rvASEO[InfinityOnly];
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvASPoleOrders[Orbit] = [(3, 2)]\n",
            "rvASGenus[Orbit] = 3\n",
            "rvASPrank[Orbit] = 1\n",
            "rvASEO[Orbit] = [1, 1, 1]\n",
            "rvASPoleOrders[Mixed] = [(1, 3), (3, 2), (5, 1)]\n",
            "rvASGenus[Mixed] = 9\n",
            "rvASPrank[Mixed] = 5\n",
            "rvASEO[Mixed] = [1, 2, 3, 4, 5, 5, 5, 5, 6]\n",
            "rvASPoleOrders[FiniteOnly] = [(5, 1)]\n",
            "rvASEO[FiniteOnly] = [0, 1]\n",
            "rvASPoleOrders[InfinityOnly] = [(5, 1)]\n",
            "rvASEO[InfinityOnly] = [0, 1]\n",
        )
    );
}

#[test]
fn pole_order_reports_stay_sparse_before_the_checked_matrix_bound() {
    let order = (1usize << (usize::BITS / 2)) + 1;
    let setup =
        format!("rvascover Huge[[0,1]] {{rvmon[[0,1]]{{{order},0}},rvmon[[0,1]]{{0,{order}}}}};");
    assert_eq!(
        stdout(interpret(&format!("{setup} rvASPoleOrders[Huge];"))),
        format!("rvASPoleOrders[Huge] = [({order}, 1)]\n")
    );
    let error = diagnostic(interpret(&format!("{setup} rvASEO[Huge];")));
    assert!(error.contains("storage limit"), "{error}");
    let error = diagnostic(interpret(&format!("{setup} show zero[ASDeRham(Huge,3)];")));
    assert!(
        error.contains("requires a one-pole polynomial cover"),
        "{error}"
    );
}
