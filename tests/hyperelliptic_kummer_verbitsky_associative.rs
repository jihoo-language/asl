mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARK: &str = "rvhodgewick B[20]{0,1}{[0,1],[1,0]};";
const EXAMPLE: &str = include_str!("../examples/hyperelliptic-kummer-verbitsky-associative.asl");

#[test]
fn noncommutative_mixed_obstructions_are_two_independent_actual_classes() {
    let text = stdout(example("hyperelliptic-kummer-verbitsky-associative"));
    let lines = text.lines().collect::<Vec<_>>();
    assert!(lines[0].ends_with("B=10, HH=4)"));
    assert_eq!(
        lines[1],
        "rvVerbitskyAssocReport[Mix] = weights=[-2, -2, 0] requested=3 associative_through=1 corrections=[] obstructions=[([0, 1, 1], w=-2, [([2, 2, 1], 1, -2), ([2, 2, 3], 3, -2)]), ([1, 0, 1], w=-2, [([2, 1, 1], 2, 2), ([3, 1, 1], 3, 2)])]"
    );
    // The disjoint coefficient supports explicitly witness both classes.
    assert!(lines[2].contains("order=2 obstructed=true classes=[(-2, [([2, 1, 1], 2, 2)"));
    assert!(lines[3].contains("order=2 obstructed=true classes=[(-2, [([2, 2, 1], 1, -2)"));
    assert!(lines[4].ends_with("order=2 obstructed=false classes=[]"));
    assert_eq!(
        &lines[5..7],
        ["rvVerbitskyHHCoeff = 1", "rvVerbitskyHHCoeff = -1"]
    );
    assert_eq!(
        lines[7],
        "rvVerbitskyAssocProduct[Mix, 1, 1] = [([1, 0, 0], [(0, 1)])]"
    );
    assert_eq!(
        lines[8],
        "rvVerbitskyAssocProduct[Mix, 2, 1] = [([0, 0, 0], [(3, 1)]), ([0, 0, 1], [(3, 1)])]"
    );
    assert_eq!(
        lines[9],
        "rvVerbitskyAssocProduct[Mix, 1, 2] = [([0, 0, 0], [(3, 1)])]"
    );
}

#[test]
fn pure_jet_products_match_all_pbw_coefficients_including_units() {
    let mut input = EXAMPLE.to_string();
    for (a, b) in (0..4).flat_map(|a| (0..4).map(move |b| (a, b))) {
        writeln!(input, "rvVerbitskyAssocProduct[Carry,{a},{b}];").unwrap();
        writeln!(input, "rvVerbitskyAssocProduct[Quantum,{a},{b}];").unwrap();
    }
    let output = stdout(interpret(&input));
    let lines = output.lines().collect::<Vec<_>>();
    assert!(lines.iter().any(
        |line| line.contains("associative_through=4 corrections=[([1, 1], [([3, 3], 0, 1)])]")
    ));
    let mut actual = lines[lines.len() - 32..].iter();
    for (a, b) in (0..4).flat_map(|a| (0..4).map(move |b| (a, b))) {
        let u = a % 2 + b % 2;
        let v = a / 2 + b / 2;
        let target = u % 2 + 2 * (v % 2);
        let carry = format!(
            "rvVerbitskyAssocProduct[Carry, {a}, {b}] = [([{}, {}], [({target}, 1)])]",
            u / 2,
            v / 2
        );
        assert_eq!(**actual.next().unwrap(), carry);
        let product = if u > 1 || v > 1 {
            "[]".to_owned()
        } else if a / 2 == 1 && b % 2 == 1 {
            format!("[([0], [({target}, 1)]), ([1], [({target}, 1)])]")
        } else {
            format!("[([0], [({target}, 1)])]")
        };
        assert_eq!(
            **actual.next().unwrap(),
            format!("rvVerbitskyAssocProduct[Quantum, {a}, {b}] = {product}")
        );
    }
}

#[test]
fn curved_exact_tangent_requires_cubic_and_quartic_primitives() {
    let output = stdout(interpret(&format!(
        "{MARK}
        let f=add(add(rvVerbitskyHHTerm[B,1,1,-1,0]{{1}},rvVerbitskyHHTerm[B,1,1,-1,0]{{2}}),
                  add(rvVerbitskyHHTerm[B,1,1,-1,1]{{3}},rvVerbitskyHHTerm[B,1,1,-1,2]{{3}}));
        let d=rvVerbitskyHHDiff(f);
        assert rvVerbitskyHHClass(d)==zero[RVVerbitskyHH(B,1,2,-1)];
        rvverbitskyassoc G[5]{{d}};
        let m3=rvVerbitskyAssocCoefficient[G]{{3}};
        let m4=rvVerbitskyAssocCoefficient[G]{{4}};
        assert rvVerbitskyHHDiff(m3)==rvVerbitskyAssocResidual[G]{{3}};
        assert rvVerbitskyHHDiff(m4)==rvVerbitskyAssocResidual[G]{{4}};
        assert rvVerbitskyAssocObstruction[G]{{3}}==zero[RVVerbitskyHH(B,1,3,-3)];
        assert m4==rvVerbitskyHHTerm[B,1,2,-4,0,4]{{3,3}};
        rvVerbitskyHHCoeff[0]{{1,3}}(m3);
        rvVerbitskyAssocReport[G];
        rvVerbitskyAssocProduct[G,3,3];
        "
    )));
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "rvVerbitskyHHCoeff = 4");
    assert!(lines[1].contains("requested=5 associative_through=5"));
    assert!(lines[1].contains("([4], [([3, 3], 0, 4)])"));
    assert_eq!(
        lines[2],
        "rvVerbitskyAssocProduct[G, 3, 3] = [([4], [(0, 4)])]"
    );
}

#[test]
fn rank_one_carry_preserves_the_flat_cubic_product() {
    let output = stdout(interpret(
        r"
        rvhodgewick B[21]{0,1}{[2]};
        let t=rvVerbitskyHHBasis[B,1,2,-3,0];
        rvverbitskyassoc J[6]{t};
        rvVerbitskyAssocReport[J];
        rvVerbitskyAssocProduct[J,0,1];
        rvVerbitskyAssocProduct[J,1,2];
        rvVerbitskyAssocProduct[J,2,2];
        assert rvVerbitskyAssocCoefficient[J]{2}==zero[RVVerbitskyHochschildCochain(B,1,2,-6)];
        ",
    ));
    assert_eq!(output, "rvVerbitskyAssocReport[J] = weights=[-3] requested=6 associative_through=6 corrections=[] obstructions=[]\nrvVerbitskyAssocProduct[J, 0, 1] = [([0], [(1, 1)])]\nrvVerbitskyAssocProduct[J, 1, 2] = [([1], [(0, 1)])]\nrvVerbitskyAssocProduct[J, 2, 2] = [([1], [(1, 1)])]\n");
}

#[test]
fn rejects_bad_tangents_uncertified_higher_coefficients_and_unbounded_orders() {
    for (query, reason) in [
        ("rvverbitskyassoc J[2]{rvVerbitskyHHTerm[B,1,2,-2,0]{1,1}};", "not a Hochschild 2-cocycle"),
        ("rvverbitskyassoc J[2]{rvVerbitskyHHBasis[B,1,0,0,0]};", "requires Hochschild 2-cochains"),
        ("rvverbitskyassoc J[2]{};", "requires at least one tangent"),
        ("rvverbitskyassoc J[0]{zero[RVVerbitskyHH(B,1,2,0)]};", "positive order"),
        ("rvverbitskyassoc J[1000000000000]{zero[RVVerbitskyHH(B,1,2,0)]};", "work or storage limit"),
        ("rvVerbitskyAssocReport[Missing];", "unknown Verbitsky associative jet"),
        ("let k=rvVerbitskyAssocCoefficient[Mix]{1,0,1};", "outside the certified positive orders"),
        ("let k=rvVerbitskyAssocObstruction[Mix]{1,1,1};", "has yet to be tested"),
        ("let k=rvVerbitskyAssocResidual[Mix]{2,0};", "wrong number of exponents"),
        ("rvVerbitskyAssocDirection[Mix]{1,0};", "wrong number of parameters"),
        ("rvVerbitskyAssocProduct[Mix,5,0];", "outside the ring basis"),
        ("rvverbitskyassoc Mix[2]{p};", "already defined"),
        ("rvhodgewick X[20]{0,1}{[0,1],[1,0]}; rvverbitskyassoc X[2]{p,zero[RVVerbitskyHH(X,1,2,0)]};", "same marking and Hilbert length"),
    ] {
        assert!(diagnostic(interpret(&format!("{EXAMPLE} {query}"))).contains(reason), "{query}");
    }
}
