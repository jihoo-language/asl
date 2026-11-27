mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn example_counts_noncommutative_products_and_the_unique_hn_flag() {
    assert_eq!(stdout(example("deligne-finite-hall")),concat!(
        "rvPartitionHallCoefficient[q=3] = 3\n",
        "rvPartitionHallCoefficient[q=3] = 1\n",
        "rvPartitionHallFlag[q=3] = 4\n",
        "rvPartitionHallTerm = 2\n",
        "rvPartitionHallHN = [(3,[(2, 1, 2)]), (-2,[(0, 1, 1)])]\n",
        "rvPartitionHallHNCheck[q=3] = flags=1, factors=[(3,[(2, 1, 2)]), (-2,[(0, 1, 1)])]\n",
        "rvPartitionHallChamber[2,0,2] = theta=[-1, 1, 2], positiveQuotients=[[1, 1], [2, 1], [2, 2], [2, 3]]\n",
        "rvPartitionHallWall[2,0,3,2] = theta=[-1, 1, 1], kingStableBrick=[2, 2]\n",
    ));
}

#[test]
fn hall_counts_over_prime_fields_and_actual_named_residue_fields() {
    let mut input = String::from("rvunram F[2,2]{1,1,1};\nrvunram E[3,1]{1,0,1};\n");
    let mut expected = String::new();
    for (j, (field, q)) in [("2", 2), ("3", 3), ("5", 5), ("F", 4), ("E", 9)]
        .into_iter()
        .enumerate()
    {
        writeln!(
            input,
            "let s{j}=rvPartitionHallModule[1,0,{field}]{{(0,1,1)}};"
        )
        .unwrap();
        writeln!(
            input,
            "let s2{j}=rvPartitionHallModule[1,0,{field}]{{(0,1,2)}};"
        )
        .unwrap();
        writeln!(
            input,
            "let j2{j}=rvPartitionHallModule[1,0,{field}]{{(0,2,1)}};"
        )
        .unwrap();
        writeln!(
            input,
            "let a{j}=rvPartitionHallBasis(s{j});\nlet z{j}=rvPartitionHallProduct(a{j},a{j});"
        )
        .unwrap();
        writeln!(input,"rvPartitionHallCoefficient(s2{j},s{j},s{j});\nrvPartitionHallCoefficient(j2{j},s{j},s{j});").unwrap();
        writeln!(
            input,
            "rvPartitionHallTerm(z{j},s2{j});\nrvPartitionHallTerm(z{j},j2{j});"
        )
        .unwrap();
        writeln!(expected, "rvPartitionHallCoefficient[q={q}] = {}", q + 1).unwrap();
        writeln!(expected, "rvPartitionHallCoefficient[q={q}] = 1").unwrap();
        writeln!(
            expected,
            "rvPartitionHallTerm = {}\nrvPartitionHallTerm = 1",
            q + 1
        )
        .unwrap();
        writeln!(
            input,
            "assert rvPartitionHallProduct(z{j},a{j}) == rvPartitionHallProduct(a{j},z{j});"
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&input)), expected);
}

#[test]
fn noncommuting_products_and_flag_associativity_use_every_extension_type() {
    let mut source = String::new();
    let mut expected = String::new();
    for q in [2, 3, 5] {
        writeln!(
            source,
            "let A{q}=rvPartitionHallModule[2,0,{q}]{{(0,1,1)}};"
        )
        .unwrap();
        writeln!(
            source,
            "let B{q}=rvPartitionHallModule[2,0,{q}]{{(2,1,1)}};"
        )
        .unwrap();
        writeln!(
            source,
            "let V{q}=rvPartitionHallModule[2,0,{q}]{{(0,2,1)}};"
        )
        .unwrap();
        writeln!(
            source,
            "let L{q}=rvPartitionHallModule[2,0,{q}]{{(0,2,1),(2,1,1)}};"
        )
        .unwrap();
        writeln!(
            source,
            "let W{q}=rvPartitionHallModule[2,0,{q}]{{(2,2,1)}};"
        )
        .unwrap();
        writeln!(source,"let a{q}=rvPartitionHallBasis(A{q});\nlet b{q}=rvPartitionHallBasis(B{q});\nlet v{q}=rvPartitionHallBasis(V{q});").unwrap();
        writeln!(source,"rvPartitionHallCoefficient(L{q},V{q},B{q});\nrvPartitionHallCoefficient(L{q},B{q},V{q});").unwrap();
        writeln!(source,"rvPartitionHallCoefficient(V{q},A{q},B{q});\nrvPartitionHallCoefficient(W{q},A{q},B{q});").unwrap();
        for (a, b, c) in [("a", "b", "a"), ("v", "b", "a"), ("b", "a", "v")] {
            writeln!(source,"assert rvPartitionHallProduct(rvPartitionHallProduct({a}{q},{b}{q}),{c}{q}) == rvPartitionHallProduct({a}{q},rvPartitionHallProduct({b}{q},{c}{q}));").unwrap();
        }
        for k in [q, 1, 1, 0] {
            writeln!(expected, "rvPartitionHallCoefficient[q={q}] = {k}").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn king_stability_and_rational_hn_slopes_agree_with_the_actual_finite_flags() {
    let output = stdout(interpret(concat!(
        "let A=rvPartitionHallModule[2,0,2]{(0,3,1)};\n",
        "let B=rvPartitionHallModule[2,0,2]{(2,3,1)};\n",
        "let C=rvPartitionHallModule[2,0,2]{(2,2,1)};\n",
        "let D=rvPartitionHallModule[2,0,2]{(2,2,2)};\n",
        "rvPartitionHallHN[-2,0,3](A);\n",
        "rvPartitionHallHNCheck[-2,0,3](A);\n",
        "rvPartitionHallHNCheck[-2,0,3](B);\n",
        "rvPartitionHallStability[-1,1,1](C);\n",
        "rvPartitionHallStability[-1,1,1](D);\n",
        "rvPartitionHallStability[0,0,0](C);\n",
        "rvPartitionHallHNCheck[-1,1,1](D);\n",
        "rvPartitionHallHNCheck[1,0,-1](zero[RVPartitionHallModule(2,0,2)]);\n",
    )));
    assert_eq!(output,concat!(
        "rvPartitionHallHN = [(1/2,[(2, 2, 1)]), (-2,[(0, 1, 1)])]\n",
        "rvPartitionHallHNCheck[q=2] = flags=1, factors=[(1/2,[(2, 2, 1)]), (-2,[(0, 1, 1)])]\n",
        "rvPartitionHallHNCheck[q=2] = flags=1, factors=[(3,[(2, 1, 1)]), (1/2,[(2, 2, 1)])]\n",
        "rvPartitionHallStability = slope=0, semistable=true, stable=true, kingSemistable=true, kingStable=true\n",
        "rvPartitionHallStability = slope=0, semistable=true, stable=false, kingSemistable=true, kingStable=false\n",
        "rvPartitionHallStability = slope=0, semistable=true, stable=false, kingSemistable=true, kingStable=false\n",
        "rvPartitionHallHNCheck[q=2] = flags=1, factors=[(0,[(2, 2, 2)])]\n",
        "rvPartitionHallHNCheck[q=2] = flags=1, factors=[]\n",
    ));
}

#[test]
fn all_tau_chamber_interiors_recover_fac_and_all_walls_recover_stable_bricks() {
    let mut source = String::new();
    let mut expected_walls = 0;
    for (n, k) in [(0, 12), (1, 24), (2, 24), (3, 16)] {
        for i in 0..k {
            writeln!(source, "rvPartitionHallChamber[2,{n},{i}];").unwrap();
        }
        let walls = stdout(interpret(&format!("rvPartitionTauWalls[2,{n}];")));
        for edge in walls.split("cover=(").skip(1) {
            let (pair, _) = edge.split_once(')').unwrap();
            writeln!(source, "rvPartitionHallWall[2,{n},{pair}];").unwrap();
            expected_walls += 1;
        }
    }
    let output = stdout(interpret(&source));
    assert_eq!(
        output.matches("positiveQuotients=").count(),
        12 + 24 + 24 + 16
    );
    assert_eq!(output.matches("kingStableBrick=").count(), expected_walls);
}

#[test]
fn unit_and_coefficient_zero_are_exact_and_work_budgets_do_not_truncate() {
    let output=stdout(interpret(concat!(
        "let M=rvPartitionHallModule[1,0,2]{(0,1,10),(0,1,2)};\n",
        "let Z=zero[RVPartitionHallModule(1,0,2)];\n",
        "let a=rvPartitionHallBasis(M);\n",
        "let u=rvPartitionHallUnit[1,0,2];\n",
        "assert rvPartitionHallBasis(Z)==u;\n",
        "assert rvPartitionHallProduct(u,a)==a;\n",
        "assert rvPartitionHallProduct(a,zero[RVPartitionHall(1,0,2)])==zero[RVPartitionHall(1,0,2)];\n",
        "assert add(a,neg(a))==zero[RVPartitionHall(1,0,2)];\n",
        "rvPartitionHallCoefficient(M,M,Z);\n",
        "rvPartitionHallCoefficient(Z,M,M);\n",
        "rvPartitionHallTerm(a,Z);\n",
    )));
    assert_eq!(output,"rvPartitionHallCoefficient[q=2] = 1\nrvPartitionHallCoefficient[q=2] = 0\nrvPartitionHallTerm = 0\n");
    for (source, expected) in [
        ("show rvPartitionHallModule[1,0,4]{};","must be prime"),
        ("show rvPartitionHallUnit[3,0,2];","limited to rank 2"),
        ("show rvPartitionHallModule[1,0,2]{(0,1,65)};","length exceeds the limit of 64"),
        ("show rvPartitionHallModule[1,0,2]{(0,3,1)};","top or length is out of range"),
        ("show rvPartitionHallProduct(rvPartitionHallUnit[1,0,2],rvPartitionHallUnit[1,0,3]);","same rank, parameter, and coefficient field"),
        ("rvunram F[2,1]{1,1,1}; rvunram E[2,1]{1,1,1}; show rvPartitionHallProduct(rvPartitionHallUnit[1,0,F],rvPartitionHallUnit[1,0,E]);","same rank, parameter, and coefficient field"),
        ("rvPartitionHallHN[0,1](rvPartitionHallModule[1,0,2]{(0,1,1)});","one entry per simple vertex"),
        ("let a=rvPartitionHallBasis(rvPartitionHallModule[1,0,2]{(0,1,4)}); show rvPartitionHallProduct(a,a);","limited to total length 6"),
        ("let S=rvPartitionHallModule[1,0,65537]{(0,1,1)}; let SS=rvPartitionHallModule[1,0,65537]{(0,1,2)}; rvPartitionHallCoefficient(SS,S,S);","exceed 50,000"),
    ] {
        let message=diagnostic(interpret(source));
        assert!(message.contains(expected),"{source}: {message}");
    }
}
