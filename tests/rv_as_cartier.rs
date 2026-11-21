mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_cartier_rank_ladders_and_a_numbers() {
    assert_eq!(
        stdout(example("rank-cartier")),
        concat!(
            "ASDiff(GenusTwo, 11){6*omega0}\n",
            "rvASCartierRank[GenusTwo, 0] = 2\n",
            "rvASCartierRank[GenusTwo, 1] = 1\n",
            "rvASCartierRank[GenusTwo, 2] = 0\n",
            "rvASaNumber[GenusTwo] = 1\n",
            "rvASCartierRank[GenusFour, 0] = 4\n",
            "rvASCartierRank[GenusFour, 1] = 2\n",
            "rvASCartierRank[GenusFour, 2] = 1\n",
            "rvASCartierRank[GenusFour, 3] = 0\n",
            "rvASaNumber[GenusFour] = 2\n",
            "rvASCartierRank[Ordinary, 1000] = 1\n",
            "rvASaNumber[Ordinary] = 0\n",
        )
    );
}

#[test]
fn extension_field_differentials_are_first_class_and_cartier_is_additive() {
    let program = r"
        rvascover C[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
        let a = rvASOmega[C,11,0];
        let b = rvASOmega[C,11,1];
        let x = add(a,rvASDiffScale[2](b));
        let y = rvASDiffScale[3](b);
        let c = \w : ASDiff(C,11) => rvASCartier(w);
        assert c(x) == rvASDiffScale[6](a);
        assert c(y) == rvASDiffScale[7](a);
        assert c(add(x,y)) == add(c(x),c(y));
        assert c(add(x,y)) == a;
        assert c(c(x)) == zero[ASDiff(C,11)];
        assert rvASDiffScale[0](x) == zero[ASDiff(C,11)];
        assert add(x,x) == zero[ASDiff(C,11)];
        show x;
        card (\w : ASDiff(C,11) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "ASDiff(C, 11){omega0 + 2*omega1}\ncard = 64 / 64\n"
    );
}

#[test]
fn stable_cartier_ranks_match_p_rank_and_the_l_polynomial_modulo_two() {
    let program = r"
        rvascover Poly[[0,1]] {rvmon[[0,1]]{0,9},rvmon[[0,1]]{9,0}};
        rvascover Swapped[[0,1]] {rvmon[[0,1]]{9,0},rvmon[[0,1]]{0,9}};
        rvascover Mixed[[0,1]] {
            rvmon[[0,1]]{5,0},
            add(rvmon[[0,1]]{0,5},add(rvmon[[0,1]]{1,4},rvmon[[0,1]]{2,3}))
        };
        rvascover Quadratic[[0,1]] {
            rvmon[[0,1]]{2,0},
            add(rvmon[[0,1]]{2,0},add(rvmon[[0,1]]{1,1},rvmon[[0,1]]{0,2}))
        };
    ";
    for (name, genus, p_rank) in [
        ("Poly", 4, 0),
        ("Swapped", 4, 0),
        ("Mixed", 3, 2),
        ("Quadratic", 1, 1),
    ] {
        let mut source = format!("{program}\nrvASPrank[{name}];\nrvASCartierRank[{name},1000];\nrvASLPolynomial[{name}];");
        for index in 0..genus {
            writeln!(source, "show rvASCartier(rvASOmega[{name},3,{index}]);").unwrap();
        }
        let rendered = stdout(interpret(&source));
        let mut lines = rendered.lines();
        assert_eq!(
            lines.next().unwrap(),
            format!("rvASPrank[{name}] = {p_rank}")
        );
        assert_eq!(
            lines.next().unwrap(),
            format!("rvASCartierRank[{name}, 1000] = {p_rank}")
        );
        let (_, l) = lines.next().unwrap().split_once(" = ").unwrap();
        let parity = l
            .trim_matches(['[', ']'])
            .split(", ")
            .enumerate()
            .map(|(degree, coefficient)| {
                ((coefficient.parse::<i128>().unwrap() & 1) as usize) << degree
            })
            .fold(0, |left, right| left | right);
        let columns: Vec<_> = lines
            .map(|line| {
                let (_, coordinates) = line.strip_suffix('}').unwrap().split_once('{').unwrap();
                if coordinates == "0" {
                    return 0;
                }
                coordinates
                    .split(" + ")
                    .map(|basis| {
                        1usize
                            << basis
                                .strip_prefix("omega")
                                .unwrap()
                                .parse::<usize>()
                                .unwrap()
                    })
                    .fold(0, |left, right| left | right)
            })
            .collect();
        assert_eq!(columns.len(), genus);
        assert_eq!(determinant(&columns, 0, 0), parity, "{name}");
    }
    let swapped = stdout(interpret(&format!(
        "{program}\nrvASCartierRank[Swapped,1]; rvASCartierRank[Swapped,2];"
    )));
    assert_eq!(
        swapped,
        "rvASCartierRank[Swapped, 1] = 2\nrvASCartierRank[Swapped, 2] = 1\n"
    );
}

// Enumerate determinant terms of (I + T*C) independently in F2[T].
fn determinant(columns: &[usize], row: usize, used: usize) -> usize {
    if row == columns.len() {
        return 1;
    }
    let mut result = 0;
    for (column, &coefficients) in columns.iter().enumerate() {
        if used >> column & 1 != 0 {
            continue;
        }
        let tail = determinant(columns, row + 1, used | (1 << column));
        if row == column {
            result ^= tail;
        }
        if coefficients >> row & 1 != 0 {
            result ^= tail << 1;
        }
    }
    result
}

#[test]
fn huge_monomial_poles_do_not_force_a_dense_cartier_matrix() {
    let h = 1usize << (usize::BITS / 2);
    let genus = h / 2;
    let program = format!(
        "rvascover Large[[0,1]] {{rvmon[[0,1]]{{0,{}}},rvmon[[0,1]]{{{},0}}}};
         rvASCartierRank[Large,1]; rvASCartierRank[Large,{}]; rvASaNumber[Large];
         show rvASCartier(rvASOmega[Large,3,{}]);",
        h + 1,
        h + 1,
        usize::MAX,
        genus - 1,
    );
    assert_eq!(stdout(interpret(&program)), format!(
        "rvASCartierRank[Large, 1] = {}\nrvASCartierRank[Large, {}] = 0\nrvASaNumber[Large] = {}\nASDiff(Large, 3){{omega{}}}\n",
        genus / 2, usize::MAX, genus / 2, (genus - 2) / 2,
    ));
}

#[test]
fn differentials_check_genus_field_coefficients_and_type() {
    let setup = "rvascover Line[[0,1]] {rvmon[[0,1]]{0,1},rvmon[[0,1]]{1,0}};";
    assert_eq!(
        stdout(interpret(&format!(
            "{setup}\nshow rvASCartier(zero[ASDiff(Line,11)]); rvASaNumber[Line];"
        ))),
        "ASDiff(Line, 11){0}\nrvASaNumber[Line] = 0\n"
    );
    for (input, reason) in [
        ("show rvASOmega[Line,3,0];", "out of range for genus 0"),
        ("show rvASCartier(rvmon[[0,1]]{0,0});", "requires an ASDiff"),
        (
            "show rvASDiffScale[8](zero[ASDiff(Line,11)]);",
            "field coefficient 8 is out of range",
        ),
        ("show zero[ASDiff(Line,5)];", "reducible over F2"),
        (
            "show zero[ASDiff(Missing,3)];",
            "unknown Artin–Schreier cover",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{input}")));
        assert!(error.contains(reason), "{input}\n{error}");
    }
}
