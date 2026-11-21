mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn quartic_models_have_exact_breaks_characters_and_weil_polynomials() {
    assert_eq!(
        stdout(example("rank-witt")),
        concat!(
            "rvWittBreaks[Small] = [1, 2]\n",
            "rvWittGenus[Small] = 1\n",
            "rvWittLPolynomial[Small] = [1, 2, 2]\n",
            "rvWittCharacterL[Small] = [(1, 0), (1, 1)]\n",
            "rvWittNewton[Small] = [1/2: 2]\n",
            "rvWittPoints[Small, 2] = 5\n",
            "rvWittPoints[Small, 4] = 25\n",
            "rvWittBreaks[CubicSecond] = [1, 3]\n",
            "rvWittGenus[CubicSecond] = 2\n",
            "rvWittLPolynomial[CubicSecond] = [1, 2, 2, 4, 4]\n",
            "rvWittCharacterL[CubicSecond] = [(1, 0), (1, -1), (0, -2)]\n",
            "rvWittNewton[CubicSecond] = [1/2: 4]\n",
            "rvWittPoints[CubicSecond, 6] = 65\n",
            "rvWittBreaks[CubicFirst] = [3, 6]\n",
            "rvWittGenus[CubicFirst] = 6\n",
            "rvWittLPolynomial[CubicFirst] = [1, 2, 2, 0, -2, -4, -8, -8, -8, 0, 32, 64, 64]\n",
            "rvWittCharacterL[CubicFirst] = [(1, 0), (1, 1), (-1, 1), (-2, 0), (0, -4), (4, -4)]\n",
            "rvWittNewton[CubicFirst] = [1/4: 4, 1/2: 4, 3/4: 4]\n",
            "rvWittPoints[CubicFirst, 8] = 209\n",
        )
    );
}

#[test]
fn the_actual_character_factors_multiply_to_the_point_count_polynomial() {
    for (m, g) in [
        (1, "zero[RVPoly([0],0)]"),
        (1, "rvmon[[0]]{3}"),
        (3, "zero[RVPoly([0],0)]"),
    ] {
        let program = format!(
            "rvwittcover W[[0]] {{rvmon[[0]]{{{m}}},{g}}};
             rvascover Quadratic[[0,1]] {{rvmon[[0,1]]{{0,{m}}},rvmon[[0,1]]{{{m},0}}}};
             rvWittLPolynomial[W]; rvWittCharacterL[W]; rvASLPolynomial[Quadratic];"
        );
        let result = stdout(interpret(&program));
        let mut lines = result.lines().map(|line| line.split_once(" = ").unwrap().1);
        let integral = parse_integral(lines.next().unwrap());
        let character: Vec<_> = lines
            .next()
            .unwrap()
            .trim_matches(['[', ']'])
            .split("), (")
            .map(|pair| {
                let (re, im) = pair.trim_matches(['(', ')']).split_once(", ").unwrap();
                (re.parse::<i128>().unwrap(), im.parse::<i128>().unwrap())
            })
            .collect();
        let conjugate: Vec<_> = character.iter().map(|&(re, im)| (re, -im)).collect();
        let quadratic = parse_integral(lines.next().unwrap());
        assert_eq!(
            convolve(&quadratic, &convolve(&character, &conjugate)),
            integral
        );
    }
}

fn parse_integral(value: &str) -> Vec<(i128, i128)> {
    value
        .trim_matches(['[', ']'])
        .split(", ")
        .map(|c| (c.parse().unwrap(), 0))
        .collect()
}

fn convolve(left: &[(i128, i128)], right: &[(i128, i128)]) -> Vec<(i128, i128)> {
    let mut result = vec![(0, 0); left.len() + right.len() - 1];
    for (i, &(a, b)) in left.iter().enumerate() {
        for (j, &(c, d)) in right.iter().enumerate() {
            result[i + j].0 += a * c - b * d;
            result[i + j].1 += a * d + b * c;
        }
    }
    result
}

#[test]
fn binary_homogenization_supports_actual_inhomogeneous_input() {
    let program = r"
        let f0 = add(rvmon[[0,1]]{3,0},add(rvmon[[0,1]]{2,1},rvmon[[0,1]]{0,3}));
        let f1 = add(rvmon[[0,1]]{0,3},add(rvmon[[0,1]]{1,2},rvmon[[0,1]]{3,0}));
        rvwittcover W0[[0]] { f0, rvmon[[]]{} };
        rvwittcover W1[[1]] { f1, rvmon[[1]]{0} };
        rvWittBreaks[W0]; rvWittBreaks[W1];
        rvWittLPolynomial[W0]; rvWittLPolynomial[W1];
        rvWittCharacterL[W0]; rvWittCharacterL[W1];
        rvWittPoints[W0,8]; rvWittPoints[W1,8];
        rvwittcover Twisted[[0]] { rvmon[[0]]{1},rvmon[[0]]{0} };
        rvWittLPolynomial[Twisted];
        rvWittCharacterL[Twisted];
    ";
    let result = stdout(interpret(program));
    let values: Vec<_> = result
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    for &[left, right] in values[..8].as_chunks::<2>().0 {
        assert_eq!(left, right);
    }
    assert_eq!(&values[8..], &["[1, -2, 2]", "[(1, 0), (-1, -1)]"]);
}

#[test]
fn witt_models_check_the_variable_pole_parity_and_budgets() {
    let setup =
        "let x = rvmon[[0]]{1}; let empty = zero[RVPoly([0],0)]; rvwittcover W[[0]] { x,empty };";
    for (statement, reason) in [
        (
            "rvwittcover K[[0,1]] {x,empty};",
            "exactly one parameter variable",
        ),
        ("rvwittcover K[[0]] {x};", "exactly two polynomials"),
        ("rvwittcover K[[1]] {x,empty};", "chosen variable"),
        (
            "rvwittcover K[[0]] {rvMul(x,x),empty};",
            "positive odd degree",
        ),
        ("rvwittcover K[[0]] {empty,x};", "positive odd degree"),
        (
            "rvwittcover K[[0]] {x,rvMul(x,x)};",
            "zero, constant, or of odd degree",
        ),
        ("rvwittcover W[[0]] {x,empty};", "already defined"),
        ("rvWittGenus[Missing];", "unknown Witt cover"),
        ("rvWittPoints[W,0];", "extension degree is out of range"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
    let error = diagnostic(interpret(&format!(
        "{setup}\nrvWittPoints[W,{}];",
        usize::BITS - 1
    )));
    assert!(error.contains("exceeds the computation limit"), "{error}");
    let error = diagnostic(interpret(&format!(
        "{setup}\nrvwittcover K[[0]] {{rvmon[[0]]{{{}}},empty}};",
        usize::MAX
    )));
    assert!(
        error.contains("upper ramification break is out of range"),
        "{error}"
    );
}
