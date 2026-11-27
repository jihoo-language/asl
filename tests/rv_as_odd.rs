mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn prime_five_and_seven_character_polygons_and_signed_normalization() {
    assert_eq!(
        stdout(example("rank-as-odd")),
        concat!(
            "rvASOddReduced[Rewrite] = [(3, 1), (6, 1)]\n",
            "rvASOddGenus[Cubic] = 10\n",
            "rvASOddPrank[Cubic] = 0\n",
            "rvASOddCharacterNewton[Pure, 1] = [1/2: 5]\n",
            "rvASOddCharacterNewton[Cubic, 1] = [3/8: 2, 1/2: 1, 5/8: 2]\n",
            "rvASOddCharacterNewton[Quartic, 1] = [1/4: 1, 1/2: 3, 3/4: 1]\n",
            "rvASOddNewton[Cubic] = [3/8: 8, 1/2: 4, 5/8: 8]\n",
            "rvASOddCharacterNewton[Septic, 1] = [1/3: 1, 2/3: 1]\n",
            "rvASOddNewton[Septic] = [1/3: 6, 2/3: 6]\n",
            "rvASOddPoints[Cubic, 1] = 11\n",
            "rvASOddPoints[Twist, 1] = 6\n",
            "rvASOddPoints[Pure, 6] = 15126\n",
            "rvASOddPoints[Cubic, 6] = 14751\n",
            "rvASOddPoints[Quartic, 6] = 15726\n",
            "rvASOddPoints[Cubic, 14] = 6103437501\n",
            "rvASOddPoints[Twist, 14] = 6103359376\n",
        )
    );
}

#[test]
fn cyclotomic_coefficients_galois_norm_and_functional_equation() {
    let program = r"
        rvasodd C[5] {(6,1),(3,1)};
        rvASOddCharacterL[C,1];
        rvASOddCharacterL[C,4];
        show rvASOddCoeff[C,1,1];
        show rvASOddCoeff[C,4,1];
        assert rvASOddCoeff[C,1,0] == rvCyclo[5,0];
        assert rvASOddCoeff[C,2,6] == zero[RVCyclo(5)];
        rvASOddLPolynomial[C];
    ";
    let output = stdout(interpret(program));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines[0], "rvASOddCharacterL[C, 1] = [[1, 0, 0, 0], [2, 1, 2, 0], [-4, -3, -2, -6], [-10, 0, -5, -10], [0, 0, -25, 0], [50, 50, 0, 25]]");
    assert_eq!(lines[1], "rvASOddCharacterL[C, 4] = [[1, 0, 0, 0], [1, -1, -1, 1], [-1, 3, -3, 1], [-10, 0, -10, -5], [0, 0, 0, -25], [0, -50, -25, -50]]");
    assert_eq!(lines[2], "RVCyclo(5){2, 1, 2}");
    assert_eq!(lines[3], "RVCyclo(5){1, -1, -1, 1}");
    let coefficients: Vec<i128> = lines[4]
        .split_once(" = ")
        .unwrap()
        .1
        .trim_matches(['[', ']'])
        .split(", ")
        .map(|entry| entry.parse().unwrap())
        .collect();
    assert_eq!(
        coefficients,
        [
            1, 5, 10, 0, 0, 125, 250, -625, -3000, -2500, 1875, -12500, -75000, -78125, 156250,
            390625, 0, 0, 3906250, 9765625, 9765625
        ]
    );
    for i in 0..=10 {
        assert_eq!(
            coefficients[20 - i],
            5i128.pow(10 - i as u32) * coefficients[i]
        );
    }
}

#[test]
fn linear_covers_and_signed_extreme_coefficients_are_exact() {
    let program = r"
        rvasodd Line[7] {(7,1),(0,3)};
        rvasodd Min[3] {(2,-170141183460469231731687303715884105728)};
        rvasodd Reference[3] {(2,1)};
        rvASOddReduced[Line];
        rvASOddGenus[Line];
        rvASOddCharacterL[Line,6];
        rvASOddNewton[Line];
        rvASOddPoints[Line,8];
        rvASOddReduced[Min];
        assert rvASOddCoeff[Min,2,1] == rvASOddCoeff[Reference,2,1];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvASOddReduced[Line] = [(0, 3), (1, 1)]\n",
            "rvASOddGenus[Line] = 0\n",
            "rvASOddCharacterL[Line, 6] = [[1, 0, 0, 0, 0, 0]]\n",
            "rvASOddNewton[Line] = []\n",
            "rvASOddPoints[Line, 8] = 5764802\n",
            "rvASOddReduced[Min] = [(2, 1)]\n",
        )
    );
}

#[test]
fn validation_separates_normalization_character_and_size_errors() {
    let setup = r"
        rvasodd Small[5] {(2,1)};
        rvasodd Large[3] {(1000,1)};
        rvASOddGenus[Large];
        assert rvASOddCoeff[Large,2,0] == rvCyclo[3,0];
    ";
    assert_eq!(stdout(interpret(setup)), "rvASOddGenus[Large] = 999\n");
    for (statement, reason) in [
        ("rvasodd Bad[2] {(1,1)};", "requires an odd prime"),
        ("rvasodd Bad[9] {(1,1)};", "9 is not prime"),
        (
            "rvasodd Bad[5] {(25,1),(1,-1),(0,1)};",
            "nonconstant reduced polynomial",
        ),
        ("rvasodd Small[5] {(1,1)};", "already defined"),
        ("rvASOddGenus[Missing];", "unknown odd Artin–Schreier cover"),
        (
            "show rvASOddCoeff[Small,0,0];",
            "character index must lie in 1..p-1",
        ),
        (
            "rvASOddCharacterNewton[Small,5];",
            "character index must lie in 1..p-1",
        ),
        (
            "rvASOddPoints[Small,0];",
            "extension degree must be positive",
        ),
        ("rvASOddPoints[Small,128];", "point count is out of range"),
        ("rvASOddCharacterL[Large,1];", "computation limit"),
        (
            "rvASOddLPolynomial[Large];",
            "cyclotomic coefficient is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn point_traces_use_the_character_even_when_the_full_norm_is_too_wide() {
    let setup = "rvasodd C[101] {(2,1)};";
    assert_eq!(
        stdout(interpret(&format!("{setup} rvASOddPoints[C,2];"))),
        "rvASOddPoints[C, 2] = 102\n",
    );
    let error = diagnostic(interpret(&format!("{setup} rvASOddLPolynomial[C];")));
    assert!(
        error.contains("cyclotomic coefficient is out of range"),
        "{error}"
    );
}
