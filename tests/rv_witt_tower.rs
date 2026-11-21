mod common;

use common::{diagnostic, example, interpret, stdout};

fn coefficients(values: &str) -> Vec<i128> {
    values
        .trim_matches(['[', ']'])
        .split(", ")
        .map(|value| value.parse().unwrap())
        .collect()
}

#[test]
fn degree_sixteen_numerator_uses_only_the_primitive_conductors() {
    let output = stdout(example("rank-witt-tower"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        &lines[..7],
        [
            "rvWittTowerBreaks[CubicLength] = [1, 2, 4]",
            "rvWittTowerGenus[CubicLength] = 7",
            "rvWittTowerCharacterL[CubicLength, 3] = [[1, 0, 0, 0], [1, 1, 0, 0], [1, 1, 1, -1], [2, 0, 2, 0]]",
            "rvWittTowerLPolynomial[CubicLength] = [1, 6, 20, 48, 90, 140, 196, 272, 392, 560, 720, 768, 640, 384, 128]",
            "rvWittTowerNewton[CubicLength] = [1/4: 4, 1/2: 6, 3/4: 4]",
            "rvWittTowerPoints[CubicLength, 8] = 177",
            "rvWittTowerBreaks[QuarticLength] = [1, 2, 4, 8]",
        ],
    );
    assert_eq!(lines[7], "rvWittTowerGenus[QuarticLength] = 35");
    let polynomial = coefficients(lines[8].split_once(" = ").unwrap().1);
    assert_eq!(polynomial.len(), 71);
    assert_eq!(
        &polynomial[..10],
        &[1, 14, 104, 544, 2244, 7752, 23288, 62496, 153034, 348092]
    );
    assert_eq!(polynomial[35], 141_811_097_600);
    assert!((0..=35).all(|i| polynomial[70 - i] == polynomial[i] * (1i128 << (35 - i))));
    assert_eq!(
        &lines[9..],
        [
            "rvWittTowerNewton[QuarticLength] = [1/8: 8, 1/4: 12, 3/8: 8, 1/2: 14, 5/8: 8, 3/4: 12, 7/8: 8]",
            "rvWittTowerPoints[QuarticLength, 1] = 17",
            "rvWittTowerPoints[QuarticLength, 2] = 17",
            "rvWittTowerPoints[QuarticLength, 3] = 17",
            "rvWittTowerPoints[QuarticLength, 4] = 17",
            "rvWittTowerPoints[QuarticLength, 5] = 17",
            "rvWittTowerPoints[QuarticLength, 8] = 273",
        ],
    );
}

#[test]
fn the_first_two_quotients_match_the_quadratic_and_quartic_implementations() {
    let program = r"
        let x3 = rvmon[[0]]{3};
        let zero_coordinate = zero[RVPoly([0],0)];
        rvascover E[[0,1]] {rvmon[[0,1]]{0,3},rvmon[[0,1]]{3,0}};
        rvwittcover W2[[0]] {x3,zero_coordinate};
        rvwitttower T1[[0]] {x3};
        rvwitttower T2[[0]] {x3,zero_coordinate};
        rvASLPolynomial[E];
        rvWittTowerCharacterL[T2,1];
        rvWittTowerLPolynomial[T1];
        rvWittLPolynomial[W2];
        rvWittTowerLPolynomial[T2];
        rvWittCharacterL[W2];
        rvWittTowerCharacterL[T2,2];
        rvWittPoints[W2,8];
        rvWittTowerPoints[T2,8];
    ";
    let out = stdout(interpret(program));
    let values: Vec<_> = out
        .lines()
        .map(|row| row.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(&values[..3], &["[1, 0, 2]", "[[1], [0], [2]]", "[1, 0, 2]"]);
    assert_eq!(values[3], values[4]);
    assert_eq!(values[5].replace('(', "[").replace(')', "]"), values[6]);
    assert_eq!(&values[7..], &["209", "209"]);
}

#[test]
fn mixed_coordinates_change_upper_jumps_and_constants_twist_characters() {
    let program = r"
        let x = rvmon[[0]]{1};
        let z = zero[RVPoly([0],0)];
        rvwitttower Mixed[[0]] {x,rvmon[[0]]{3},z};
        rvwitttower Twist[[0]] {x,z,rvmon[[]]{}};
        rvWittTowerBreaks[Mixed];
        rvWittTowerGenus[Mixed];
        rvWittTowerCharacterL[Mixed,3];
        rvWittTowerNewton[Mixed];
        rvWittTowerPoints[Mixed,8];
        rvWittTowerLPolynomial[Twist];
        rvWittTowerCharacterL[Twist,3];
        rvWittTowerPoints[Twist,7];
        rvWittTowerPoints[Twist,8];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvWittTowerBreaks[Mixed] = [1, 3, 6]\n",
            "rvWittTowerGenus[Mixed] = 12\n",
            "rvWittTowerCharacterL[Mixed, 3] = [[1, 0, 0, 0], [1, 0, 0, 1], [1, 0, -1, 2], [0, 2, -2, 2], [-2, 2, -2, 2], [-4, 0, -4, 0]]\n",
            "rvWittTowerNewton[Mixed] = [1/4: 8, 1/2: 8, 3/4: 8]\n",
            "rvWittTowerPoints[Mixed, 8] = 209\n",
            "rvWittTowerLPolynomial[Twist] = [1, -2, 4, -8, 10, -20, 36, -48, 72, -80, 80, -128, 128, -128, 128]\n",
            "rvWittTowerCharacterL[Twist, 3] = [[1, 0, 0, 0], [-1, -1, 0, 0], [1, 1, 1, -1], [-2, 0, -2, 0]]\n",
            "rvWittTowerPoints[Twist, 7] = 169\n",
            "rvWittTowerPoints[Twist, 8] = 177\n",
        ),
    );
}

#[test]
fn validation_keeps_large_genus_and_point_queries_available() {
    let setup = r"
        let x = rvmon[[0]]{1};
        let z = zero[RVPoly([0],0)];
        rvwitttower W[[0]] {x,z,z,z,z,z,z,z};
        rvwitttower Base[[0]] {x};
        rvWittTowerBreaks[W];
        rvWittTowerGenus[W];
        rvWittTowerPoints[W,1];
        rvWittTowerCharacterL[W,1];
        rvWittTowerLPolynomial[Base];
        rvWittTowerNewton[Base];
    ";
    assert_eq!(
        stdout(interpret(setup)),
        concat!(
            "rvWittTowerBreaks[W] = [1, 2, 4, 8, 16, 32, 64, 128]\n",
            "rvWittTowerGenus[W] = 10795\n",
            "rvWittTowerPoints[W, 1] = 257\n",
            "rvWittTowerCharacterL[W, 1] = [[1]]\n",
            "rvWittTowerLPolynomial[Base] = [1]\n",
            "rvWittTowerNewton[Base] = []\n",
        )
    );
    for (statement, reason) in [
        ("rvwitttower Q[[0]] {};", "between 1 and 8"),
        ("rvwitttower Q[[0]] {x,z,z,z,z,z,z,z,z};", "between 1 and 8"),
        ("rvwitttower Q[[0,1]] {x};", "exactly one parameter"),
        ("rvwitttower Q[[0]] {z,x};", "coordinate 0 must have odd"),
        (
            "rvwitttower Q[[0]] {x,rvMul(x,x)};",
            "coordinate 1 must have odd",
        ),
        ("rvwitttower W[[0]] {x};", "already defined"),
        (
            "rvWittTowerPoints[W,0];",
            "extension degree is out of range",
        ),
        ("rvWittTowerPoints[W,19];", "exceeds the storage limit"),
        (
            "rvWittTowerCharacterL[W,9];",
            "quotient height is out of range",
        ),
        ("rvWittTowerCharacterL[W,8];", "character degree exceeds"),
        ("rvWittTowerLPolynomial[W];", "coefficient is out of range"),
        ("rvWittTowerGenus[Missing];", "unknown Witt tower"),
    ] {
        let result = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(result.contains(reason), "{statement}\n{result}");
    }
}
