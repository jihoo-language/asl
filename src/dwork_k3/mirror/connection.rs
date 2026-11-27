use crate::error::Result;

use super::{integrity, Operator};

pub(super) type Polynomial = [i128; 3];
const ZERO: Polynomial = [0; 3];
pub(super) const D: Polynomial = [1, -256, 0];

pub(super) fn report() -> Result<String> {
    let [c0, c1, c2, c3] = Operator::PicardFuchs.coefficients();
    let [m0, m1, m2, _] = Operator::Clausen.coefficients();
    if z(c3) != subtract([1, 0, 0], D) || m2 != c3 {
        return Err(integrity());
    }
    // Monic rank-two equation theta^2 + a theta + b, with a=A2/D, b=B2/D.
    let (a, b) = (z(-m1), z(-m0));
    let da = quotient_derivative_numerator(a);
    let db = quotient_derivative_numerator(b);
    let square = [
        add(scale(db, 2), scale(product(a, b), 4)),
        add(add(da, scale(product(b, D), 4)), scale(product(a, a), 2)),
        scale(product(a, D), 3),
    ];
    if square != [c0, c1, c2].map(|c| product(z(-c), D)) {
        return Err(integrity());
    }

    // For (y,theta y,theta^2 y) the companion matrix is A/D. B comes
    // from the discriminant of the corresponding binary quadratic.
    let a_num = companion_three();
    let pairing = polarization();
    for i in 0..3 {
        for j in 0..3 {
            let mut derivative = product(D, theta(pairing[i][j]));
            for k in 0..3 {
                derivative = add(derivative, product(a_num[k][i], pairing[k][j]));
                derivative = add(derivative, product(pairing[i][k], a_num[k][j]));
            }
            if derivative != ZERO {
                return Err(integrity());
            }
        }
    }
    Ok(format!(
        "rvDworkMirrorConnection = D=1-256z, A=D^-1*{}, B={}, det B=D^3, certified symmetric square and horizontal pairing",
        display_matrix(a_num), display_matrix(pairing),
    ))
}

pub(super) fn companion_two() -> [[Polynomial; 2]; 2] {
    let [c0, c1, _, _] = Operator::Clausen.coefficients();
    [[ZERO, D], [z(c0), z(c1)]]
}

pub(super) fn companion_three() -> [[Polynomial; 3]; 3] {
    let [c0, c1, c2, _] = Operator::PicardFuchs.coefficients();
    [[ZERO, D, ZERO], [ZERO, ZERO, D], [z(c0), z(c1), z(c2)]]
}

pub(super) fn polarization() -> [[Polynomial; 3]; 3] {
    let [c0, c1, _, _] = Operator::Clausen.coefficients();
    [
        [z(-4 * c0), z(-c1), D],
        [z(-c1), scale(D, -1), ZERO],
        [D, ZERO, ZERO],
    ]
}

fn z(coefficient: i128) -> Polynomial {
    [0, coefficient, 0]
}

fn add(a: Polynomial, b: Polynomial) -> Polynomial {
    std::array::from_fn(|i| a[i] + b[i])
}

fn subtract(a: Polynomial, b: Polynomial) -> Polynomial {
    std::array::from_fn(|i| a[i] - b[i])
}

fn scale(a: Polynomial, constant: i128) -> Polynomial {
    a.map(|x| x * constant)
}

fn product(a: Polynomial, b: Polynomial) -> Polynomial {
    // All multiplicands in the polynomial checks have degree at most one.
    debug_assert!(a[2] == 0 && b[2] == 0);
    [a[0] * b[0], a[0] * b[1] + a[1] * b[0], a[1] * b[1]]
}

fn theta(a: Polynomial) -> Polynomial {
    [0, a[1], 2 * a[2]]
}

fn quotient_derivative_numerator(numerator: Polynomial) -> Polynomial {
    subtract(product(theta(numerator), D), product(numerator, theta(D)))
}

fn display_matrix(matrix: [[Polynomial; 3]; 3]) -> String {
    let rows = matrix.map(|row| {
        let row = row.map(|[c, a, quadratic]| {
            debug_assert_eq!(quadratic, 0);
            match (c, a) {
                (_, 0) => c.to_string(),
                (0, _) => format!("{a}z"),
                _ => format!("{c}{a:+}z"),
            }
        });
        format!("[{}]", row.join(", "))
    });
    format!("[{}]", rows.join(", "))
}
