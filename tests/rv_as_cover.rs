mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn exact_l_polynomials_and_newton_slopes_distinguish_the_covers() {
    assert_eq!(
        stdout(example("rank-artin-schreier")),
        concat!(
            "rvASGenus[Supersingular] = 1\n",
            "rvASPrank[Supersingular] = 0\n",
            "rvASLPolynomial[Supersingular] = [1, 0, 2]\n",
            "rvASNewton[Supersingular] = [1/2: 2]\n",
            "rvASPoints[Supersingular, 4] = 9\n",
            "rvASGenus[Ordinary] = 1\n",
            "rvASPrank[Ordinary] = 1\n",
            "rvASLPolynomial[Ordinary] = [1, 1, 2]\n",
            "rvASNewton[Ordinary] = [0: 1, 1: 1]\n",
            "rvASPoints[Ordinary, 5] = 44\n",
            "rvASGenus[GenusTwo] = 2\n",
            "rvASLPolynomial[GenusTwo] = [1, 0, 0, 0, 4]\n",
            "rvASNewton[GenusTwo] = [1/2: 4]\n",
            "rvASPoints[GenusTwo, 4] = 33\n",
            "rvASGenus[QuadraticPoles] = 1\n",
            "rvASPrank[QuadraticPoles] = 1\n",
            "rvASLPolynomial[QuadraticPoles] = [1, -1, 2]\n",
            "rvASNewton[QuadraticPoles] = [0: 1, 1: 1]\n",
            "rvASPoints[QuadraticPoles, 2] = 8\n",
        )
    );
}

#[test]
fn coordinate_and_artin_schreier_changes_preserve_the_counts() {
    let program = r"
        let s3 = rvmon[[0,1]]{3,0};
        let t3 = rvmon[[0,1]]{0,3};
        rvascover AtInfinity[[0,1]] {t3,s3};
        rvascover AtZero[[0,1]] {s3,t3};
        rvascover Shifted[[0,1]] {
            add(t3,add(rvmon[[0,1]]{1,2},rvmon[[0,1]]{2,1})),s3
        };
        rvASLPolynomial[AtInfinity];
        rvASLPolynomial[AtZero];
        rvASLPolynomial[Shifted];
        rvASPoints[AtInfinity,6];
        rvASPoints[AtZero,6];
        rvASPoints[Shifted,6];
        rvASPoints[Shifted,6];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvASLPolynomial[AtInfinity] = [1, 0, 2]\n",
            "rvASLPolynomial[AtZero] = [1, 0, 2]\n",
            "rvASLPolynomial[Shifted] = [1, 0, 2]\n",
            "rvASPoints[AtInfinity, 6] = 81\n",
            "rvASPoints[AtZero, 6] = 81\n",
            "rvASPoints[Shifted, 6] = 81\n",
            "rvASPoints[Shifted, 6] = 81\n",
        )
    );
}

#[test]
fn genus_zero_and_a_cover_with_repeated_and_nonrational_poles() {
    let program = r"
        let s = rvmon[[0,1]]{1,0};
        let t = rvmon[[0,1]]{0,1};
        rvascover Line[[0,1]] { t,add(s,t) };
        rvASGenus[Line];
        rvASPrank[Line];
        rvASLPolynomial[Line];
        rvASNewton[Line];
        rvASPoints[Line,5];

        rvascover Mixed[[0,1]] {
            rvmon[[0,1]]{5,0},
            add(rvmon[[0,1]]{0,5},add(rvmon[[0,1]]{1,4},rvmon[[0,1]]{2,3}))
        };
        rvASGenus[Mixed];
        rvASPrank[Mixed];
        rvASPoints[Mixed,1];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvASGenus[Line] = 0\n",
            "rvASPrank[Line] = 0\n",
            "rvASLPolynomial[Line] = [1]\n",
            "rvASNewton[Line] = []\n",
            "rvASPoints[Line, 5] = 33\n",
            "rvASGenus[Mixed] = 3\n",
            "rvASPrank[Mixed] = 2\n",
            "rvASPoints[Mixed, 1] = 3\n",
        )
    );
}

#[test]
fn malformed_even_and_geometrically_non_coprime_poles_are_rejected() {
    let setup = r"
        let s = rvmon[[0,1]]{1,0};
        let t = rvmon[[0,1]]{0,1};
        let s2 = rvMul(s,s);
        let t2 = rvMul(t,t);
        let q = add(add(s2,t2),rvMul(s,t));
        rvascover C[[0,1]] {t,s};
    ";
    for (statement, reason) in [
        ("rvascover K[[0,1]] {};", "numerator and a denominator"),
        (
            "rvascover K[[0]] {rvmon[[0]]{1},rvmon[[0]]{1}};",
            "exactly two",
        ),
        ("rvascover K[[0,1]] {s,t2};", "same positive degree"),
        (
            "rvascover K[[0,1]] {s,zero[RVPoly([0,1],1)]};",
            "denominator must be nonzero",
        ),
        (
            "rvascover K[[0,1]] {s2,rvMul(s,t)};",
            "common projective zero at [0:1]",
        ),
        (
            "rvascover K[[0,1]] {rvMul(q,s),rvMul(q,t)};",
            "common projective zeros on the chart",
        ),
        (
            "rvascover K[[0,1]] {t2,s2};",
            "pole at [0:1] must have odd order",
        ),
        (
            "rvascover K[[0,1]] {s2,t2};",
            "finite poles must all have odd order",
        ),
        (
            "rvascover K[[0,1]] {rvmon[[0,1]]{5,0},rvMul(rvMul(q,q),add(s,t))};",
            "finite poles must all have odd order",
        ),
        ("rvascover C[[0,1]] {t,s};", "already defined"),
        ("rvASGenus[Missing];", "unknown Artin–Schreier cover"),
        (
            "rvASPoints[C,0];",
            "point-count extension degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn sparse_high_poles_and_point_work_limits_are_checked() {
    let h = 1usize << (usize::BITS / 2);
    let prefix = format!(
        "rvascover Large[[0,1]] {{rvmon[[0,1]]{{{},0}},rvmon[[0,1]]{{0,{}}}}};",
        h + 1,
        h + 1,
    );
    assert_eq!(
        stdout(interpret(&format!(
            "{prefix}\nrvASGenus[Large]; rvASPrank[Large]; rvASPoints[Large,1];"
        ))),
        format!(
            "rvASGenus[Large] = {}\nrvASPrank[Large] = 0\nrvASPoints[Large, 1] = 3\n",
            h / 2
        )
    );
    let error = diagnostic(interpret(&format!("{prefix}\nrvASLPolynomial[Large];")));
    assert!(
        error.contains("exceeds the point-count extension range"),
        "{error}"
    );
    let error = diagnostic(interpret(&format!(
        "{prefix}\nrvASPoints[Large,{}];",
        usize::BITS - 1
    )));
    assert!(error.contains("exceeds the computation limit"), "{error}");
}
