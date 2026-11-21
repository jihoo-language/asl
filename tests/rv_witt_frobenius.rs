mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn primitive_character_slopes_and_frobenius_points() {
    assert_eq!(
        stdout(example("rank-witt-frobenius")),
        concat!(
            "rvWittTowerCharacterNewton[T, 4] = [1/8: 1, 1/4: 1, 3/8: 1, 1/2: 1, 5/8: 1, 3/4: 1, 7/8: 1]\n",
            "rvWittTowerNewton[T] = [1/8: 8, 1/4: 12, 3/8: 8, 1/2: 14, 5/8: 8, 3/4: 12, 7/8: 8]\n",
            "rvWittTowerCharacterBaseL[T, 3, 2] = [[1, 0, 0, 0], [1, 0, 1, -2], [-2, 0, -2, -4], [0, 0, -8, 0]]\n",
            "rvWittTowerCharacterBaseL[T, 3, 4] = [[1, 0, 0, 0], [-4, -4, -2, -4], [-16, 16, 8, 16], [64, 0, 0, 0]]\n",
            "rvWittTowerPoints[T, 8] = 273\n",
            "rvWittTowerFrobPoints[T, 8] = 273\n",
            "rvWittTowerFrobPoints[T, 32] = 4295365505\n",
            "rvWittTowerFrobPoints[T, 48] = 281474971197953\n",
            "rvWittTowerFrobPoints[T, 96] = 79228162514257559424621051905\n",
        )
    );
}

#[test]
fn odd_base_changes_keep_a_quadratic_twist_and_even_base_changes_remove_it() {
    let program = r"
        let x = rvmon[[0]]{1};
        let z = zero[RVPoly([0],0)];
        let unit = rvmon[[]]{};
        rvwitttower W[[0]] {x,z,z};
        rvwitttower Twist[[0]] {x,z,unit};
        rvWittTowerCharacterL[W,2];
        rvWittTowerCharacterBaseL[W,2,1];
        rvWittTowerCharacterBaseL[W,2,3];
        rvWittTowerCharacterBaseL[W,2,4];
        rvWittTowerCharacterBaseL[W,3,2];
        rvWittTowerCharacterBaseL[Twist,3,2];
        rvWittTowerCharacterBaseL[W,3,3];
        rvWittTowerCharacterBaseL[Twist,3,3];
        rvWittTowerCharacterNewton[W,1];
        rvWittTowerCharacterBaseL[W,1,10000];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvWittTowerCharacterL[W, 2] = [[1, 0], [1, 1]]\n",
            "rvWittTowerCharacterBaseL[W, 2, 1] = [[1, 0], [1, 1]]\n",
            "rvWittTowerCharacterBaseL[W, 2, 3] = [[1, 0], [-2, 2]]\n",
            "rvWittTowerCharacterBaseL[W, 2, 4] = [[1, 0], [4, 0]]\n",
            "rvWittTowerCharacterBaseL[W, 3, 2] = [[1, 0, 0, 0], [1, 0, 1, -2], [-2, 0, -2, -4], [0, 0, -8, 0]]\n",
            "rvWittTowerCharacterBaseL[Twist, 3, 2] = [[1, 0, 0, 0], [1, 0, 1, -2], [-2, 0, -2, -4], [0, 0, -8, 0]]\n",
            "rvWittTowerCharacterBaseL[W, 3, 3] = [[1, 0, 0, 0], [1, -3, 3, 1], [4, -4, 8, -8], [-16, 0, 16, 0]]\n",
            "rvWittTowerCharacterBaseL[Twist, 3, 3] = [[1, 0, 0, 0], [-1, 3, -3, -1], [4, -4, 8, -8], [16, 0, -16, 0]]\n",
            "rvWittTowerCharacterNewton[W, 1] = []\n",
            "rvWittTowerCharacterBaseL[W, 1, 10000] = [[1]]\n",
        )
    );
}

#[test]
fn the_trace_query_handles_wide_counts_with_explicit_bounds() {
    let setup = r"
        let x = rvmon[[0]]{1};
        let z = zero[RVPoly([0],0)];
        rvwitttower W[[0]] {x,z,z,z};
    ";
    assert_eq!(
        stdout(interpret(&format!(
            "{setup} rvWittTowerFrobPoints[W,64]; rvWittTowerFrobPoints[W,120];"
        ))),
        concat!(
            "rvWittTowerFrobPoints[W, 64] = 18446744062836078593\n",
            "rvWittTowerFrobPoints[W, 120] = 1329227995784915885637038743763288065\n",
        )
    );
    for (statement, reason) in [
        (
            "rvWittTowerCharacterBaseL[W,4,48];",
            "coefficient is out of range",
        ),
        (
            "rvWittTowerCharacterBaseL[W,4,0];",
            "extension degree must be positive",
        ),
        (
            "rvWittTowerFrobPoints[W,0];",
            "extension degree must be positive",
        ),
        (
            "rvWittTowerFrobPoints[W,128];",
            "extension degree exceeds 127",
        ),
        (
            "rvWittTowerCharacterBaseL[W,0,2];",
            "quotient height is out of range",
        ),
        (
            "rvWittTowerCharacterNewton[W,5];",
            "quotient height is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
