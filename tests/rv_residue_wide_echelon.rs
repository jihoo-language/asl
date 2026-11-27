mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn two_split_primes_retain_the_signed_angular_ratio_and_drop_dependent_rows() {
    assert_eq!(
        stdout(example("rv-residue-wide-echelon")),
        concat!(
            "RVResSpace(7, 3){[1, 4, 0], [0, 0, 1]}\n",
            "rvResSupport = [[0, 2], [1, 2]]\n",
            "RVResSpace(7, 3){[1, 4, 0], [0, 0, 1]}\n",
        )
    );
}

#[test]
fn weighted_pivoting_uses_differences_even_at_both_signed_endpoints() {
    let program = r"
        rvqpadic K[7,2,3];
        let z = zero[RVQpQuad(K)];
        let one = rvQpQuad[K,1,0,1];
        let x = rvQpQuadMul(rvQpQuad[K,2,0,1],rvQpQuadPow[46](rvQpQuad[K,3,-1,1]));
        assert rvResSpecializeQuad[K,4,[9223372036854775807,9223372036854775761,-9223372036854775808,0]]{
            [x,one,z,z],[z,z,one,one],[x,one,one,one]
        } == rvResSpace[7,4]{[1,4,0,0],[0,0,0,1]};
        assert rvResSpecializeQuad[K,4,[9223372036854775761,9223372036854775807,0,-9223372036854775808]]{
            [one,x,z,z],[z,z,one,one],[one,x,one,one]
        } == rvResSpace[7,4]{[1,2,0,0],[0,0,1,0]};
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn padded_square_chart_has_full_rank_without_enumerating_bases() {
    let weights = [46]
        .into_iter()
        .chain(std::iter::repeat_n(0, 63))
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let rows = (0..64)
        .map(|i| {
            let row = (0..64)
                .map(|j| match (i, j) {
                    (0, 0) => "x",
                    (0, 1) => "one",
                    (i, j) if i == j => "one",
                    _ => "z",
                })
                .collect::<Vec<_>>()
                .join(",");
            format!("[{row}]")
        })
        .collect::<Vec<_>>()
        .join(",");
    let program = format!(
        "rvqpadic K[7,2,3];
         let z = zero[RVQpQuad(K)];
         let one = rvQpQuad[K,1,0,1];
         let x = rvQpQuadPow[46](rvQpQuad[K,3,-1,1]);
         let S = rvResSpecializeQuad[K,64,[{weights}]]{{{rows}}};
         rvResDim(S);"
    );
    assert_eq!(stdout(interpret(&program)), "rvResDim = 64\n");
}

#[test]
fn matrix_shape_and_denominator_errors_keep_their_diagnostics() {
    for (program, expected) in [
        (
            "show rvResSpecialize[7,2,[0]]{[1,0],[0,1]};",
            "needs 2 weights",
        ),
        (
            "show rvResSpecialize[7,2,[0,0]]{[1,0],[(0,0),1]};",
            "denominator must be positive",
        ),
        (
            "show rvResSpecialize[7,2,[0,0]]{[1,0],[1]};",
            "rows of length 2",
        ),
    ] {
        let output = diagnostic(interpret(program));
        assert!(output.contains(expected), "{program}\n{output}");
    }
}
