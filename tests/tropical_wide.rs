mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn both_places_and_integer_rank_three_preserve_exact_support() {
    let output = stdout(example("rv-tropical-wide"));
    let valuations: Vec<_> = output
        .lines()
        .map(|line| line.strip_prefix("rvValPlucker = ").unwrap())
        .collect();
    assert_eq!(
        valuations,
        ["46", "46", "infinity", "0", "0", "0", "infinity", "0", "0", "0", "1", "0"],
    );
}

#[test]
fn wide_quadratic_minors_allow_unrepresentable_products_and_genuine_row_denominators() {
    let setup = r"
        rvqpadic K[7,2,3];
        let x = rvQpQuadPow[46](rvQpQuad[K,3,-1,1]);
        let z = zero[RVQpQuad(K)];
        let o = rvQpQuad[K,1,0,1];
        let two = rvQpQuad[K,2,0,1];
    ";
    let result = stdout(interpret(&format!(
        "{setup}
        let Q = rvValMatQuad[K,3]{{[x,x,z],[x,neg(x),o]}};
        rvValPlucker[[0,1]](Q);
        rvValPlucker[[0,2]](Q);
        rvValPlucker[[1,2]](Q);
        let Base = rvValMatQuad[K,3]{{[x,o,z],[x,two,o]}};
        let cx = rvQpQuadMul(x,rvQpQuad[K,1,0,11]);
        let dx = rvQpQuadMul(x,rvQpQuad[K,1,0,147]);
        let Scaled = rvValMatQuad[K,3]{{
            [cx,rvQpQuad[K,1,0,11],z],
            [dx,rvQpQuad[K,2,0,147],rvQpQuad[K,1,0,147]]
        }};
        assert Base == Scaled;
        assert Base == rvValMatQuad[K,3]{{[x,two,o],[x,o,z]}};
        let Rational = rvValMatQuad[K,3]{{[x,z,o],[o,o,z]}};
        assert Rational == rvValMatQuad[K,3]{{[o,o,z],[x,z,o]}};
        rvValPlucker[[0,1]](Rational);
        rvValPlucker[[0,2]](Rational);
        rvValPlucker[[1,2]](Rational);"
    )));
    assert_eq!(
        result,
        concat!(
            "rvValPlucker = 46
rvValPlucker = 0
rvValPlucker = 0
",
            "rvValPlucker = 46
rvValPlucker = 0
rvValPlucker = 0
",
        )
    );
    let error = diagnostic(interpret(&format!("{setup}show rvQpQuadInv(x);")));
    assert!(
        error.contains("exact integer arithmetic is out of range"),
        "{error}"
    );
}

#[test]
fn integer_projective_shear_cancels_before_valuation_and_respects_column_permutations() {
    let output = stdout(interpret(
        r"
        let N = rvValMat[5,3]{
           [10000000000000000000000000,10000000000000000000000001,0],
           [9999999999999999999999999,10000000000000000000000000,1]
        };
        let P = rvValMat[5,3]{
           [1,9999999999999999999999999,10000000000000000000000000],
           [0,10000000000000000000000000,10000000000000000000000001]
        };
        rvValPlucker[[0,1]](N); rvValPlucker[[0,2]](N); rvValPlucker[[1,2]](N);
        rvValPlucker[[1,2]](P); rvValPlucker[[0,1]](P); rvValPlucker[[0,2]](P);
    ",
    ));
    assert_eq!(
        output,
        "rvValPlucker = 0
rvValPlucker = 25
rvValPlucker = 0
"
        .repeat(2)
    );
}

#[test]
fn a_large_rank_with_a_noninvertible_pivot_keeps_the_normalized_full_basis() {
    let r = 64;
    let mut inputs = Vec::new();
    for i in 0..r {
        let row: Vec<_> = (0..r)
            .map(|j| match (i, j) {
                (0 | 1, 0) => "x",
                (1, 1) => "two",
                (0, 1) => "one",
                (i, j) if i == j => "one",
                _ => "z",
            })
            .collect();
        inputs.push(format!("[{}]", row.join(",")));
    }
    let result = stdout(interpret(&format!(
        "rvqpadic K[7,2,3]; let x=rvQpQuadPow[46](rvQpQuad[K,3,-1,1]);
        let z=zero[RVQpQuad(K)]; let one=rvQpQuad[K,1,0,1]; let two=rvQpQuad[K,2,0,1];
        rvValSupport(rvValMatQuad[K,{r}]{{{}}});",
        inputs.join(",")
    )));
    assert_eq!(
        result,
        format!(
            "rvValSupport = [{:?}]
",
            (0..r).collect::<Vec<_>>()
        )
    );
}
