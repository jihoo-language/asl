mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn example_matches_the_quadratic_field_middle_gram_matrix() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-hodge-wick")),
        concat!(
            "rvHodgeTrace[E] = 1\n",
            "rvHodgeTrace[E] = 4/3\n",
            "rvHodgePair = 12\n",
            "rvHodgePair = 44\n",
            "rvHodgePair = 22\n",
            "rvHodgePair = 528\n",
            "rvHodgePair = 264\n",
            "rvHodgePair = 187\n",
            "rvHodgePair = 5720\n",
            "rvHodgeIntegral = 15525\n",
        )
    );
}

fn combine((a, b): (i128, i128), (c, d): (i128, i128)) -> (i128, i128) {
    (a * c + 5 * b * d, a * d + b * c)
}

fn moment(rank: i128, power: usize) -> i128 {
    (0..power).map(|j| rank + 2 * j as i128).product()
}

fn independent_quadratic_moment(inputs: &[(i128, i128)]) -> i128 {
    // Evaluate directly after splitting E, using independent chi-square
    // variables on the two eigenspaces. Coefficients are pairs for a+b√5.
    let mut mixed = vec![(1, 0)];
    for &(a, b) in inputs {
        let mut next = vec![(0, 0); mixed.len() + 1];
        for (j, &value) in mixed.iter().enumerate() {
            let x = combine(value, (a, b));
            let y = combine(value, (a, -b));
            next[j + 1].0 += x.0;
            next[j + 1].1 += x.1;
            next[j].0 += y.0;
            next[j].1 += y.1;
        }
        mixed = next;
    }
    let (mut result, mut irrational) = (0, 0);
    for (j, (a, b)) in mixed.into_iter().enumerate() {
        let moments = moment(11, j) * moment(11, inputs.len() - j);
        result += a * moments;
        irrational += b * moments;
    }
    assert_eq!(irrational, 0);
    result
}

#[test]
fn mixed_rational_traces_and_cumulants_match_a_direct_eigenvalue_calculation() {
    let input = [(1, 3), (-2, 1), (0, -2), (5, -2), (1, 0)];
    let mut program = "rvhodgewick E[11]{-5,0,1}{[(3,2)]};\n".to_owned();
    let mut expected = Vec::new();
    for n in 2..=input.len() {
        let factors: Vec<_> = input[..n]
            .iter()
            .map(|(a, b)| format!("rvHodgeQ[E,{n}]{{{a},{b}}}"))
            .collect();
        let product = factors
            .into_iter()
            .reduce(|left, right| format!("rvHodgeCup({left},{right})"))
            .unwrap();
        writeln!(program, "rvHodgeIntegral({product});").unwrap();
        expected.push(format!(
            "rvHodgeIntegral = {}",
            independent_quadratic_moment(&input[..n])
        ));
    }
    program.push_str(
        r"
        rvHodgeTrace[E]{(1,2),(1,3)};
        rvHodgePair(rvHodgeQ[E,2]{0,1},rvHodgeQ[E,2]{0,1});
        rvHodgePair(rvHodgeQ[E,2]{(1,2),(1,3)},rvHodgeQ[E,2]{(1,2),(1,3)});
        rvHodgePair(rvHodgeCup(rvHodgeD[E,2,0],rvHodgeD[E,2,0]),rvHodgeQ[E,2]{1});
        ",
    );
    expected.extend(
        [
            "rvHodgeTrace[E] = 1",
            "rvHodgePair = 220",
            "rvHodgePair = 1408/9",
            "rvHodgePair = 33",
        ]
        .map(str::to_owned),
    );
    assert_eq!(
        stdout(interpret(&program)).lines().collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn inverse_beauville_powers_and_surface_normalization_include_nondiagonal_gram() {
    let mut program = "rvhodgewick A[21]{0,1}{[0,(3,2)],[(3,2),(1,3)]};\n".to_owned();
    let mut expected = Vec::new();
    for n in 2..=6 {
        let q = format!("rvHodgeInverseBB[A,{n}]");
        let product = std::iter::repeat_n(q, n)
            .reduce(|a, b| format!("rvHodgeCup({a},{b})"))
            .unwrap();
        writeln!(program, "rvHodgeIntegral({product});").unwrap();
        expected.push(format!("rvHodgeIntegral = {}", moment(23, n)));
    }
    program.push_str(
        r"
        rvhodgewick S[21]{0,1}{[2]};
        rvHodgePair(rvHodgeD[S,1,0],rvHodgeD[S,1,0]);
        rvHodgeIntegral(rvHodgeQ[S,1]{1});
        rvHodgeIntegral(rvHodgeInverseBB[S,1]);
        ",
    );
    expected.extend(
        [
            "rvHodgePair = 2",
            "rvHodgeIntegral = 21",
            "rvHodgeIntegral = 22",
        ]
        .map(str::to_owned),
    );
    assert_eq!(
        stdout(interpret(&program)).lines().collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn real_etale_idempotents_and_typed_rational_operations_are_exact() {
    let text = stdout(interpret(
        r"
        rvhodgewick A[11]{-1,0,1}{[2]};
        let p=rvHodgeQ[A,2]{(1,2),(1,2)};
        let q=rvHodgeQ[A,2]{(1,2),(-1,2)};
        let identity=rvHodgeUnit[A,2];
        assert add(p,q)==rvHodgeQ[A,2]{1};
        assert rvHodgeCup(identity,p)==p;
        assert add(p,neg(p))==zero[RVHodgeTensor(A,2,4)];
        assert rvHodgeScale[2,3](p)==rvHodgeQ[A,2]{(1,3),(1,3)};
        assert rvHodgeCup(p,q)==rvHodgeCup(q,p);
        assert rvHodgeCup(rvHodgeCup(p,q),p)==zero[RVHodgeTensor(A,2,12)];
        rvHodgeTrace[A]{(1,2),(1,2)};
        rvHodgePair(p,p);
        rvHodgePair(q,q);
        rvHodgePair(q,p);
        rvHodgePair(rvHodgeCup(p,q),identity);
        ",
    ));
    assert_eq!(
        text,
        "rvHodgeTrace[A] = 1\nrvHodgePair = 143\nrvHodgePair = 143\nrvHodgePair = 121\nrvHodgePair = 121\n"
    );
}

#[test]
fn validates_etale_algebra_rational_gram_and_nominal_grades() {
    let inputs = [
        ("rvhodgewick B[11]{1,0,1}{[2]};", "all roots real"),
        ("rvhodgewick B[11]{0,0,1}{[2]};", "square-free"),
        ("rvhodgewick B[11]{-1,0,2}{[2]};", "monic"),
        ("rvhodgewick B[9]{-1,0,1}{[2]};", "total dimension 22 or 23"),
        ("rvhodgewick B[21]{0,1}{[1,0],[1,1]};", "symmetric"),
        ("rvhodgewick B[21]{0,1}{[1,1],[1,1]};", "singular"),
        ("rvhodgewick B[22]{0,1}{[(1,0)]};", "denominator must be positive"),
        ("rvhodgewick B[11]{-1,0,1}{[2]}; rvhodgewick B[11]{-1,0,1}{[2]};", "already defined"),
        ("rvhodgewick B[11]{-1,0,1}{[2]}; rvHodgeIntegral(rvHodgeQ[B,1]{1});", "requires total Beauville rank 22"),
        ("rvhodgewick B[11]{-1,0,1}{[2]}; rvHodgeTrace[B]{1,2,3};", "too many power-basis coordinates"),
        ("rvhodgewick B[11]{-1,0,1}{[2]}; rvHodgePair(rvHodgeD[B,2,1],rvHodgeD[B,2,0]);", "outside the marked Gram"),
        ("rvhodgewick B[11]{-1,0,1}{[2]}; rvHodgePair(rvHodgeQ[B,2]{1},rvHodgeD[B,2,0]);", "requires complementary degrees"),
        ("rvhodgewick B[11]{-1,0,1}{[2]}; rvhodgewick C[11]{-1,0,1}{[2]}; rvHodgePair(rvHodgeQ[B,2]{1},rvHodgeQ[C,2]{1});", "same marking"),
    ];
    for (program, reason) in inputs {
        assert!(diagnostic(interpret(program)).contains(reason), "{program}");
    }
}

#[test]
fn a_compact_radial_moment_can_exceed_the_machine_integer_range() {
    let mut program = "rvhodgewick A[23]{0,1}{};\n".to_owned();
    let mut product = "rvHodgeInverseBB[A,30]".to_owned();
    for _ in 1..30 {
        product = format!("rvHodgeCup({product},rvHodgeInverseBB[A,30])");
    }
    writeln!(program, "rvHodgeIntegral({product});").unwrap();
    assert_eq!(
        stdout(interpret(&program)),
        "rvHodgeIntegral = 469988167894565540026942788743061790463241080859375\n"
    );
}
