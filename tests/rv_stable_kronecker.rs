mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const HUGE: u128 = u128::MAX;
const SYMS: &str = r"
    let a = rvStableSchurCharacter[[1]];
    let b = rvStableSchurCharacter[[2]];
    let c = rvStableSchurCharacter[[1,1]];
    let d = rvStableSchurCharacter[[2,1]];
    let u = rvStableSchurCharacter[[]];
";

#[test]
fn example_straightens_the_same_reduced_square_at_three_ranks() {
    assert_eq!(
        stdout(example("witt-stable-kronecker")),
        concat!(
            "RVStableSchurCharacter{([2],1), ([1, 1],1), ([1],1), ([],1)}\n",
            "rvStableSchurCharacterCoeff[[2]] = 1\n",
            "RVSchurCharacter(2){([2],1)}\n",
            "RVSchurCharacter(3){([3],1), ([2, 1],1), ([1, 1, 1],1)}\n",
            "RVSchurCharacter(4){([4],1), ([3, 1],1), ([2, 2],1), ([2, 1, 1],1)}\n",
            "rvStableSchurCharacterCoeff[[2, 1]] = 9\n",
        )
    );
}

#[test]
fn stable_powers_specialize_to_finite_internal_products_on_both_sides_of_the_bound() {
    let mut input = format!(
        r"{SYMS}
        let classes = rvStableSchurKronecker(a,a);
        assert classes == add(add(u,a),add(b,c));
        let dd = rvStableSchurKronecker(d,d);
        let bc = rvStableSchurKronecker(b,c);
        let ac = rvStableSchurKronecker(a,c);
    "
    );
    for (left, right, product, min, max) in [
        ("a", "a", "classes", 2, 8),
        ("d", "d", "dd", 5, 12),
        ("b", "c", "bc", 4, 8),
        ("a", "c", "ac", 3, 6),
    ] {
        for n in min..=max {
            writeln!(input,"assert rvStableSchurSpecialize[{n}]({product}) == rvSchurCharacterKronecker(rvStableSchurSpecialize[{n}]({left}),rvStableSchurSpecialize[{n}]({right}));").unwrap();
        }
    }
    let output = stdout(interpret(&input));
    assert!(output.is_empty());
}

#[test]
fn littlewood_triples_and_virtual_ring_identities_agree() {
    let input = format!(
        r"{SYMS}
        let bb = rvStableSchurKronecker(b,b);
        let dd = rvStableSchurKronecker(d,d);
        rvStableSchurCharacterCoeff[[2]](bb);
        rvStableSchurCharacterCoeff[[2,1]](dd);
        rvStableSchurCharacterCoeff[[3]](dd);
        rvStableSchurCharacterCoeff[[1,1,1]](dd);
        rvStableSchurCharacterCoeff[[7]](dd);
        let x = add(b,rvStableSchurCharacterScale[-2](c));
        let y = add(a,neg(d));
        assert rvStableSchurKronecker(x,y) == rvStableSchurKronecker(y,x);
        assert rvStableSchurKronecker(rvStableSchurKronecker(x,y),a) == rvStableSchurKronecker(x,rvStableSchurKronecker(y,a));
        assert rvStableSchurKronecker(x,y) == add(rvStableSchurKronecker(b,y),rvStableSchurCharacterScale[-2](rvStableSchurKronecker(c,y)));
        assert add(x,neg(x)) == zero[RVStableSchurCharacter];
        assert rvStableSchurKronecker(x,zero[RVStableSchurCharacter]) == zero[RVStableSchurCharacter];
        assert rvStableSchurKronecker(u,x) == x;
        assert rvStableSchurKronecker(neg(u),x) == neg(x);
    "
    );
    let output = stdout(interpret(&input));
    // In Littlewood's σ[XYZ+XY+XZ+YZ], tails of weight three on all
    // axes use either the degree-three ordinary Kronecker term or
    // degree-one terms on all four factors. The latter contributes
    // 2³ for (21,21,21) and 2² for (21,21,3) or (21,21,111).
    let values: Vec<_> = output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1)
        .collect();
    assert_eq!(values, ["2", "9", "5", "5", "0"]);
}

fn partitions(n: usize) -> Vec<Vec<usize>> {
    fn visit(n: usize, max: usize, prefix: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if n == 0 {
            out.push(prefix.clone());
            return;
        }
        for k in (1..=n.min(max)).rev() {
            prefix.push(k);
            visit(n - k, k, prefix, out);
            prefix.pop();
        }
    }
    let mut out = Vec::new();
    visit(n, n, &mut Vec::new(), &mut out);
    out
}

fn dimension(partition: &[usize]) -> usize {
    let numerator: usize = (1..=partition.iter().sum()).product();
    let hooks: usize = partition
        .iter()
        .enumerate()
        .flat_map(|(i, &width)| {
            (0..width)
                .map(move |j| width - j + partition[i + 1..].iter().filter(|&&row| row > j).count())
        })
        .product();
    numerator / hooks
}

#[test]
fn ordinary_products_satisfy_character_orthogonality_and_the_regular_character_identity() {
    let n = 6;
    let terms = partitions(n);
    let mut input = String::new();
    let mut regular = Vec::new();
    for (i, shape) in terms.iter().enumerate() {
        writeln!(input, "let s{i} = rvSchurCharacter[{shape:?}];").unwrap();
        regular.push(format!("rvSchurCharacterScale[{}](s{i})", dimension(shape)));
    }
    let regular = regular
        .into_iter()
        .reduce(|a, b| format!("add({a},{b})"))
        .unwrap();
    writeln!(input, "let regular = {regular};").unwrap();
    for (i, shape) in terms.iter().enumerate() {
        writeln!(
            input,
            "assert rvSchurCharacterKronecker(s{i},regular) == rvSchurCharacterScale[{}](regular);",
            dimension(shape)
        )
        .unwrap();
        for j in 0..terms.len() {
            writeln!(
                input,
                "rvSchurCharacterCoeff[[{n}]](rvSchurCharacterKronecker(s{i},s{j}));"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&input));
    let coefficients: Vec<usize> = output
        .lines()
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect();
    let expected: Vec<_> = (0..terms.len())
        .flat_map(|i| (0..terms.len()).map(move |j| usize::from(i == j)))
        .collect();
    assert_eq!(coefficients, expected);
}

#[test]
fn straightening_handles_collisions_negative_rows_and_wide_cancellations() {
    let almost = HUGE - 1;
    let input = format!(
        r"{SYMS}
        assert rvStableSchurSpecialize[0](u) == rvSchurCharacter[[]];
        assert rvStableSchurSpecialize[0](a) == neg(rvSchurCharacter[[]]);
        assert rvStableSchurSpecialize[0](c) == rvSchurCharacter[[]];
        assert rvStableSchurSpecialize[0](b) == zero[RVSchurCharacter(0)];
        assert rvStableSchurSpecialize[1](a) == zero[RVSchurCharacter(1)];
        assert rvStableSchurSpecialize[2](c) == zero[RVSchurCharacter(2)];
        assert rvStableSchurSpecialize[3](rvStableSchurCharacter[[2,1,1]]) == zero[RVSchurCharacter(3)];
        assert rvStableSchurSpecialize[1](b) == neg(rvSchurCharacter[[1]]);
        assert rvStableSchurSpecialize[0](rvStableSchurKronecker(a,a)) == rvSchurCharacter[[]];
        assert rvStableSchurSpecialize[{HUGE}](a) == rvSchurCharacter[[{almost},1]];
        let h = rvStableSchurCharacter[[{HUGE}]];
        assert rvStableSchurKronecker(h,u) == h;
        assert rvStableSchurSpecialize[{HUGE}](h) == neg(rvSchurCharacter[[{almost},1]]);
        assert rvSchurCharacterKronecker(rvSchurCharacter[[{HUGE}]],rvSchurCharacter[[{almost},1]]) == rvSchurCharacter[[{almost},1]];
        let large = add(rvStableSchurCharacterScale[-170141183460469231731687303715884105728](a),neg(u));
        assert rvStableSchurSpecialize[0](large) == rvSchurCharacterScale[170141183460469231731687303715884105727](rvSchurCharacter[[]]);
        let at3 = rvStableSchurSpecialize[3](rvStableSchurKronecker(a,a));
        let f = rvSchurPlethysm(rvSchurCharacter[[1]],at3);
        assert rvSchurCharacterEvaluate(f,rvWittSpectrumLine[2]) == rvWittSpectrumLine[8];
    "
    );
    assert!(stdout(interpret(&input)).is_empty());
}

#[test]
fn invalid_types_and_exact_limits_are_checked() {
    for (query, expected) in [
        ("show rvSchurCharacterKronecker(rvSchurCharacter[[2]],rvSchurCharacter[[1]]);", "RVSchurCharacter(2)"),
        ("show rvStableSchurSpecialize[2](rvSchurCharacter[[2]]);", "RVStableSchurCharacter"),
        ("show rvStableSchurKronecker(a,rvSchurCharacter[[1]]);", "RVStableSchurCharacter"),
        ("rvStableSchurCharacterCoeff[[1]](rvSchurCharacter[[1]]);", "RVStableSchurCharacter"),
        ("show rvStableSchurCharacter[[1,2]];", "weakly decreasing"),
        ("show rvStableSchurKronecker(rvStableSchurCharacter[[9]],rvStableSchurCharacter[[8]]);", "combined tail limit"),
        ("show rvSchurCharacterKronecker(rvSchurCharacter[[16,1]],rvSchurCharacter[[16,1]]);", "weight limit"),
        ("show rvStableSchurSpecialize[0](rvStableSchurCharacterScale[-170141183460469231731687303715884105728](a));", "signed 128-bit range"),
    ] {
        let error = diagnostic(interpret(&format!("{SYMS} {query}")));
        assert!(error.contains(expected),"{error}");
    }
}
