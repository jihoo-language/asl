mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn marking() -> &'static str {
    include_str!("../examples/hyperelliptic-ogrady-sixfold.asl")
        .split_once("let theta")
        .unwrap()
        .0
}

const STANDARD: &str = r"
    let w=rvAbMukai[B,1,-1]{0,0,0,0,0,0};
    rvog6 M(w);
";

#[test]
fn example_checks_the_two_exceptional_normalizations_and_intersections() {
    assert_eq!(
        stdout(example("hyperelliptic-ogrady-sixfold")),
        concat!(
            "rvAbMukaiPair = 2\n",
            "rvK3LatticePair = 4\n",
            "rvOG6Lattice[Six] = n=3, dimension=6, rank=8, signature=(3,5), determinant=-4, discriminant=(Z/2)^2, q=diag(-1/2,-1/2) mod 2\n",
            "rvOG6Pair = -2\n",
            "rvOG6Pair = -8\n",
            "rvOG6Div = 2\n",
            "rvOG6Fujiki[0] = -480\n",
            "rvOG6Fujiki[0] = -30720\n",
            "rvOG6Fujiki[2] = -96\n",
            "rvOG6RR = 4\n",
            "rvOG6RR = 0\n",
            "rvOG6RR = -4\n",
        )
    );
}

#[test]
fn all_abelian_basis_pairings_and_middle_embedding_use_the_marked_orientation() {
    let mut source = marking().to_owned();
    let mut expected = String::new();
    for i in 0..8 {
        writeln!(
            source,
            "let x{i}=rvAbMukaiBasis[B,{i}]; assert add(x{i},neg(x{i}))==zero[RVAbMukai(B)];"
        )
        .unwrap();
    }
    let signs = [-1, -1, 1, -1, -1, 1, -1, -1];
    for (i, &sign) in signs.iter().enumerate() {
        for j in i..8 {
            let value = if i + j == 7 { sign } else { 0 };
            writeln!(source, "rvAbMukaiPair(x{i},x{j});").unwrap();
            writeln!(expected, "rvAbMukaiPair = {value}").unwrap();
            if (1..=6).contains(&i) && (1..=6).contains(&j) {
                writeln!(
                    source,
                    "rvK3LatticePair(rvAbMukaiKummerH2(x{i}),rvAbMukaiKummerH2(x{j}));"
                )
                .unwrap();
                writeln!(expected, "rvK3LatticePair = {}", 2 * value).unwrap();
            }
        }
    }
    source.push_str(
        r"
        assert rvAbMukaiKummerH2(add(x0,x7))==zero[RVK3Lattice(B)];
        assert rvAbMukaiKummerH2(rvAbMukai[B,5,-1]{1,-2,3,-4,5,-6})
            ==rvK3Lattice[B]{[2,-4,6,-8,10,-12],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
        ",
    );
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn independent_integral_gram_determinant_for_a_nontrivial_mukai_vector() {
    // The perpendicular covector of w has gcd 1 without a unit entry.
    let mut source = format!(
        r"{}
        let w=rvAbMukai[B,2,3]{{2,3,0,0,5,4}};
        rvAbMukaiPair(w,w);
        rvog6 M(w);
        let e=rvOG6ExcHalf[M];
        let exc=rvOG6Exc[M];
        rvOG6Div(e);
        rvOG6Div(exc);
        let u=rvAbMukai[B,0,0]{{0,0,1,0,0,0}};
        let u2=rvAbMukai[B,0,0]{{0,0,0,-1,0,0}};
        let x=rvOG6[M,0](add(u,u2));
        rvOG6Pair(x,e);
        rvOG6RR(x);
        rvOG6RR(add(e,exc));
        ",
        marking()
    );
    for i in 0..8 {
        writeln!(source, "let a{i}=rvOG6Basis[M,{i}];").unwrap();
    }
    for i in 0..8 {
        for j in i..8 {
            writeln!(source, "rvOG6Pair(a{i},a{j});").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut numbers = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1.parse::<i128>().unwrap());
    assert_eq!(
        numbers.by_ref().take(6).collect::<Vec<_>>(),
        [2, 2, 4, 0, 16, -224]
    );

    let mut upper = [[0i128; 8]; 8];
    for (i, row) in upper.iter_mut().enumerate() {
        for entry in row.iter_mut().skip(i) {
            *entry = numbers.next().unwrap();
        }
        assert_eq!(row[i] % 2, 0);
    }
    assert_eq!(numbers.next(), None);
    let gram: Vec<Vec<_>> = (0..8)
        .map(|i| (0..8).map(|j| upper[i.min(j)][i.max(j)]).collect())
        .collect();
    assert_eq!(gram[7], [0, 0, 0, 0, 0, 0, 0, -2]);
    assert_eq!(determinant(&gram), -4);
}

fn determinant(matrix: &[Vec<i128>]) -> i128 {
    let Some(first) = matrix.first() else {
        return 1;
    };
    first
        .iter()
        .enumerate()
        .filter(|(_, entry)| **entry != 0)
        .map(|(column, &entry)| {
            let minor: Vec<Vec<_>> = matrix[1..]
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter_map(|(index, &x)| (index != column).then_some(x))
                        .collect()
                })
                .collect();
            let sign = if column % 2 == 0 { 1 } else { -1 };
            sign * entry * determinant(&minor)
        })
        .sum()
}

#[test]
fn exact_wide_cancellation_keeps_valid_large_classes() {
    let t = 1i128 << 126;
    assert_eq!(
        stdout(interpret(&format!(
            r"{}{STANDARD}
            let u=rvAbMukai[B,0,0]{{0,{t},0,0,{t},0}};
            let x=rvOG6[M,{t}](u);
            assert add(x,neg(x))==zero[RVOG6(M)];
            rvOG6Pair(x,x);
            rvOG6Pair(x,neg(x));
            rvOG6Fujiki[2](x,x);
            rvOG6RR(x);
            ",
            marking()
        ))),
        "rvOG6Pair = 0\nrvOG6Pair = 0\nrvOG6Fujiki[2] = 0\nrvOG6RR = 4\n"
    );
}

#[test]
fn invalid_types_norms_and_arithmetic_ranges_are_diagnosed() {
    let setup = format!("{}{STANDARD}", marking());
    for (source, expected) in [
        ("rvog6 Bad(add(w,w));", "requires a primitive abelian Mukai vector of square two"),
        ("rvog6 Bad(rvMukaiIdeal[B,2]);", "an RVAbMukai vector"),
        ("rvog6 M(w);", "is already defined"),
        ("rvOG6Lattice[Missing];", "unknown OG6 lattice"),
        ("show rvOG6Basis[M,8];", "between 0 and 7"),
        ("show rvAbMukaiBasis[B,8];", "between 0 and 7"),
        ("show rvAbMukai[B,0,0]{0,0,0,0,0};", "six wedge coordinates"),
        ("show rvOG6[M,0](w);", "perpendicular to the defining vector"),
        ("rvog6 Twin(w); rvOG6Pair(rvOG6Exc[M],rvOG6Exc[Twin]);", "expected type RVOG6(M)"),
        ("rvmukaimoduli K(rvMukaiIdeal[B,2]); rvOG6Pair(rvOG6Exc[M],rvOG10Exc[K]);", "expected type RVOG6(M)"),
        ("rvjactorsionbasis Twin[C,2](neg(o)){e1,e2,f1,f2}; show rvOG6[M,0](zero[RVAbMukai(Twin)]);", "expected type RVAbMukai(B)"),
        ("rvjactorsionbasis Twin[C,2](neg(o)){e1,e2,f1,f2}; rvAbMukaiPair(w,zero[RVAbMukai(Twin)]);", "expected type RVAbMukai(B)"),
        ("show rvAbMukaiKummerH2(rvAbMukai[B,0,0]{0,85070591730234615865843651857942052864,0,0,0,0});", "Kummer wedge coordinate exceeds"),
        ("rvOG6Fujiki[7](rvOG6ExcHalf[M],rvOG6ExcHalf[M]);", "between zero and 2n"),
        ("show neg(rvOG6[M,-170141183460469231731687303715884105728](zero[RVAbMukai(B)]));", "half-exceptional coefficient exceeds"),
        ("let x=rvOG6[M,85070591730234615865843651857942052864](zero[RVAbMukai(B)]); rvOG6Pair(x,x);", "OG6 pairing exceeds"),
        ("let x=rvOG6[M,85070591730234615865843651857942052864](zero[RVAbMukai(B)]); rvOG6Div(x);", "OG6 divisibility exceeds"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
