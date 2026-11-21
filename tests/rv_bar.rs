mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn both_nonzero_binary_contractions_cancel_in_the_square_of_the_bar_differential() {
    assert_eq!(
        stdout(example("rank-bar")),
        concat!(
            "rvbardim[P, 0, 3, 6] = 1\n",
            "rvbardim[P, 0, 3, 5, 2] = 2\n",
            "rvbardim[P, 0, 3, 4] = 1\n",
            "RVBar(P, 0, 3, 5){[(0, 1, 1):g0 | (1, 3, 2):g0] + [(0, 2, 2):g0 | (2, 3, 1):g0]}\n"
        )
    );
}

#[test]
fn a_staircase_has_only_two_bar_words_and_its_differential_is_the_top_operation() {
    for n in [3, 4, 6] {
        let program = staircase(n);
        let lines = (0..=3 * n)
            .map(|b| {
                let rank = usize::from(b == 2 * n || b == 2 * n - 1);
                format!("rvbardim[P, 0, {n}, {b}] = {rank}\n")
            })
            .collect::<String>();
        assert_eq!(stdout(interpret(&program)), lines, "arity {n}");
    }
}

#[test]
fn multilinearity_expands_class_supports_and_respects_zero() {
    let source = format!(
        r"{INDEPENDENT}
        rvbarpath P[J] {{ [1,1,1], [1,1,1] }};
        let a = rvmtorgen[J,2,[1,1,1],0];
        let b = rvmtorgen[J,2,[1,1,1],1];
        let z = zero[RVMTor(J,2,[1,1,1])];
        let x = rvBarWord[P,[0,1,2]](add(a,b),add(a,b));
        let aa = rvBarWord[P,[0,1,2]](a,a);
        let ab = rvBarWord[P,[0,1,2]](a,b);
        let ba = rvBarWord[P,[0,1,2]](b,a);
        let bb = rvBarWord[P,[0,1,2]](b,b);
        assert x == add(add(aa,ab),add(ba,bb));
        assert rvBarWord[P,[0,1,2]](a,z) == zero[RVBar(P,0,2,6)];
        assert rvBarD(x) == zero[RVBar(P,0,2,5)];
        card \x : RVBar(P,0,2,6) => S.hi;
        "
    );
    assert_eq!(stdout(interpret(&source)), "card = 16 / 16\n");
}

#[test]
fn sparse_path_dp_counts_large_spaces_without_tabulation() {
    let steps = vec!["[1,1,1]"; 20].join(",");
    let source = format!(
        "{INDEPENDENT} rvbarpath P[J] {{{steps}}};
         rvbardim[P,0,20,60];
         rvbardim[P,0,20,60,20];
         rvbardim[P,0,20,60,19];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvbardim[P, 0, 20, 60] = 1048576\n",
            "rvbardim[P, 0, 20, 60, 20] = 1048576\n",
            "rvbardim[P, 0, 20, 60, 19] = 0\n"
        )
    );
    let error = diagnostic(interpret(&format!(
        "{source} card \\x:RVBar(P,0,20,60)=>S.hi;"
    )));
    assert!(error.contains("enumeration limit"), "{error}");
}

const INDEPENDENT: &str = "
    rvideal J[[0,1,2]] {
        rvmon[[0,1,2]]{1,1,0}, rvmon[[0,1,2]]{1,0,1}, rvmon[[0,1,2]]{0,1,1}
    };
";

#[test]
fn bar_declarations_and_words_validate_paths_intervals_and_classes() {
    let base = format!(
        "{INDEPENDENT} rvbarpath P[J] {{[1,1,1],[1,1,1]}};
        let a = rvmtorgen[J,2,[1,1,1],0];"
    );
    for (source, expected) in [
        ("rvbarpath Q[J] {[0,0,0]};", "nonzero multidegree"),
        ("rvbarpath Q[J] {[1,1]};", "one exponent"),
        ("rvbarpath P[J] {[1,1,1]};", "already defined"),
        ("rvbardim[Unknown,0,1,2];", "unknown bar path"),
        ("rvbardim[P,1,1,2];", "invalid ordered bar interval"),
        (
            "show rvBarWord[P,[0,1,2]](a);",
            "one class per cut interval",
        ),
        (
            "show rvBarWord[P,[1,0]](a);",
            "invalid ordered bar interval",
        ),
        ("show rvBarWord[P,[0,2]](a);", "wrong ideal or multidegree"),
        ("show rvBarWord[P,[0,1]](One.pt);", "requires RVMTor"),
        ("show rvBarD(zero[RVBar(P,0,1,0)]);", "positive suspended"),
    ] {
        let error = diagnostic(interpret(&format!("{base}\n{source}")));
        assert!(error.contains(expected), "{source}: {error}");
    }
}

fn staircase(n: usize) -> String {
    let indices: Vec<_> = (1..=2 * n).collect();
    let powers = |i, j| -> Vec<_> {
        indices
            .iter()
            .map(|&k| usize::from(k == i || k == j))
            .collect()
    };
    let monomial = |i, j| {
        format!(
            "rvmon[{indices:?}]{{{}}}",
            powers(i, j)
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        )
    };
    let ideal = (1..=n)
        .flat_map(|j| {
            (1..=j)
                .filter(move |&k| (j, k) != (n, 1))
                .map(move |k| monomial(j, n + k))
        })
        .collect::<Vec<_>>()
        .join(", ");
    let mut source = format!("rvideal Stair[{indices:?}] {{{ideal}}};\n");
    let steps = (1..=n)
        .map(|j| format!("{:?}", powers(j, n + j)))
        .collect::<Vec<_>>()
        .join(",");
    writeln!(source, "rvbarpath P[Stair] {{{steps}}};").unwrap();
    for j in 1..=n {
        writeln!(
            source,
            "let a{j} = rvmtorgen[Stair,1,{:?},0];",
            powers(j, n + j)
        )
        .unwrap();
    }
    let args = (1..=n)
        .map(|j| format!("a{j}"))
        .collect::<Vec<_>>()
        .join(",");
    let cuts: Vec<_> = (0..=n).collect();
    let wedge = indices
        .iter()
        .copied()
        .filter(|&k| k != n && k != n + 1)
        .collect::<Vec<_>>();
    writeln!(source, "let raw = rvBarWord[P,{cuts:?}]({args});").unwrap();
    writeln!(
        source,
        "let expected = rvTorClass(rvmkterm[Stair,{wedge:?}]({}));",
        monomial(n, n + 1)
    )
    .unwrap();
    writeln!(source, "let top = rvBarWord[P,[0,{n}]](expected);").unwrap();
    writeln!(source, "assert rvBarD(raw) == top;").unwrap();
    writeln!(source, "assert rvBarLength[1](rvBarD(raw)) == top;").unwrap();
    writeln!(
        source,
        "assert rvBarD(top) == zero[RVBar(P,0,{n},{})];",
        2 * n - 2
    )
    .unwrap();
    for b in 0..=3 * n {
        writeln!(source, "rvbardim[P,0,{n},{b}];").unwrap();
    }
    source
}

#[test]
fn squarefree_interval_ranks_keep_true_higher_nonfaces_and_forbidden_vertices() {
    let source = "
        rvideal Ring[[0,1,2,3]] {rvmon[[0,1,2,3]]{1,1,1,0},rvmon[[0,1,2,3]]{0,0,0,1}};
        rvbarpath P[Ring] {[1,1,1,1]};
        rvbardim[P,0,1,3];
        rvbardim[P,0,1,4];
        assert rvBarWord[P,[0,1]](rvmtorgen[Ring,2,[1,1,1,1],0]) ==
               rvBarLength[1](rvBarWord[P,[0,1]](rvmtorgen[Ring,2,[1,1,1,1],0]));
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvbardim[P, 0, 1, 3] = 1\nrvbardim[P, 0, 1, 4] = 0\n"
    );
}

#[test]
fn overflowing_suffix_counts_do_not_obscure_a_zero_degree() {
    let steps = vec!["[1,1,1]"; 130].join(",");
    let source = format!(
        "{INDEPENDENT} rvbarpath P[J] {{{steps}}};
         rvbardim[P,0,130,391];"
    );
    assert_eq!(stdout(interpret(&source)), "rvbardim[P, 0, 130, 391] = 0\n");
    let error = diagnostic(interpret(&format!("{source} rvbardim[P,0,130,390];")));
    assert!(error.contains("dimension is out of range"), "{error}");
}

#[test]
fn binary_and_ternary_replacements_cancel_in_the_length_one_part_of_d_squared() {
    let source = "
        rvideal FlagSquare[[1,2,3,4,5,6,7]] {
            rvmon[[1,2,3,4,5,6,7]]{1,1,0,0,0,0,0},
            rvmon[[1,2,3,4,5,6,7]]{0,0,1,1,0,0,0},
            rvmon[[1,2,3,4,5,6,7]]{0,0,0,0,1,1,0},
            rvmon[[1,2,3,4,5,6,7]]{1,0,1,0,0,0,0},
            rvmon[[1,2,3,4,5,6,7]]{0,0,0,1,1,0,0},
            rvmon[[1,2,3,4,5,6,7]]{0,0,0,0,0,0,2}
        };
        rvbarpath P[FlagSquare] {
            [1,1,0,0,0,0,0], [0,0,0,0,1,1,0], [0,0,1,1,0,0,0], [0,0,0,0,0,0,2]
        };
        let a = rvmtorgen[FlagSquare,1,[1,1,0,0,0,0,0],0];
        let c = rvmtorgen[FlagSquare,1,[0,0,0,0,1,1,0],0];
        let b = rvmtorgen[FlagSquare,1,[0,0,1,1,0,0,0],0];
        let d = rvmtorgen[FlagSquare,1,[0,0,0,0,0,0,2],0];
        let x = rvBarWord[P,[0,1,2,3,4]](a,c,b,d);
        let dx = rvBarD(x);
        let ternary = rvBarLength[2](dx);
        let binary = rvBarLength[3](dx);
        let surviving = rvmtorgen[FlagSquare,5,[1,1,1,1,1,1,2],0];
        let endpoint = rvBarWord[P,[0,4]](surviving);
        assert rvBarLength[1](rvBarD(ternary)) == endpoint;
        assert rvBarLength[1](rvBarD(binary)) == endpoint;
        assert rvBarD(dx) == zero[RVBar(P,0,4,6)];
        ";
    assert!(stdout(interpret(source)).is_empty());
}
