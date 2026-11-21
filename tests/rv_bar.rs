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

const FLAG_SQUARE: &str = "
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

#[test]
fn binary_and_ternary_replacements_cancel_in_the_length_one_part_of_d_squared() {
    assert!(stdout(interpret(FLAG_SQUARE)).is_empty());
}

#[test]
fn the_first_staircase_page_differential_is_d_n_minus_one_on_both_actual_classes() {
    for n in [3, 4, 6] {
        let mut source = staircase(n);
        for r in 0..n {
            writeln!(source, "let e{r} = rvBarClass[{r},{n}](raw);").unwrap();
            writeln!(
                source,
                "assert e{r} == rvbaregen[P,0,{n},{r},{},{n},0];",
                2 * n
            )
            .unwrap();
            writeln!(source, "assert rvBarRep(e{r}) == raw;").unwrap();
            if r + 1 < n {
                writeln!(
                    source,
                    "assert rvBarD(e{r}) == zero[RVBarE(P,0,{n},{r},{},{})];",
                    2 * n - 1,
                    n - r
                )
                .unwrap();
            }
            if r > 0 {
                writeln!(source, "assert rvBarPrev(e{r}) == e{};", r - 1).unwrap();
                writeln!(source, "assert rvBarNext(e{}) == e{r};", r - 1).unwrap();
            }
        }
        writeln!(source, "let target = rvBarClass[{},1](top);", n - 1).unwrap();
        writeln!(
            source,
            "assert target == rvbaregen[P,0,{n},{},{},1,0];",
            n - 1,
            2 * n - 1
        )
        .unwrap();
        writeln!(source, "assert rvBarD(e{}) == target;", n - 1).unwrap();
        writeln!(
            source,
            "assert rvBarD(target) == zero[RVBarE(P,0,{n},{},{},{})];",
            n - 1,
            2 * n - 2,
            2isize - n as isize
        )
        .unwrap();
        writeln!(
            source,
            "assert rvBarNext(target) == zero[RVBarE(P,0,{n},{n},{},1)];",
            2 * n - 1
        )
        .unwrap();
        writeln!(source, "rvbaredim[P,0,{n},{n},{},{n}];", 2 * n).unwrap();
        writeln!(source, "rvbaredim[P,0,{n},{n},{},1];", 2 * n - 1).unwrap();
        assert!(stdout(interpret(&source)).ends_with(&format!(
            "rvbaredim[P, 0, {n}, {n}, {}, {n}] = 0\nrvbaredim[P, 0, {n}, {n}, {}, 1] = 0\n",
            2 * n,
            2 * n - 1
        )));
    }
}

#[test]
fn page_projection_uses_the_lower_length_correction_and_checks_cycles() {
    let source = format!(
        "{FLAG_SQUARE}
        assert rvBarClass[3,3](dx) == zero[RVBarE(P,0,4,3,7,3)];
        assert rvBarClass[2,3](binary) == zero[RVBarE(P,0,4,2,7,3)];
        "
    );
    assert!(stdout(interpret(&source)).is_empty());
    let error = diagnostic(interpret(&format!(
        "{FLAG_SQUARE} show rvBarClass[3,3](binary);"
    )));
    assert!(error.contains("page cycle condition"), "{error}");
    let error = diagnostic(interpret(&format!(
        "{FLAG_SQUARE} show rvBarClass[1,2](dx);"
    )));
    assert!(error.contains("above the requested"), "{error}");
}

#[test]
fn the_first_binary_page_checks_next_and_previous_types_and_enumerates_classes() {
    let squares = include_str!("../examples/rank-bar.asl");
    let error = diagnostic(interpret(&format!("{squares} show rvBarNext(e1);")));
    assert!(error.contains("requires a cycle for d_1"), "{error}");
    let error = diagnostic(interpret(&format!("{squares} show rvBarPrev(e0);")));
    assert!(error.contains("positive page number"), "{error}");
    let result = stdout(interpret(&format!(
        r"{squares}
        card \x:RVBarE(P,0,3,1,6,3)=>S.hi;
        rvbaredim[P,0,3,1,6,-1];
        rvbaredim[P,0,3,1,6,4];
        "
    )));
    assert!(
        result.ends_with(
            "card = 2 / 2\nrvbaredim[P, 0, 3, 1, 6, -1] = 0\nrvbaredim[P, 0, 3, 1, 6, 4] = 0\n"
        ),
        "{result}"
    );
}

#[test]
fn dual_differential_is_the_transpose_and_the_cut_product_obeys_leibniz() {
    assert_eq!(
        stdout(example("rank-path-ext")),
        concat!(
            "rvpathextdim[P, 0, 1, 2] = 1\n",
            "rvpathextdim[P, 0, 3, 4] = 0\n",
            "rvpathextdim[P, 0, 3, 5] = 0\n",
            "rvpathextdim[P, 0, 3, 6] = 0\n",
            "RVBarDual(P, 0, 3, 5){[(0, 1, 1):g0 | (1, 3, 2):g0] + [(0, 2, 2):g0 | (2, 3, 1):g0]}\n"
        )
    );
}

#[test]
fn nonzero_path_ext_products_are_associative_and_expand_class_sums() {
    let input = format!(
        r"{INDEPENDENT}
        rvbarpath P[J] {{[1,1,1],[1,1,1],[1,1,1]}};
        let u=rvmtorgen[J,2,[1,1,1],0];
        let v=rvmtorgen[J,2,[1,1,1],1];
        let raw_a=rvBarWord[P,[0,1]](u);
        let a=rvPathClass(rvBarDualize(raw_a));
        let b=rvpathextgen[P,0,1,3,1];
        let c=rvpathextgen[P,1,2,3,0];
        let d=rvpathextgen[P,2,3,3,0];
        assert rvPathRep(a)==rvBarDualize(raw_a);
        let left=rvPathCup(rvPathCup(add(a,b), c), d);
        assert left==rvPathCup(add(a,b),rvPathCup(c,d));
        assert left==add(rvPathCup(a,rvPathCup(c,d)),rvPathCup(b,rvPathCup(c,d)));
        assert left==add(rvpathextgen[P,0,3,9,0],rvpathextgen[P,0,3,9,4]);
        assert rvPathClass(rvPathRep(left))==left;
        card \s:RVBarDual(P,0,2,6)=>S.hi;
        rvpathextdim[P,0,3,9];
        "
    );
    assert_eq!(
        stdout(interpret(&input)),
        "card = 16 / 16\nrvpathextdim[P, 0, 3, 9] = 8\n"
    );
}

#[test]
fn staircase_proper_cups_survive_and_the_full_cup_has_the_prescribed_primitive() {
    for n in [3, 4, 6] {
        let mut source = staircase(n);
        for j in 1..=n {
            writeln!(
                source,
                "let dc{j} = rvBarDualize(rvBarWord[P,[{},{j}]](a{j}));",
                j - 1
            )
            .unwrap();
            writeln!(source, "let ext{j}_{j} = rvPathClass(dc{j});").unwrap();
            writeln!(
                source,
                "assert ext{j}_{j} == rvpathextgen[P,{},{j},2,0];",
                j - 1
            )
            .unwrap();
        }
        writeln!(source, "let cup1 = dc1;").unwrap();
        for j in 2..=n {
            writeln!(source, "let cup{j} = rvPathCup(cup{},dc{j});", j - 1).unwrap();
        }
        for length in 2..=n {
            for i in 1..=n - length + 1 {
                let j = i + length - 1;
                writeln!(
                    source,
                    "let ext{i}_{j}=rvPathCup(ext{i}_{},ext{j}_{j});",
                    j - 1
                )
                .unwrap();
                if length < n {
                    writeln!(
                        source,
                        "assert ext{i}_{j} == rvpathextgen[P,{},{j},{},0];",
                        i - 1,
                        2 * length
                    )
                    .unwrap();
                }
            }
        }
        writeln!(source, "let primitive = rvBarDualize(top);").unwrap();
        writeln!(source, "assert rvBarDualD(primitive)==cup{n};").unwrap();
        writeln!(source, "assert cup{n}==rvBarDualize(raw);").unwrap();
        writeln!(source, "assert rvBarPair(cup{n},raw)==S.hi;").unwrap();
        writeln!(source, "assert rvPathClass(cup{n})==ext1_{n};").unwrap();
        writeln!(
            source,
            "assert ext1_{n}==zero[RVPathExt(P,0,{n},{})];",
            2 * n
        )
        .unwrap();
        writeln!(source, "rvpathextdim[P,0,{n},{}];", 2 * n).unwrap();
        assert!(stdout(interpret(&source))
            .ends_with(&format!("rvpathextdim[P, 0, {n}, {}] = 0\n", 2 * n)));
        let error = diagnostic(interpret(&format!("{source} show rvPathClass(primitive);")));
        assert!(
            error.contains("requires a closed dual bar cochain"),
            "{error}"
        );
    }
}

#[test]
fn duality_products_and_cohomology_validate_compatible_variances_and_degrees() {
    let source = include_str!("../examples/rank-path-ext.asl");
    for (expr, expected) in [
        ("rvBarPair(word,word)", "requires an RVBarDual"),
        ("rvBarPair(first,left)", "RVBar(P, 0, 1, 2)"),
        ("rvPathCup(first,first)", "sharing a cut"),
        (
            "rvPathCup(first,rvPathClass(first))",
            "two RVBarDual cochains or two RVPathExt classes",
        ),
        ("rvPathRep(top)", "requires an RVPathExt"),
        ("rvBarDualD(word)", "requires an RVBarDual"),
        ("rvPathClass(f)", "requires a closed dual bar cochain"),
    ] {
        let error = diagnostic(interpret(&format!("{source} show {expr};")));
        assert!(error.contains(expected), "{expr}: {error}");
    }
    let error = diagnostic(interpret(&format!(
        "{source} show rvBarDualD(zero[RVBarDual(P,0,1,{})]);",
        usize::MAX
    )));
    assert!(error.contains("degree is out of range"), "{error}");
}
