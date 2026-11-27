mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn higher_word_ranks_and_simultaneous_invariants_separate_the_genus_five_examples() {
    assert_eq!(
        stdout(example("rank-as-words")),
        concat!(
            "rvASWordRank[A, [F, F, F]] = 2\n",
            "rvASWordRank[B, [F, F, F]] = 1\n",
            "rvASWordRank[A, [V, V, V]] = 2\n",
            "rvASWordRank[B, [V, V, V]] = 1\n",
            "rvASCommonKer[A, [[F, F], [V]]] = 3\n",
            "rvASCommonKer[B, [[F, F], [V]]] = 4\n",
            "rvASSumImages[A, [[V, V], [F]]] = 7\n",
            "rvASSumImages[B, [[V, V], [F]]] = 6\n",
        )
    );
}

#[test]
fn words_agree_with_a_number_p_rank_and_composition() {
    let setup = include_str!("../examples/rank-as-words.asl");
    let rules = r"
        rvASaNumber[A];
        rvASPrank[A];
        rvASWordRank[A,[]];
        rvASWordRank[A,[F,V,F]];
        rvASWordRank[A,[V,F]];
        rvASWordRank[A,[F,F,F,F,F,F,F,F,F,F,F,F]];
        rvASWordRank[A,[V,V,V,V,V,V,V,V,V,V,V,V]];
        rvASCommonKer[A,[[F],[V],[F,F],[V,F]]];
        rvASSumImages[A,[[F],[V],[V,V],[F,V]]];
        rvASCommonKer[A,[]];
        rvASSumImages[A,[]];
        rvASCommonKer[A,[[F],[],[V]]];
        rvASSumImages[A,[[F],[],[V]]];
    ";
    let result = stdout(interpret(&format!("{setup}\n{rules}")));
    assert!(
        result.ends_with(concat!(
            "rvASaNumber[A] = 2\n",
            "rvASPrank[A] = 1\n",
            "rvASWordRank[A, []] = 10\n",
            "rvASWordRank[A, [F, V, F]] = 0\n",
            "rvASWordRank[A, [V, F]] = 0\n",
            "rvASWordRank[A, [F, F, F, F, F, F, F, F, F, F, F, F]] = 1\n",
            "rvASWordRank[A, [V, V, V, V, V, V, V, V, V, V, V, V]] = 1\n",
            "rvASCommonKer[A, [[F], [V], [F, F], [V, F]]] = 2\n",
            "rvASSumImages[A, [[F], [V], [V, V], [F, V]]] = 8\n",
            "rvASCommonKer[A, []] = 10\n",
            "rvASSumImages[A, []] = 0\n",
            "rvASCommonKer[A, [[F], [], [V]]] = 0\n",
            "rvASSumImages[A, [[F], [], [V]]] = 10\n",
        )),
        "{result}"
    );
}

#[test]
fn malformed_words_are_diagnosed_and_trivial_words_need_no_matrix() {
    let order: usize = (1usize << (usize::BITS / 2)) + 1;
    let setup =
        format!("rvascover Huge[[0,1]] {{rvmon[[0,1]]{{{order},0}},rvmon[[0,1]]{{0,{order}}}}};");
    assert_eq!(
        stdout(interpret(&format!(
            "{setup} rvASWordRank[Huge,[]]; rvASWordRank[Huge,[F,V]];"
        ))),
        format!(
            "rvASWordRank[Huge, []] = {}\nrvASWordRank[Huge, [F, V]] = 0\n",
            order - 1
        ),
    );
    for (statement, reason) in [
        ("rvASWordRank[Huge,[F]];", "storage limit"),
        ("rvASWordRank[Huge,[Q]];", "expected F or V"),
        (
            "rvASCommonKer[Missing,[[F],[V]]];",
            "unknown Artin–Schreier cover",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
