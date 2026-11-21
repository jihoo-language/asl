mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn p_adic_and_invariant_ideal_examples_detect_different_pages() {
    assert_eq!(
        stdout(example("bp-pages")),
        concat!(
            "bpBockdim[2, 1, 3, 0, 8, 0] = 1\n",
            "bpBockdim[2, 1, 4, 0, 8, 0] = 0\n",
            "bpBockdim[2, 2, 4, 0, 8, 0] = 1\n",
            "bpBockdim[2, 2, 5, 0, 8, 0] = 0\n",
            "bpNovdim[2, 2, 2, 1, 16, 1] = 1\n",
            "bpNovdim[2, 2, 2, 2, 16, 3] = 1\n",
            "bpNovdim[2, 2, 3, 1, 16, 1] = 0\n",
            "bpNovdim[2, 2, 3, 2, 16, 3] = 0\n",
        )
    );
}

#[test]
fn dimensions_match_exact_filtered_lattices_at_two_and_three() {
    let mut program = String::new();
    let mut output = String::new();
    for (cmd, p, n, a, t, f, ranks) in [
        ("bpBockdim", 2, 2, 0, 8, 0, [2, 1, 1, 1, 1, 0]),
        ("bpBockdim", 2, 2, 1, 8, 4, [7, 3, 1, 1, 1, 0]),
        ("bpBockdim", 2, 1, 1, 8, 3, [4, 3, 1, 1, 0, 0]),
        ("bpNovdim", 2, 2, 1, 16, 1, [8, 1, 1, 0, 0, 0]),
        ("bpNovdim", 2, 2, 2, 16, 3, [78, 1, 1, 0, 0, 0]),
        ("bpNovdim", 2, 2, 1, 16, 2, [14, 1, 1, 0, 0, 0]),
        ("bpNovdim", 2, 2, 2, 16, 4, [84, 1, 1, 0, 0, 0]),
        ("bpBockdim", 3, 2, 0, 12, 0, [1, 1, 1, 0, 0, 0]),
        ("bpBockdim", 3, 2, 1, 12, 2, [3, 2, 1, 0, 0, 0]),
        ("bpNovdim", 3, 2, 0, 12, 3, [1, 0, 0, 0, 0, 0]),
    ] {
        for (r, rank) in ranks.into_iter().enumerate() {
            let query = format!("{cmd}[{p}, {n}, {r}, {a}, {t}, {f}]");
            writeln!(program, "{query};").unwrap();
            writeln!(output, "{query} = {rank}").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&program)), output);
}

#[test]
fn dimensions_check_bounds_and_remove_common_filtration_powers() {
    let max = usize::MAX;
    let program = format!("bpBockdim[3, 0, {max}, 0, 0, {max}];\nbpNovdim[3, 0, 1, 1, 0, {max}];");
    assert_eq!(
        stdout(interpret(&program)),
        format!("bpBockdim[3, 0, {max}, 0, 0, {max}] = 1\nbpNovdim[3, 0, 1, 1, 0, {max}] = 0\n")
    );
    assert!(diagnostic(interpret("bpNovdim[4, 2, 0, 0, 0, 0];")).contains("requires a prime"));
    assert!(diagnostic(interpret("bpNovdim[3, 9, 0, 0, 0, 0];")).contains("height exceeds"));
}
