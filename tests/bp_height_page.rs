mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn integral_height_differential_and_extension_are_visible() {
    assert_eq!(
        stdout(example("bp-height-pages")),
        concat!(
            "bpHeightPage[2, 2, 0, 2, 1, 6] = free 1, torsion []\n",
            "bpHeightImage[2, 2, 0, 2, 0, 6] = free 1, order infinite\n",
            "bpHeightPage[2, 2, 1, 2, 1, 6] = free 0, torsion [2]\n",
            "bpHeightPage[2, 2, 1, 1, 2, 6] = free 0, torsion [2]\n",
            "bpHeightImage[2, 2, 1, 2, 1, 6] = order 2\n",
            "bpHeightPage[2, 2, 2, 2, 1, 6] = free 0, torsion []\n",
            "bpHeightPage[2, 2, 2, 1, 2, 6] = free 0, torsion []\n",
            "bpHeightPage[2, 2, 3, 1, 1, 8] = free 0, torsion [8]\n",
            "bpHeightPage[2, 2, 3, 2, 1, 8] = free 0, torsion [2]\n",
            "bpext[2, 2, 1, 8] = free 0, torsion [16]\n",
        )
    );
}

#[test]
fn very_late_pages_keep_boundaries_from_above_the_ambient_height() {
    let late = usize::MAX;
    let source = format!(
        "bpHeightPage[2,2,{late},1,2,6];\n\
         bpHeightPage[2,2,{late},2,1,8];\n\
         bpHeightImage[2,2,{late},2,1,8];\n\
         bpHeightPage[2,0,{late},0,0,0];\n\
         bpHeightPage[2,0,{late},0,1,0];\n\
         bpHeightPage[2,6,{late},5,1,0];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "bpHeightPage[2, 2, {late}, 1, 2, 6] = free 0, torsion []\n\
             bpHeightPage[2, 2, {late}, 2, 1, 8] = free 0, torsion [2]\n\
             bpHeightImage[2, 2, {late}, 2, 1, 8] = order 1\n\
             bpHeightPage[2, 0, {late}, 0, 0, 0] = free 1, torsion []\n\
             bpHeightPage[2, 0, {late}, 0, 1, 0] = free 0, torsion []\n\
             bpHeightPage[2, 6, {late}, 5, 1, 0] = free 0, torsion []\n"
        ),
    );
}

#[test]
fn page_one_agrees_with_relative_groups_at_both_primes() {
    let source = r"
        bpHeightPage[3,2,1,2,1,16];
        bprelext[3,1,2,1,16];
        bpHeightImage[3,2,1,2,1,16];
        bpHeightPage[3,2,2,1,2,16];
        bpHeightPage[3,2,2,2,1,16];
        bpHeightPage[2,3,1,1,1,8];
        bpext[2,1,1,8];
        bpHeightPage[2,2,2,1,2,8];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "bpHeightPage[3, 2, 1, 2, 1, 16] = free 0, torsion [3]\n",
            "bprelext[3, 1, 2, 1, 16] = free 0, torsion [3]\n",
            "bpHeightImage[3, 2, 1, 2, 1, 16] = order 3\n",
            "bpHeightPage[3, 2, 2, 1, 2, 16] = free 0, torsion []\n",
            "bpHeightPage[3, 2, 2, 2, 1, 16] = free 0, torsion []\n",
            "bpHeightPage[2, 3, 1, 1, 1, 8] = free 0, torsion [8]\n",
            "bpext[2, 1, 1, 8] = free 0, torsion [8]\n",
            "bpHeightPage[2, 2, 2, 1, 2, 8] = free 0, torsion [2, 2]\n",
        ),
    );
}

#[test]
fn height_page_statements_validate_the_ambient_grade() {
    for (source, reason) in [
        (
            "bpHeightPage[2,1,1,2,1,6];",
            "BP height filtration exceeds the ambient height",
        ),
        ("bpHeightImage[4,1,1,1,1,6];", "BP requires a prime"),
        (
            "bpHeightPage[2,3,1,1,17,6];",
            "BP arrow arity exceeds the limit",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
