mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn transverse_cycles_are_hit_on_distinct_actual_pages() {
    let output = stdout(example("annular-khovanov-pages"));
    for line in [
        "khPageDim[Cross, 3, 1, 0, 3, 2] = 1",
        "khPageDim[Cross, 3, 2, 0, 3, 2] = 0",
        "khPageDim[Trefoil, 2, 1, 2, 7, 0] = 1",
        "khPageDim[Trefoil, 2, 1, 3, 7, 0] = 1",
        "khPageDim[Trefoil, 3, 1, 2, 7, 0] = 0",
        "khPageDim[Trefoil, 3, 1, 3, 7, 0] = 0",
        "KhPage(A, 2, 2, -1, -1, 0){g0: 1}",
        "KhPage(B, 2, 1, -1, -1, -2){g0: 1}",
        "KhPage(BPlus, 2, 2, -1, -1, -1){g0: 1}",
    ] {
        assert!(
            output.lines().any(|actual| actual == line),
            "{line}\n{output}"
        );
    }
}

#[test]
fn page_zero_odd_signs_and_finite_arithmetic_use_actual_coefficients() {
    let source = r"
        khbraid P[2,[1]];
        khbraid Mixed[3,[2,1,-2]];
        let raw = khMod[3,0](khTerm[P,[],[0],-1]);
        let x = khPageClass[0,0](raw);
        let y = khPageD(x);
        assert khPageClass[0,0](khPageRep(x)) == x;
        assert khPageRep(y) == khMod[3,0](khTerm[P,[0],[0],2]);
        assert khD(khPageRep(khPageBoundary(y))) == khPageRep(y);
        assert khPageNext(y) == zero[KhPage(P,3,1,1,1,0)];
        assert add(raw,khScale[2](raw)) == zero[KhPC(P,3,0,1,0)];
        assert add(x,neg(x)) == zero[KhPage(P,3,0,0,1,0)];
        let f = \v : KhPage(P,3,0,0,1,0) => S.hi;
        card f;
        let c = \v : KhPC(P,3,0,1,0) => S.hi;
        card c;
        let mixed = khMod[3,1](khTerm[Mixed,[0],[1,2]]);
        assert khMod[3,-1](khD(mixed)) == khMod[3,-1](khTerm[Mixed,[0,1],[0,1],-1]);
        assert khD(khD(mixed)) == zero[KhPC(Mixed,3,2,0,1)];
        khPageDim[P,5,0,0,1,0];
        khPageDim[P,5,1,0,1,0];
        khPageDim[P,5,99,0,1,0];
        khPageDim[P,5,99,-9223372036854775808,1,0];
        khPageDim[P,5,99,0,1,-9223372036854775808];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 9 / 9\n",
            "card = 9 / 9\n",
            "khPageDim[P, 5, 0, 0, 1, 0] = 2\n",
            "khPageDim[P, 5, 1, 0, 1, 0] = 1\n",
            "khPageDim[P, 5, 99, 0, 1, 0] = 1\n",
            "khPageDim[P, 5, 99, -9223372036854775808, 1, 0] = 0\n",
            "khPageDim[P, 5, 99, 0, 1, -9223372036854775808] = 0\n",
        )
    );
}

#[test]
fn lifts_validate_lower_residuals_and_pages_validate_witnesses() {
    for (query, expected) in [
        (
            "show khPageClass[2,2](khMod[3,2](khTerm[P,[],[0]]));",
            "differential above its required cutoff",
        ),
        (
            "show khPageClass[0,0](khMod[3,2](khTerm[P,[],[]]));",
            "term above its filtration",
        ),
        (
            "show khPageBoundary(khPageClass[1,-4](khMod[2,-4](khPsi[A])));",
            "is not a d_1 boundary",
        ),
        (
            "show khPageNext(khPageGen[P,3,1,0,3,2,0]);",
            "has nonzero d_1",
        ),
        (
            "show khPagePrev(khPageGen[P,3,0,0,3,2,0]);",
            "page index is out of range",
        ),
        (
            "show khPageD(zero[KhPage(P,3,9223372036854775807,0,1,0)]);",
            "page filtration is out of range",
        ),
        (
            "show khPageGen[P,3,1,0,3,2,1];",
            "generator index is out of range",
        ),
        (
            "show khMod[2,0](khMod[3,0](khTerm[P,[],[0]]));",
            "requires the same characteristic",
        ),
        ("khPageDim[P,9,1,0,3,2];", "9 is not prime"),
    ] {
        let source = format!("khbraid P[2,[1]]; khbraid A[4,[3,-2,-2,3,3,2,-3,-1,2,1,1]]; {query}");
        let error = diagnostic(interpret(&source));
        assert!(error.contains(expected), "{query}\n{error}");
    }
}
