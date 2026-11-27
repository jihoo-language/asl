mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn actual_columns_and_page_one_connes_map_agree_in_odd_and_even_weights() {
    assert_eq!(
        stdout(example("rank-hoch-mixed")),
        concat!(
            "MixedE(Odd, 1, -3, 4){g0}\n",
            "MixedE(Odd, 1, -4, 5){g0}\n",
            "MixedE(Even, 2, 3, 0){g0}\n",
        )
    );
}

#[test]
fn nonformal_flag_fiber_preserves_higher_products_and_allows_no_later_differentials() {
    let source = format!(
        r"{}
        rvHHMixed F[Flag,[1,1,1,1,1,1]];
        let original = rvHHToMixed[F](w);
        assert rvMixedToHH(original) == w;
        assert rvMixedToHH(mixb(original)) == add(left,ends);
        let cyclic = rvHHWord(a,b,c);
        assert mixB(rvHHToMixed[F](cyclic)) == rvHHToMixed[F](rvConnesB(cyclic));

        let C4 = rvHHToMixed[F](coefficient);
        assert C4 == mixgen[F,1,4,0,0];
        assert rvMixedToHH(C4) == coefficient;
        assert rvMixedToHH(mixD(C4)) == rvConnesB(coefficient);

        mixdim[F,1,4,0];
        mixdim[F,1,5,0];
        mixdim[F,1,6,0];
        mixdim[F,2,4,0];
        mixdim[F,2,5,0];
        mixdim[F,2,6,0];
        mixdim[F,2,2,1];
        mixdim[F,2,3,1];
        mixdim[F,2,4,1];
        ",
        include_str!("../examples/rank-hochhom.asl"),
    );
    assert!(stdout(interpret(&source)).ends_with(concat!(
        "mixdim[F, 1, 4, 0] = 1\n",
        "mixdim[F, 1, 5, 0] = 3\n",
        "mixdim[F, 1, 6, 0] = 2\n",
        "mixdim[F, 2, 4, 0] = 0\n",
        "mixdim[F, 2, 5, 0] = 1\n",
        "mixdim[F, 2, 6, 0] = 2\n",
        "mixdim[F, 2, 2, 1] = 0\n",
        "mixdim[F, 2, 3, 1] = 0\n",
        "mixdim[F, 2, 4, 1] = 0\n",
    )));
}

#[test]
fn imported_values_interoperate_with_finite_field_literals_and_predicates() {
    let source = r"
        rvideal J[[0,1,2]] { rvmon[[0,1,2]]{1,1,0}, rvmon[[0,1,2]]{1,0,1}, rvmon[[0,1,2]]{0,1,1} };
        rvHHMixed C[J,[1,1,1]];
        let a = rvHHgen[J,2,[1,1,1],0];
        let b = rvHHgen[J,2,[1,1,1],1];
        let M = rvHHToMixed[C](add(a,b));
        assert M == add(rvHHToMixed[C](a),rvHHToMixed[C](b));
        assert rvMixedToHH(M) == add(a,b);
        assert rvMixedToHH(mixD(M)) == rvConnesB(add(a,b));
        card \x : MixedChain(C,2) => S.hi;
        card \x : MixedE(C,1,2,0) => S.hi;
        card \x : MixedChain(C,-1) => S.hi;
        ";
    assert_eq!(
        stdout(interpret(source)),
        "card = 4 / 4\ncard = 4 / 4\ncard = 1 / 1\n"
    );
}

#[test]
fn the_bridge_checks_origin_weight_page_and_materialization_limits() {
    let source = include_str!("../examples/rank-hoch-mixed.asl");
    for (statement, message) in [
        (
            "show rvHHToMixed[Even](oddLow);",
            "requires the ideal and weight",
        ),
        (
            "show rvHHToMixed[Odd,0](oddLow);",
            "takes only a mixed complex name",
        ),
        (
            "show rvMixedToHH(mixgen[Odd,0,5,0,0]);",
            "page-1 MixedE class",
        ),
        ("show rvMixedToHH(mixU[0](raw));", "an imported MixedChain"),
        (
            "mixedcomplex M[2] {x:0} {} {}; show rvMixedToHH(mixchain[M,0][(x,1)]);",
            "imported by rvHHMixed",
        ),
        ("rvHHMixed Odd[Exterior,[3]];", "already defined"),
        ("rvHHMixed TooLarge[Exterior,[100000]];", "storage limit"),
    ] {
        let error = diagnostic(interpret(&format!("{source} {statement}")));
        assert!(error.contains(message), "{statement}: {error}");
    }
}
