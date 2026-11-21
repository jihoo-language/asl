mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn witt_cartier_example_gives_the_pole_basis_and_rank_ladder() {
    assert_eq!(
        stdout(example("rank-witt-cartier")),
        concat!(
            "rvWittGaps[Elliptic] = [1]\n",
            "rvWittBasis[GenusTwo] = [(0, 0, 0), (0, 1, 0)]\n",
            "rvWittGaps[GenusSix] = [1, 2, 3, 5, 7, 11]\n",
            "rvWittBasis[GenusSix] = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 0, 0), (0, 0, 1), (1, 1, 0)]\n",
            "rvWittCartierRank[GenusSix, 0] = 6\n",
            "rvWittCartierRank[GenusSix, 1] = 4\n",
            "rvWittCartierRank[GenusSix, 2] = 3\n",
            "rvWittCartierRank[GenusSix, 3] = 2\n",
            "rvWittCartierRank[GenusSix, 4] = 1\n",
            "rvWittCartierRank[GenusSix, 5] = 0\n",
            "rvWittaNumber[GenusSix] = 2\n",
            "WittDiff(GenusSix, 11){6*omega5}\n",
        )
    );
}

#[test]
fn the_second_equation_cancels_two_even_poles_before_cartier() {
    let program = r"
        let f = add(rvmon[[0,1]]{3,0},add(rvmon[[0,1]]{2,1},rvmon[[0,1]]{0,3}));
        let g = add(rvmon[[0,1]]{9,0},add(rvmon[[0,1]]{7,2},rvmon[[0,1]]{0,9}));
        rvwittcover W[[0]] { f,g };
        rvWittBasis[W];
        rvWittGaps[W];
        let out = rvWittCartier(rvWittOmega[W,3,7]);
        assert out == add(rvWittOmega[W,3,8],add(rvWittOmega[W,3,5],
                          add(rvWittOmega[W,3,3],rvWittOmega[W,3,1])));
        show out;
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvWittBasis[W] = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 0, 0), (1, 1, 0), (3, 0, 0), (2, 1, 0), (0, 0, 1), (4, 0, 0)]\n",
            "rvWittGaps[W] = [1, 2, 3, 5, 7, 9, 11, 13, 17]\n",
            "WittDiff(W, 3){omega1 + omega3 + omega5 + omega8}\n",
        ),
    );
}

#[test]
fn extension_field_witt_differentials_are_first_class_and_nominal() {
    let program = r"
        rvwittcover W[[0]] {rvmon[[0]]{1},rvmon[[0]]{3}};
        rvascover A[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
        let a = rvWittOmega[W,11,0];
        let b = rvWittOmega[W,11,1];
        let apply = \v : WittDiff(W,11) => rvWittCartier(v);
        let x = add(a,rvWittDiffScale[2](b));
        let y = rvWittDiffScale[3](b);
        assert apply(x) == rvWittDiffScale[6](a);
        assert apply(y) == rvWittDiffScale[7](a);
        assert apply(add(x,y)) == add(apply(x),apply(y));
        assert apply(add(x,y)) == a;
        assert apply(apply(x)) == zero[WittDiff(W,11)];
        assert rvWittDiffScale[0](x) == zero[WittDiff(W,11)];
        assert add(x,x) == zero[WittDiff(W,11)];
        assert rvASCartier(rvASOmega[A,11,1]) == rvASOmega[A,11,0];
        card (\v : WittDiff(W,11) => S.hi);
    ";
    assert_eq!(stdout(interpret(program)), "card = 64 / 64\n");

    for (statement, reason) in [
        (
            "show rvWittCartier(rvASOmega[A,11,1]);",
            "requires a WittDiff",
        ),
        (
            "show rvASCartier(rvWittOmega[W,11,1]);",
            "requires an ASDiff",
        ),
        (
            "show rvASDiffScale[1](rvWittOmega[W,11,1]);",
            "requires an ASDiff",
        ),
        (
            "show rvWittDiffScale[8](rvWittOmega[W,11,1]);",
            "out of range for F8",
        ),
        ("show rvWittOmega[W,11,2];", "out of range"),
        ("show zero[WittDiff(Missing,11)];", "unknown Witt cover"),
    ] {
        let result = diagnostic(interpret(&format!("{program}\n{statement}")));
        assert!(result.contains(reason), "{statement}\n{result}");
    }
}

#[test]
fn isolated_cartier_columns_do_not_enumerate_the_ambient_basis() {
    let setup = "
        rvwittcover W[[0]] {rvmon[[0]]{200001},zero[RVPoly([0],0)]};
        rvWittGenus[W];
        show rvWittCartier(rvWittOmega[W,3,1]);
        rvWittCartierRank[W,1000000];
    ";
    assert_eq!(
        stdout(interpret(setup)),
        concat!(
            "rvWittGenus[W] = 500001\n",
            "WittDiff(W, 3){omega0}\n",
            "rvWittCartierRank[W, 1000000] = 0\n",
        )
    );
    for query in [
        "rvWittGaps[W];",
        "rvWittBasis[W];",
        "rvWittCartierRank[W,1];",
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{query}")));
        assert!(
            error.contains("exceeds the storage limit"),
            "{query}\n{error}"
        );
    }
}
