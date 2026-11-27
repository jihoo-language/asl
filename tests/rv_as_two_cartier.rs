mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn stable_ordinary_cartier_distinguishes_constant_twists() {
    assert_eq!(
        stdout(example("rank-as-two-cartier")),
        concat!(
            "rvASTwoBasis[K] = [(-1, 0), (-1, 1), (-1, 2), (-1, 3)]\n",
            "ASTwoDiff(K){2*omega(-1,0) + omega(-1,2)}\n",
            "ASTwoDiff(K){omega(-1,1) + omega(-1,3)}\n",
            "rvASTwoNewton[K] = [0: 4, 1: 4]\n",
            "rvASTwoNewton[Twist] = [0: 4, 1: 4]\n",
            "rvASTwoOrdinaryRanks[K] = T-1 [4, 3, 2, 1, 0], C-1 [4, 2, 0]\n",
            "rvASTwoOrdinaryRanks[Twist] = T-1 [4, 3, 2, 1, 0], C-1 [4, 3, 2, 1, 0]\n",
            "rvASTwoHasseWitt[K] = [1, 1, 1, 1, 1]\n",
            "rvASTwoHasseWitt[Twist] = [1, 1, 1, 1, 1]\n",
            "rvASTwoCartierRank[Wild, 0] = 12\n",
            "rvASTwoCartierRank[Wild, 1] = 6\n",
            "rvASTwoCartierRank[Wild, 2] = 5\n",
            "rvASTwoStableRank[Wild] = 4\n",
            "rvASTwoaNumber[Wild] = 6\n",
            "rvASTwoOrdinaryRanks[Wild] = T-1 [4, 3, 2, 1, 0], C-1 [4, 1, 0]\n",
            "rvASTwoHasseWitt[Wild] = [1, 1, 1, 1, 1]\n",
            "ASTwoDiff(Wild){omega(-1,2)}\n",
        )
    );
}

#[test]
fn cartier_commutes_with_deck_action_across_the_signed_basis() {
    let mut source = String::from("rvastwo W[5] {(4,1),(-2,1)}; rvastwo Inv[5] {(-4,1),(2,1)};\n");
    for (j, (low, high)) in [(-2, 2), (-2, 1), (-1, 0), (-1, -1)]
        .into_iter()
        .enumerate()
    {
        for i in low..=high {
            let name = format!("x{j}a{}", i + 2);
            writeln!(&mut source, "let {name} = rvASTwoOmega[W,{i},{j}];").unwrap();
            writeln!(&mut source, "assert rvASTwoCartier(rvASTwoTranslate[1]({name})) == rvASTwoTranslate[1](rvASTwoCartier({name}));").unwrap();
            writeln!(
                &mut source,
                "assert rvASTwoTranslate[-1](rvASTwoTranslate[6]({name})) == {name};"
            )
            .unwrap();
            writeln!(
                &mut source,
                "assert rvASTwoCartierIter[15]({name}) == rvASTwoCartierIter[1000000000]({name});"
            )
            .unwrap();
            writeln!(
                &mut source,
                "assert rvASTwoCartierIter[0]({name}) == {name};"
            )
            .unwrap();
        }
    }
    source.push_str(
        r"
        let x = rvASTwoOmega[W,0,2];
        assert rvASTwoCartier(x) == pscale[3](rvASTwoOmega[W,0,1]);
        // Pullback under inversion sends the source to -x^-2 y^2 dx.
        let pullback = pscale[4](rvASTwoOmega[Inv,-2,2]);
        assert rvASTwoCartier(pullback) == pscale[2](rvASTwoOmega[Inv,-2,1]);
        assert rvASTwoCartier(pscale[4](x)) == pscale[4](rvASTwoCartier(x));
        assert rvASTwoCartier(zero[ASTwoDiff(W)]) == zero[ASTwoDiff(W)];
        rvASTwoCartierRank[W,1000000000];
        ",
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvASTwoCartierRank[W, 1000000000] = 4\n"
    );
}

#[test]
fn regularity_ordinary_basis_and_nominal_types_are_checked() {
    let setup = "rvastwo K[3] {(1,1),(-1,1)}; rvastwo B[3] {(1,1),(-1,1)}; rvasodd O[3] {(2,1)};";
    assert_eq!(
        stdout(interpret(&format!(
            r"{setup} card (\x : ASTwoDiff(K) => S.hi); rvASTwoCartierRank[K,1];"
        ))),
        "card = 9 / 9\nrvASTwoCartierRank[K, 1] = 2\n",
    );
    for (statement, reason) in [
        (
            "show rvASTwoOmega[K,-2,0];",
            "outside the two-pole regular differential basis",
        ),
        (
            "show rvASTwoOmega[K,-1,2];",
            "outside the two-pole regular differential basis",
        ),
        ("show rvASTwoOrdinaryGen[K,2];", "outside dimension 2"),
        (
            "show add(zero[ASTwoDiff(K)],zero[ASTwoDiff(B)]);",
            "expected type ASTwoDiff(K)",
        ),
        (
            "show rvASTwoCartier(rvASOddOmega[O,0,0]);",
            "an ASTwoDiff differential",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
