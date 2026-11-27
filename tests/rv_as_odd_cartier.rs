mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn genus_ten_cartier_iterates_distinguish_the_same_a_number() {
    assert_eq!(
        stdout(example("rank-as-odd-cartier")),
        concat!(
            "rvASOddBasis[Cubic] = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (0, 3)]\n",
            "rvASOddCartierRank[Pure, 1] = 0\n",
            "rvASOddaNumber[Cubic] = 6\n",
            "rvASOddaNumber[Quartic] = 6\n",
            "rvASOddCartierRank[Cubic, 0] = 10\n",
            "rvASOddCartierRank[Cubic, 1] = 4\n",
            "rvASOddCartierRank[Cubic, 2] = 1\n",
            "rvASOddCartierRank[Cubic, 3] = 0\n",
            "rvASOddCartierRank[Quartic, 0] = 10\n",
            "rvASOddCartierRank[Quartic, 1] = 4\n",
            "rvASOddCartierRank[Quartic, 2] = 2\n",
            "rvASOddCartierRank[Quartic, 3] = 1\n",
            "rvASOddCartierRank[Quartic, 4] = 0\n",
            "rvASOddCharacterNewton[Cubic, 1] = [3/8: 2, 1/2: 1, 5/8: 2]\n",
            "rvASOddCharacterNewton[Quartic, 1] = [1/4: 1, 1/2: 3, 3/4: 1]\n",
            "ASOddDiff(Quartic){2*omega(2,0) + 2*omega(0,2)}\n",
            "ASOddDiff(Quartic){omega(0,1)}\n",
            "ASOddDiff(Quartic){4*omega(0,0)}\n",
            "ASOddDiff(Quartic){omega(0,0) + 3*omega(0,1) + 3*omega(0,2) + omega(0,3)}\n",
        )
    );
}

#[test]
fn signed_cartier_and_translation_are_linear_and_commute() {
    let input = r"
        rvasodd Cubic[5] {(6,1),(3,1)};
        let a = rvASOddOmega[Cubic,0,3];
        let b = rvASOddOmega[Cubic,1,2];
        let c = rvASOddOmega[Cubic,1,0];
        let d = rvASOddOmega[Cubic,1,1];
        let e = rvASOddOmega[Cubic,0,1];
        assert rvASOddCartier(a) == add(pscale[4](c),d);
        assert rvASOddCartier(b) == pscale[3](e);
        assert rvASOddCartier(add(a,b)) == add(rvASOddCartier(a),rvASOddCartier(b));
        assert rvASOddCartier(pscale[3](a)) == pscale[3](rvASOddCartier(a));
        assert rvASOddCartierIter[0](a) == a;
        assert rvASOddCartierIter[2](a) == rvASOddCartier(rvASOddCartier(a));
        assert rvASOddCartierIter[1000](a) == zero[ASOddDiff(Cubic)];
        assert rvASOddTranslate[1](a) == rvASOddTranslate[2](rvASOddTranslate[4](a));
        assert rvASOddTranslate[-1](a) == rvASOddTranslate[4](a);
        assert rvASOddTranslate[5](a) == a;
        assert rvASOddTranslate[1](rvASOddTranslate[1](rvASOddTranslate[1](rvASOddTranslate[1](rvASOddTranslate[1](a))))) == a;
        assert rvASOddCartierIter[2](rvASOddTranslate[2](a)) == rvASOddTranslate[2](rvASOddCartierIter[2](a));
    ";
    assert_eq!(stdout(interpret(input)), "");
}

#[test]
fn empty_and_large_regular_bases_follow_their_sparse_bounds() {
    let program = r"
        rvasodd Line[7] {(7,1)};
        rvasodd Elliptic[3] {(2,1)};
        rvASOddBasis[Line];
        rvASOddCartierRank[Line,100000];
        rvASOddaNumber[Line];
        card (\x : ASOddDiff(Line) => S.hi);
        card (\x : ASOddDiff(Elliptic) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvASOddBasis[Line] = []\n",
            "rvASOddCartierRank[Line, 100000] = 0\n",
            "rvASOddaNumber[Line] = 0\n",
            "card = 1 / 1\n",
            "card = 3 / 3\n",
        )
    );
    let degree = (1usize << (usize::BITS / 2)) + 1;
    let setup = format!("rvasodd Huge[3] {{({degree},1)}};");
    assert_eq!(
        stdout(interpret(&format!(
            "{setup} rvASOddCartierRank[Huge,0]; show rvASOddCartier(rvASOddOmega[Huge,2,0]);"
        ))),
        format!(
            "rvASOddCartierRank[Huge, 0] = {}\nASOddDiff(Huge){{omega(0,0)}}\n",
            degree - 1,
        )
    );
    for query in ["rvASOddBasis[Huge];", "rvASOddCartierRank[Huge,1];"] {
        let error = diagnostic(interpret(&format!("{setup} {query}")));
        assert!(error.contains("storage limit"), "{error}");
    }
}

#[test]
fn invalid_monomials_and_mismatched_differentials_have_explicit_diagnostics() {
    let setup = "rvasodd A[5] {(6,1)}; rvasodd B[5] {(6,1)};";
    for (statement, reason) in [
        (
            "show rvASOddOmega[A,4,0];",
            "outside the regular differential basis",
        ),
        (
            "show rvASOddOmega[A,0,4];",
            "outside the regular differential basis",
        ),
        (
            "show rvASOddCartier(rvASOddCoeff[A,1,0]);",
            "an ASOddDiff differential",
        ),
        ("show rvASOddCartierIter[1](S.lo);", "rvASOddCartierIter"),
        (
            "show add(zero[ASOddDiff(A)],zero[ASOddDiff(B)]);",
            "expected type ASOddDiff(A)",
        ),
        (
            "show zero[ASOddDiff(Missing)];",
            "unknown odd Artin–Schreier cover",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
