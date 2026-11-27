mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn long_torus_knot_levels_and_exact_rational_lifts() {
    assert_eq!(
        stdout(example("lee-rasmussen")),
        concat!(
            "leeS[Trefoil] = 2\n",
            "leeS[Mirror] = -2\n",
            "leeS[FigureEight] = 0\n",
            "leeLevels[Long, 0] = [13, 15]\n",
            "leeRank[Long, 0, 16] = 0\n",
            "leeRank[Long, 0, 14] = 1\n",
            "leeRank[Long, 0, 13] = 2\n",
            "leeS[Long] = 14\n",
            "leeFiltration = 15\n",
            "leeFiltration = 13\n",
            "leeIndependent = true\n",
            "leeFiltration = 1\n",
            "leeFiltration = boundary\n",
        )
    );
}

#[test]
fn genuine_lee_saddles_and_grade_four_chain_perturbation() {
    let program = r"
        khbraid Pos[2,[1]];
        let ax = khTerm[Pos,[],[0,1]];
        let m = leePhi(ax);
        assert m == khTerm[Pos,[0],[]];
        assert khD(ax) == zero[KhC(Pos,1,-1)];
        assert leeD(leeLift(ax)) == leeFilter[-1](leeLift(m));
        assert leePhi(m) == zero[KhC(Pos,2,7)];
        assert leeD(leeD(leeLift(ax))) == zero[LeeC(Pos,2,-1)];
        show leeScale[-2,6](leeD(leeLift(ax)));

        khbraid Neg[2,[-1]];
        let x = khTerm[Neg,[],[0]];
        let a = khD(x);
        let b = leePhi(x);
        assert a == khTerm[Neg,[0],[0,1]];
        assert b == khTerm[Neg,[0],[]];
        assert leeD(leeLift(x)) == add(leeLift(a),leeFilter[-3](leeLift(b)));
        let z = leeRep[Neg,0,1,0];
        let w = leeRep[Neg,0,-1,0];
        leeIndependent(z,w);
        leeIndependent(z,leeFilter[-1](leeScale[1,2](z)));
        leeIndependent();
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "LeeC(Pos, 1, -1){([0], []): -1/3}\n",
            "leeIndependent = true\n",
            "leeIndependent = false\n",
            "leeIndependent = true\n",
        )
    );
}

#[test]
fn both_markov_stabilizations_and_the_long_mirror() {
    let program = r"
        khbraid Plus[3,[1,1,1,2]];
        khbraid Minus[3,[1,1,1,-2]];
        khbraid Mirror[3,[-1,-2,-1,-2,-1,-2,-1,-2,-1,-2,-1,-2,-1,-2,-1,-2]];
        leeS[Plus];
        leeS[Minus];
        leeLevels[Mirror,0];
        leeS[Mirror];
        leeRank[Mirror,0,-14];
        leeRank[Mirror,0,-15];
        let g = leeRep[Mirror,0,-15,0];
        let h = leeRep[Mirror,0,-13,0];
        leeFiltration(g);
        leeFiltration(h);
        leeIndependent(g,h);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "leeS[Plus] = 2\n",
            "leeS[Minus] = 2\n",
            "leeLevels[Mirror, 0] = [-15, -13]\n",
            "leeS[Mirror] = -14\n",
            "leeRank[Mirror, 0, -14] = 1\n",
            "leeRank[Mirror, 0, -15] = 2\n",
            "leeFiltration = -15\n",
            "leeFiltration = -13\n",
            "leeIndependent = true\n",
        )
    );
}

#[test]
fn full_link_filtration_has_multiple_homological_degrees_and_multiplicities() {
    let program = r"
        khbraid Hopf[2,[1,1]];
        khbraid Unlink[2,[]];
        leeLevels[Hopf,0];
        leeLevels[Hopf,1];
        leeLevels[Hopf,2];
        leeLevels[Unlink,0];
        leeRank[Unlink,0,0];
        let g = leeRep[Unlink,0,0,0];
        let h = leeRep[Unlink,0,0,1];
        leeIndependent(g,h);
        assert leeD(g) == zero[LeeC(Unlink,1,0)];
        assert leeNull(zero[LeeC(Unlink,0,2)]) == zero[LeeC(Unlink,-1,-2)];
        card(\x : LeeC(Unlink,0,12) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "leeLevels[Hopf, 0] = [0, 2]\n",
            "leeLevels[Hopf, 1] = []\n",
            "leeLevels[Hopf, 2] = [4, 6]\n",
            "leeLevels[Unlink, 0] = [-2, 0, 0, 2]\n",
            "leeRank[Unlink, 0, 0] = 3\n",
            "leeIndependent = true\n",
            "card = 1 / 1\n",
        )
    );
}

#[test]
fn lee_operations_enforce_rational_domains_cycles_and_filtrations() {
    let initial = r"
        khbraid Hopf[2,[1,1]];
        khbraid Neg[2,[-1]];
        let z = leeRep[Neg,0,1,0];
        let x = leeLift(khTerm[Neg,[],[0]]);
        let raw = leeLift(khTerm[Neg,[0],[0,1]]);
    ";
    for (input, message) in [
        ("leeS[Hopf];", "defined here for knots only"),
        ("leeFiltration(x);", "not a cycle"),
        ("show leeNull(z);", "not a boundary"),
        ("show leeFilter[1](raw);", "term below filtration 1"),
        ("show leeRep[Neg,0,0,0];", "index is out of range"),
        ("show leeScale[1,0](z);", "denominator must be positive"),
        ("show leePhi(z);", "a homogeneous integral Khovanov chain"),
        ("show leeD(khTerm[Neg,[],[0]]);", "a Lee cochain"),
        (
            "leeIndependent(z,zero[LeeC(Hopf,0,1)]);",
            "one braid and homological degree",
        ),
        (
            "leeIndependent(z,zero[LeeC(Neg,1,1)]);",
            "one braid and homological degree",
        ),
        (
            "card(\\x : LeeC(Neg,0,1) => S.hi);",
            "cannot enumerate a nonzero rational Lee cochain",
        ),
    ] {
        let output = diagnostic(interpret(&format!("{initial}\n{input}")));
        assert!(output.contains(message), "{input}\n{output}");
    }
}
