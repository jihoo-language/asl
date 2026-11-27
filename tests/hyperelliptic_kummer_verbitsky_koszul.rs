mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn genuine_betti_numbers_detect_the_two_cancelling_numerator_terms() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-verbitsky-koszul")),
        concat!(
            "rvVerbitskyTorTable[B, 2] = betti=[(0, 0, 1), (1, 3, 4), (1, 4, 1), (2, 4, 1), (2, 5, 4), (3, 8, 1)] numerator=[(0, 1), (3, -4), (5, 4), (8, -1)]\n",
            "rvVerbitskyTorBetti[B, 2, 1, 4] = 1\n",
            "rvVerbitskyTorBetti[B, 2, 2, 4] = 1\n",
            "rvVerbitskyMinimalRelation[3] = {3*D[0]*D[1]^2, 1*D[0]^3}\n",
        )
    );
    assert_eq!(
        stdout(interpret("rvhodgewick E[11]{-1,-1,1}{[2]}; rvVerbitskyTorTable[E,2];")),
        "rvVerbitskyTorTable[E, 2] = betti=[(0, 0, 1), (1, 3, 2), (1, 4, 3), (2, 5, 3), (2, 6, 2), (3, 9, 1)] numerator=[(0, 1), (3, -2), (4, -3), (5, 3), (6, 2), (9, -1)]\n",
    );
    assert_eq!(
        stdout(interpret("rvhodgewick B[21]{0,1}{[2,0],[0,-2]}; rvVerbitskyTorTable[B,3];")),
        "rvVerbitskyTorTable[B, 3] = betti=[(0, 0, 1), (1, 4, 5), (2, 6, 5), (3, 10, 1)] numerator=[(0, 1), (4, -5), (6, 5), (10, -1)]\n",
    );
}

#[test]
fn exterior_differential_and_cycle_projection_kill_actual_boundaries() {
    let result = stdout(interpret(
        r"
        rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
        let u=rvVerbitskyBasis[B,2,0,0];
        let x=rvVerbitskyKoszulTerm{0}(u);
        let y=rvVerbitskyKoszulTerm{1}(u);
        let dx=rvVerbitskyKoszulDiff(x);
        let dy=rvVerbitskyKoszulDiff(y);
        let xy=rvVerbitskyKoszulCup(x,y);
        assert xy==neg(rvVerbitskyKoszulCup(y,x));
        assert rvVerbitskyKoszulCup(x,x)==zero[RVVerbitskyKoszul(B,2,2,2)];
        let dxy=rvVerbitskyKoszulDiff(xy);
        assert dxy==add(rvVerbitskyKoszulCup(dx,y),neg(rvVerbitskyKoszulCup(x,dy)));
        assert rvVerbitskyKoszulDiff(dxy)==zero[RVVerbitskyKoszul(B,2,0,2)];
        assert rvVerbitskyTorClass(dxy)==zero[RVVerbitskyTor(B,2,1,2)];
        let a=rvVerbitskyTorBasis[B,2,1,3,0];
        let c=rvVerbitskyTorRepresentative(a);
        assert rvVerbitskyKoszulDiff(c)==zero[RVVerbitskyKoszul(B,2,0,3)];
        assert rvVerbitskyTorClass(c)==a;
        assert add(a,neg(a))==zero[RVVerbitskyTor(B,2,1,3)];
        assert rvVerbitskyTorCup(rvVerbitskyTorBasis[B,2,0,0,0],a)==a;
        rvVerbitskyTorBetti[B,2,-1,0];
        ",
    ));
    assert_eq!(result, "rvVerbitskyTorBetti[B, 2, -1, 0] = 0\n");
}

#[test]
fn every_homological_one_basis_element_lifts_a_concrete_minimal_relation() {
    let report = stdout(interpret(
        r"
        rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
        rvVerbitskyMinimalRelation(rvVerbitskyTorBasis[B,2,1,3,0]);
        rvVerbitskyMinimalRelation(rvVerbitskyTorBasis[B,2,1,3,1]);
        rvVerbitskyMinimalRelation(rvVerbitskyTorBasis[B,2,1,3,2]);
        rvVerbitskyMinimalRelation(rvVerbitskyTorBasis[B,2,1,3,3]);
        rvVerbitskyMinimalRelation(rvVerbitskyTorBasis[B,2,1,4,0]);
        rvhodgewick S[22]{0,1}{};
        rvVerbitskyMinimalRelation(rvVerbitskyTorBasis[S,1,1,4,0]);
        ",
    ));
    // With x=D[0], y=D[1], z=Q[0], the weight-three lifts
    // span the four cubic generators from the diagonal Fujiki
    // calculation. The weight-four lift is z²−(161/4)x⁴ modulo
    // products of the cubic relations, after scaling by 23/2.
    assert_eq!(
        report,
        concat!(
            "rvVerbitskyMinimalRelation[3] = {-2/7*D[0]*Q[0], 1*D[0]^3}\n",
            "rvVerbitskyMinimalRelation[3] = {1*D[0]*D[1]^2, 2/21*D[0]*Q[0]}\n",
            "rvVerbitskyMinimalRelation[3] = {1*D[0]^2*D[1], -2/21*D[1]*Q[0]}\n",
            "rvVerbitskyMinimalRelation[3] = {2/7*D[1]*Q[0], 1*D[1]^3}\n",
            "rvVerbitskyMinimalRelation[4] = {1*D[1]^2*Q[0], 2/23*Q[0]^2}\n",
            "rvVerbitskyMinimalRelation[4] = {1*Q[0]^2}\n",
        )
    );
}

#[test]
fn tor_products_pair_nontrivially_with_the_gorenstein_top_class() {
    let program = r"
        rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
        let x=rvVerbitskyBasis[B,2,2,0];
        let y=rvVerbitskyBasis[B,2,2,1];
        let z=rvVerbitskyClass(rvHodgeQ[B,2]{1});
        let x2=rvVerbitskyCup(x,x);
        let y2=rvVerbitskyCup(y,y);
        let x3=rvVerbitskyCup(x,x2);
        let a=rvVerbitskyTorClass(add(rvVerbitskyKoszulTerm{2}(rvVerbitskyScale[4,1](z)),rvVerbitskyKoszulTerm{0}(rvVerbitskyScale[-161,1](x3))));
        let b=rvVerbitskyTorClass(rvVerbitskyKoszulTerm{0,1}(add(rvVerbitskyScale[8,1](z),add(rvVerbitskyScale[-21,1](x2),rvVerbitskyScale[21,1](y2)))));
        let top=rvVerbitskyTorScale[-2100,1](rvVerbitskyTorBasis[B,2,3,8,0]);
        assert rvVerbitskyTorCup(a,b)==top;
        assert rvVerbitskyTorCup(b,a)==top;
        assert rvVerbitskyTorCup(a,a)==zero[RVVerbitskyTor(B,2,2,8)];
        ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn validates_cycles_marks_indices_and_gradings() {
    let marking = "rvhodgewick B[21]{0,1}{[2,0],[0,-2]}; ";
    for (query, message) in [
        ("let bad=rvVerbitskyTorClass(rvVerbitskyKoszulTerm{0}(rvVerbitskyBasis[B,2,0,0]));", "requires a Koszul cycle"),
        ("let bad=rvVerbitskyKoszulTerm{1,0}(rvVerbitskyBasis[B,2,0,0]);", "strictly increasing"),
        ("let bad=rvVerbitskyKoszulTerm{3}(rvVerbitskyBasis[B,2,0,0]);", "in range"),
        ("let bad=rvVerbitskyTorBasis[B,2,1,4,1];", "outside the Betti number"),
        ("rvVerbitskyMinimalRelation(rvVerbitskyTorBasis[B,2,3,8,0]);", "homological degree one"),
        ("rvVerbitskyTorTable[B,1000000000000000];", "storage limit"),
        ("rvVerbitskyTorTable[B,170141183460469231731687303715884105727];", "storage limit"),
        ("rvVerbitskyTorBetti[A,2,1,1];", "unknown Hodge Wick marking"),
        ("let bad=rvVerbitskyTorCup(rvVerbitskyTorBasis[B,2,0,0,0],rvVerbitskyTorBasis[B,3,0,0,0]);", "same marking and Hilbert length"),
    ] {
        assert!(diagnostic(interpret(&format!("{marking}{query}"))).contains(message), "{query}");
    }
}
