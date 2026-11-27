mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_legendre_frobenius_roots_and_two_chosen_filtrations() {
    assert_eq!(
        stdout(example("rv-quadratic-phodge")),
        concat!(
            "rvQpQuadVal = 0\n",
            "rvQpQuadVal = 1\n",
            "rvQpQuadResidue = 169738\n",
            "legendreDwork[5, 2, 8] mod 390625 = z=280182, B_8=135706, B_7=5837, ratio=169738, unit=169738, trace=-2\n",
            "rvQpQuadVal = 8\n",
            "rvQpQuadResidue = 0\n",
            "rvQpQuadResidue = 390625\n",
            "rvPHodgeWeak[Axis] = destabilized by [0] (tH 1, tN 0)\n",
            "rvPHodgeHN[Axis] = [([0], 1), ([0, 1], -1)]\n",
            "rvPHodgeWeak[Tilted] = weakly admissible\n",
            "rvPHodgeHN[Tilted] = [([0, 1], 0)]\n",
            "rvQpQuadVal = 1\n",
            "rvQpQuadVal = 0\n",
        )
    );
}

#[test]
fn exact_quadratic_cancellation_changes_the_filtration_intersection() {
    let input = r"
        rvqpadic K[5,-1,2];
        let z = zero[RVQpQuad(K)];
        let one = rvQpQuad[K,1,0,1];
        let i = rvQpQuad[K,0,1,1];
        let c = rvQpQuad[K,1,1,1];
        let p = rvQpQuad[K,5,0,1];
        let p2 = rvQpQuadPow[2](p);
        rvphodgequad D[K] {
            {c,rvQpQuadMul(p,c),rvQpQuadMul(p2,c)},
            {(0,1,i),(1,2,rvQpQuadInv(i))},
            {(2,[z,one,i]),(1,[one,i,neg(one)]),(0,[z,z,one])}
        };
        rvPHodgeData[D,[0]];
        rvPHodgeData[D,[0,1]];
        rvPHodgeHN[D];
        rvphodgequad Changed[K] {
            {c,rvQpQuadMul(p,c),rvQpQuadMul(p2,c)},
            {(0,1,i),(1,2,rvQpQuadInv(i))},
            {(2,[z,p,rvQpQuadMul(p,i)]),(1,[one,i,neg(one)]),(0,[z,z,one])}
        };
        rvPHodgeData[Changed,[0]];
        rvPHodgeHN[Changed];
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvPHodgeData[D, [0]] = tH 1, tN 0, degree 1, Hodge [1], Newton [0]\n",
            "rvPHodgeData[D, [0, 1]] = tH 1, tN 1, degree 0, Hodge [0, 1], Newton [0, 1]\n",
            "rvPHodgeHN[D] = [([0], 1), ([0, 1, 2], -1/2)]\n",
            "rvPHodgeData[Changed, [0]] = tH 1, tN 0, degree 1, Hodge [1], Newton [0]\n",
            "rvPHodgeHN[Changed] = [([0], 1), ([0, 1, 2], -1/2)]\n",
        )
    );
}

#[test]
fn residues_detect_deep_cancellation_even_with_a_nonunit_denominator() {
    let input = r"
        rvqpadic K[5,-1,2];
        let beta = rvQpQuad[K,-1,-2,1];
        let q = rvQpQuad[K,-1,-2,5];
        assert add(q,neg(q)) == zero[RVQpQuad(K)];
        assert rvQpQuadMul(rvQpQuad[K,5,0,1],q) == beta;
        assert rvQpQuad[K,0,4,2] == rvQpQuad[K,0,2,1];
        assert rvQpQuadPow[0](q) == rvQpQuad[K,1,0,1];
        let large = rvQpQuad[K,1000000000000000000000000000000,0,1];
        assert rvQpQuadMul(large,rvQpQuadInv(large)) == rvQpQuad[K,1,0,1];
        rvQpQuadVal(q);
        rvQpQuadVal(zero[RVQpQuad(K)]);
        rvQpQuadVal(rvQpQuad[K,-1,2,5]);
        rvQpQuadResidue[8](q);
        rvQpQuadResidue[40](rvQpQuad[K,0,1,1]);
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvQpQuadVal = 0\n",
            "rvQpQuadVal = infinity\n",
            "rvQpQuadVal = -1\n",
            "rvQpQuadResidue = 278552\n",
            "rvQpQuadResidue = 2224618918409236552857702057\n",
        )
    );
}

#[test]
fn embedding_and_filtered_operator_errors_are_reported() {
    let setup = r"
        rvqpadic K[5,-1,2];
        rvqpadic L[5,-1,3];
        let z = zero[RVQpQuad(K)];
        let one = rvQpQuad[K,1,0,1];
        let two = rvQpQuad[K,2,0,1];
        let other = rvQpQuad[L,1,0,1];
    ";
    for (statement, reason) in [
        ("rvqpadic B[2,-1,1];", "requires an odd prime"),
        ("rvqpadic B[5,9,2];", "radicand must be nonsquare"),
        ("rvqpadic B[5,5,0];", "requires a nonzero root"),
        ("rvqpadic B[5,-1,1];", "requires a nonzero root"),
        (
            "rvQpQuadResidue[8](rvQpQuad[K,-1,2,5]);",
            "requires an integral element",
        ),
        ("rvQpQuadResidue[0](one);", "precision must be positive"),
        (
            "rvQpQuadResidue[1000000000](one);",
            "residue modulus at the required precision is out of range",
        ),
        ("show rvQpQuad[K,1,1,0];", "denominator must be positive"),
        ("show rvQpQuadInv(z);", "zero has no inverse"),
        ("show rvQpQuadMul(one,other);", "RVQpQuad(K)"),
        ("rvphodgequad A[K] {{other},{},{(0,[one])}};", "RVQpQuad(K)"),
        (
            "rvphodgequad A[K] {{one,one},{},{(0,[one,z]),(1,[z,one])}};",
            "distinct nonzero",
        ),
        (
            "rvphodgequad A[K] {{one,two},{(0,1,one)},{(0,[one,z]),(1,[z,one])}};",
            "fails N phi = p phi N",
        ),
        (
            "rvphodgequad A[K] {{one,two},{},{(0,[one,one]),(1,[two,two])}};",
            "basis over characteristic zero",
        ),
    ] {
        let output = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(output.contains(reason), "{statement}\n{output}");
    }
}
