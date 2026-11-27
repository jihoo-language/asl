mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn residues_distinguish_cross_ratios_and_detect_a_split_quadratic_wall() {
    assert_eq!(
        stdout(example("rv-residue-specialization")),
        concat!(
            "RVResSpace(7, 4){[1, 0, 1, 1], [0, 1, 1, 2]}\n",
            "RVResSpace(7, 4){[1, 0, 1, 1], [0, 1, 1, 3]}\n",
            "rvResSupport = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]\n",
            "rvResSupport = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]\n",
            "rvResLeq = false\n",
            "rvResContains = true\n",
            "rvResLeq = true\n",
            "RVResSpace(5, 4){[1, 4, 0, 4], [0, 0, 1, 1]}\n",
            "rvResSupport = [[0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]\n",
            "RVResSpace(5, 4){[0, 0, 1, 0], [0, 0, 0, 1]}\n",
            "rvResSupport = [[2, 3]]\n",
        )
    );
}

#[test]
fn dependent_rows_basis_changes_duals_and_wide_weight_arithmetic() {
    let input = r"
        let a = rvResSpecialize[3,3,[0,-1,0]]{[1,2,3],[0,1,1],[1,3,4]};
        let b = rvResSpecialize[3,3,[9223372036854775807,9223372036854775806,9223372036854775807]]{
            [0,-1,-1],[1,(7,3),(10,3)]
        };
        assert a == b;
        assert rvResOrthogonal(a) == rvResSpecialize[3,3,[0,1,0]]{[1,1,-1]};
        assert rvResOrthogonal(rvResOrthogonal(a)) == a;

        // -v_p(3) is an actual pole in the second coordinate.
        let line = rvResSpecialize[3,2,[0,-1]]{[1,(2,3)]};
        assert line == rvResSpace[3,2]{[1,2]};
        rvResContains[[2,1]](line);
        rvResContains[[1,0]](line);

        let e = rvResSpecialize[3,2,[-9223372036854775807,9223372036854775807]]{[(1,3),1]};
        assert e == rvResSpace[3,2]{[0,1]};
        rvResDim(e);
    ";
    assert_eq!(
        stdout(interpret(input)),
        "rvResContains = true\nrvResContains = false\nrvResDim = 1\n"
    );
}

#[test]
fn sum_intersection_and_bounded_finite_enumeration() {
    let input = r"
        let a = rvResSpace[2,3]{[1,1,0],[0,1,1]};
        let b = rvResSpace[2,3]{[1,0,1],[1,0,0]};
        let c = rvResIntersection(a,b);
        assert c == rvResSpace[2,3]{[1,0,1]};
        assert rvResOrthogonal(rvResSum(a,b)) == zero[RVResSpace(2,3)];
        assert rvResOrthogonal(c) == rvResSum(rvResOrthogonal(a),rvResOrthogonal(b));
        assert rvResSpecialize[2,0,[]]{[]} == zero[RVResSpace(2,0)];
        assert rvResSpecialize[2,2,[0,0]]{[2,0],[0,4]} == rvResSpace[2,2]{[1,0],[0,1]};
        rvResSupport(zero[RVResSpace(7,0)]);
        card(\a : RVResSpace(2,3) => S.hi);
        card(\a : RVResSpace(7,2) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(input)),
        "rvResSupport = [[]]\ncard = 16 / 16\ncard = 10 / 10\n"
    );
}

#[test]
fn specialization_and_subspace_diagnostics() {
    for (source, reason) in [
        ("show rvResSpace[7,3]{[0,1]};", "rows of length 3"),
        (
            "show rvResSpecialize[7,1,[0]]{[(0,0)]};",
            "denominator must be positive",
        ),
        ("show rvResSpecialize[7,2,[0]]{[0,1]};", "needs 2 weights"),
        ("show rvResSpace[3,65]{};", "ground-set limit of 64"),
        ("show rvResSpace[9,1]{};", "9 is not prime"),
        (
            "show rvResIntersection(rvResSpace[2,2]{},rvResSpace[2,3]{});",
            "RVResSpace(2, 2)",
        ),
        (
            "rvResLeq(rvResSpace[2,2]{},rvResSpace[5,2]{});",
            "RVResSpace(2, 2)",
        ),
        (
            "rvResContains[[1]](rvResSpace[2,2]{});",
            "vector must have 2 coordinates",
        ),
        ("rvResDim(rvValMat[3,1]{[1]});", "an RVResSpace"),
        (
            "card(\\a : RVResSpace(7,7) => S.hi);",
            "exceeds the enumeration limit",
        ),
    ] {
        let output = diagnostic(interpret(source));
        assert!(output.contains(reason), "{source}\n{output}");
    }
}

#[test]
fn the_chosen_quadratic_flag_and_its_orthogonal_specialize_at_opposite_weights() {
    let input = r"
        rvqpadic K[5,-1,2];
        let z = zero[RVQpQuad(K)];
        let one = rvQpQuad[K,1,0,1];
        let beta = rvQpQuad[K,-1,-2,1];
        let u = add(one,rvQpQuadPow[19](beta));
        rvphodgequad Flag[K] {
            {one,rvQpQuad[K,2,0,1],rvQpQuad[K,3,0,1],rvQpQuad[K,4,0,1]}, {},
            {(3,[one,z,one,one]),(2,[z,one,one,u]),(0,[z,z,one,z]),(0,[z,z,z,one])}
        };
        let filtered = rvResFil[Flag,2,[0,0,19,19]];
        assert filtered == rvResSpecializeQuad[K,4,[0,0,19,19]]{[one,z,one,one],[z,one,one,u]};
        let dual = rvResSpecializeQuad[K,4,[0,0,-19,-19]]{[neg(one),neg(one),one,z],[neg(one),neg(u),z,one]};
        assert rvResOrthogonal(filtered) == dual;
        assert rvResFil[Flag,4,[0,0,19,19]] == zero[RVResSpace(5,4)];
        rvResLeq(rvResFil[Flag,3,[0,0,19,19]],filtered);
        show filtered;

        rvqpadic J[5,-1,3];
        let jz = zero[RVQpQuad(J)];
        let jo = rvQpQuad[J,1,0,1];
        let ju = add(jo,rvQpQuadPow[19](rvQpQuad[J,-1,-2,1]));
        show rvResSpecializeQuad[J,4,[0,0,19,19]]{[jo,jz,jo,jo],[jz,jo,jo,ju]};
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "rvResLeq = true\n",
            "RVResSpace(5, 4){[1, 4, 0, 2], [0, 0, 1, 1]}\n",
            "RVResSpace(5, 4){[0, 0, 1, 0], [0, 0, 0, 1]}\n",
        )
    );
}
