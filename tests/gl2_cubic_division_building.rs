mod common;

use std::collections::HashSet;
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

const BASE: &str = "
    rvunram B[3,7]{0,1};
    rvcyclotomiclocal K[B,1];
";

fn setup(invariant: usize, precision: usize) -> String {
    format!(
        "{BASE}
        rvcubicdiv D[K,{invariant},{precision}];
        let one = rvCubicDivElem[D]([1]);
        let nul = zero[RVCubicDiv(D)];
        let th = rvCubicDivTheta[D];
        let pi = rvCubicDivPi[D];
        let o = rvCubicTreeRoot[D];
        let g = rvCubicTreeMatrix(one,th,pi,one);
        let sigma = rvCubicTreePi[D];"
    )
}

#[test]
fn spheres_incidence_and_radial_hecke_coefficients() {
    let output = stdout(example("gl2-cubic-division-building"));
    for expected in [
        "rvCubicTreeSphere[D,1] = size=28",
        "rvCubicTreeSphere[D,2] = size=756",
        "rvCubicTreeSphere[D,3] = size=20412",
        "rvCubicTreeLevel = 2",
        "rvCubicTreePath = radius=0, root -> radius=1,",
        "rvCubicTreeHecke[1] = distance=0, paths=28",
        "rvCubicTreeHecke[1] = distance=2, paths=1",
        "rvCubicTreeHecke[2] = distance=1, paths=27",
    ] {
        assert!(output.contains(expected), "{expected}\n{output}");
    }
    assert_eq!(
        output.matches("rvCubicTreeDistance = 2").count(),
        2,
        "{output}"
    );

    let boundary = stdout(interpret(&format!(
        "{}
        let end = rvCubicTreeAt[D,2,755];
        let parent = rvCubicTreeAncestor[1](end);
        assert rvCubicTreeChild[26](parent) == end;
        rvCubicTreeHecke[3](end,o);
        rvCubicTreeHecke[2](end,o);
        rvCubicTreeHecke[0](end,parent);",
        setup(1, 2)
    )));
    assert_eq!(
        boundary,
        concat!(
            "rvCubicTreeHecke[3] = distance=2, paths=27\n",
            "rvCubicTreeHecke[2] = distance=2, paths=0\n",
            "rvCubicTreeHecke[0] = distance=1, paths=1\n",
        )
    );
}

#[test]
fn full_shell_permutation_and_normalizer_fixed_points_over_both_invariants() {
    for t in [1, 2] {
        let mut program = setup(t, 2);
        writeln!(
            program,
            "
            let sig2 = rvCubicTreeCompose(sigma,sigma);
            let id = rvCubicTreeMatrix(one,nul,nul,one);
            assert rvCubicTreeCompose(sig2,sigma) == id;
            let gi = rvCubicTreeInverse(g);
            assert rvCubicTreeCompose(g,gi) == id;
            assert rvCubicTreeCompose(gi,g) == id;
            let j = rvCubicTreeMatrix(nul,one,one,nul);
            assert rvCubicTreeInverse(j) == j;
            assert rvCubicTreeFixes(g,rvCubicTreeAt[D,1,0]) == S.hi;
            assert rvCubicTreeFixes(g,rvCubicTreeAt[D,2,0]) == S.lo;
            let x = rvCubicTreeAt[D,2,1];
            let xcoord = rvCubicDivMul(th,th);
            assert x == rvCubicTreeLine[2](th,rvCubicDivMul(xcoord,th));
            let expected = rvCubicDivElem[D]([1],[2],[1],[{}],[{}],[1]);
            assert rvCubicTreeApply(g,x) == rvCubicTreeLine[2](one,expected);
            let gp = rvCubicTreeCompose(g,sigma);
            assert rvCubicTreeApply(gp,x) == rvCubicTreeApply(g,rvCubicTreeApply(sigma,x));
            assert rvCubicTreeApply(rvCubicTreeInverse(gp),rvCubicTreeApply(gp,x)) == x;",
            if t == 1 { 2 } else { 0 },
            if t == 1 { 0 } else { 2 }
        )
        .unwrap();
        for i in 0..28 {
            writeln!(
                program,
                "show rvCubicTreeFixes(sigma,rvCubicTreeAt[D,1,{i}]);"
            )
            .unwrap();
        }
        let small = stdout(interpret(&program));
        assert_eq!(small.lines().filter(|x| *x == "S.hi").count(), 4, "{small}");
        assert_eq!(small.lines().count(), 28, "{small}");

        let mut shell = HashSet::new();
        let mut image = HashSet::new();
        let mut fixed = 0;
        for start in (0..756).step_by(126) {
            let mut program = setup(t, 2);
            for i in start..start + 126 {
                writeln!(
                    program,
                    "
                    let x{i} = rvCubicTreeAt[D,2,{i}];
                    show x{i}; show rvCubicTreeApply(g,x{i}); show rvCubicTreeFixes(sigma,x{i});"
                )
                .unwrap();
            }
            let output = stdout(interpret(&program));
            let rows: Vec<_> = output.lines().collect();
            assert_eq!(rows.len(), 126 * 3);
            for triple in rows.chunks_exact(3) {
                assert!(shell.insert(triple[0].to_string()));
                assert!(image.insert(triple[1].to_string()));
                fixed += usize::from(triple[2] == "S.hi");
            }
        }
        assert_eq!(shell.len(), 756);
        assert!(shell == image);
        assert_eq!(fixed, 12);
    }
}

#[test]
fn translation_moves_outward_or_inward_by_the_correct_right_lattice_formulas() {
    for t in [1, 2] {
        let output = stdout(interpret(&format!(
            "{}
            let thp = rvCubicDivMul(th,pi);
            let left = rvCubicTreeLine[1](one,nul);
            let right = rvCubicTreeLine[1](nul,one);
            let inside = rvCubicTreeLine[2](one,thp);
            let across = rvCubicTreeLine[2](one,th);
            let outside = rvCubicTreeLine[2](thp,one);
            assert rvCubicTreeTranslate(o) == right;
            assert rvCubicTreeTranslateBack(o) == left;
            assert rvCubicTreeTranslate(left) == o;
            assert rvCubicTreeTranslateBack(right) == o;
            assert rvCubicTreeTranslate(inside) == rvCubicTreeLine[1](one,th);
            assert rvCubicTreeTranslate(across) ==
                rvCubicTreeLine[3](rvCubicDivMul(pi,rvCubicDivInverse(th)),one);
            assert rvCubicTreeTranslate(outside) == rvCubicTreeLine[3](rvCubicDivMul(pi,thp),one);
            assert rvCubicTreeTranslateBack(rvCubicTreeTranslate(outside)) == outside;
            assert rvCubicTreeTranslateBack(rvCubicTreeTranslate(inside)) == inside;
            assert rvCubicTreeTranslate(rvCubicTreeTranslateBack(across)) == across;
            assert rvCubicTreeTranslateBack(rvCubicTreeTranslate(across)) == across;
            let transformed = rvCubicTreeApply(g,across);
            let shiftedpair = rvCubicTreeTranslate(transformed);
            rvCubicTreeDistance(inside,transformed);
            rvCubicTreeDistance(rvCubicTreeTranslate(inside),shiftedpair);
            rvCubicTreeDistance(across,outside);
            rvCubicTreeDistance(rvCubicTreeApply(g,across),rvCubicTreeApply(g,outside));",
            setup(t, 4)
        )));
        let lines: Vec<_> = output.lines().collect();
        assert_eq!(lines[0], lines[1]);
        assert_eq!(lines[2], lines[3]);
    }
}

#[test]
fn geometric_types_and_precision_failures_are_specific() {
    let setup = setup(2, 2);
    for (statement, fragment) in [
        (
            "rvcubicdiv E[K,2,2]; show rvCubicTreeApply(g,rvCubicTreeRoot[E]);",
            "RVCubicTreeVertex(D)",
        ),
        (
            "rvcubicdiv E[K,2,2]; show rvCubicTreeLca(o,rvCubicTreeRoot[E]);",
            "RVCubicTreeVertex(D)",
        ),
        (
            "rvcubicdiv E[K,2,2]; show rvCubicTreeCompose(g,rvCubicTreePi[E]);",
            "RVCubicTreeOp(D)",
        ),
        (
            "show rvCubicTreeMatrix(one,th,one,th);",
            "invertible reduction",
        ),
        ("show rvCubicTreeLine[1](pi,nul);", "primitive right column"),
        ("show rvCubicTreeAt[D,3,0];", "radius exceeds"),
        (
            "show rvCubicTreeAt[D,2,756];",
            "sphere index is out of range",
        ),
        ("show rvCubicTreeAncestor[2](o);", "ancestor radius exceeds"),
        (
            "show rvCubicTreeChild[27](rvCubicTreeAt[D,1,4]);",
            "nonroot child index",
        ),
        (
            "show rvCubicTreeTranslate(rvCubicTreeAt[D,2,729]);",
            "one more outward digit",
        ),
        (
            "show rvCubicTreeTranslateBack(rvCubicTreeAt[D,2,0]);",
            "one more outward digit",
        ),
    ] {
        let err = diagnostic(interpret(&format!("{setup} {statement}")));
        assert!(err.contains(fragment), "{statement}\n{err}");
    }
}
