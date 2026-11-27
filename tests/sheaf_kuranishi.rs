mod common;

use common::{diagnostic, example, interpret, stdout};

const ARMS: &str = include_str!("../examples/sheaf-kuranishi.asl");

#[test]
fn unequal_parallel_arms_give_one_equation_and_a_nonclosed_higher_defect() {
    let output = stdout(example("sheaf-kuranishi"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        lines[0],
        "shkuranishiDims[K] = tangent=6, obstruction=1, order=4"
    );
    assert_eq!(lines[1], "shkuranishiEq[K, 0] = x0*x1 + 4*x2*x3*x4");
    assert_eq!(lines[2], "ShHomC(F, F, 3){[a<l<b<c]#0->0}");
    assert_eq!(lines[3], "ShHomC(F, F, 1){[a<c]#0->0}");
    assert_eq!(lines[4], "shkuranishiEq[Diagonal, 0] = x0^2");
}

#[test]
fn the_two_path_relation_and_ramified_arc_work_in_small_characteristics() {
    for p in [2, 3, 5] {
        let source = ARMS.replace("sheaf F[X,5]", &format!("sheaf F[X,{p}]"));
        let output = stdout(interpret(&source));
        assert!(output.contains("ShHomC(F, F, 1){[a<c]#0->0}"));
        assert!(output
            .lines()
            .nth(2)
            .unwrap()
            .starts_with("ShHomC(F, F, 3){"));
        assert_ne!(output.lines().nth(2).unwrap(), "ShHomC(F, F, 3){0}");

        let blocked = diagnostic(interpret(&format!(
            "{source} show shKurArc[K,[1,1,1,1,1,1],[1,1,1,1,1,1],4];"
        )));
        assert!(blocked.contains("obstruction in order 2"), "{blocked}");
    }
}

#[test]
fn an_unobstructed_matrix_arc_agrees_with_direct_multiplication() {
    for p in [2, 3, 5] {
        let source = format!(
            r"
            space M {{a<b<c}};
            sheaf E[M,{p}] {{a:2,b:2,c:2}} {{a<b:[1,0,0,0],b<c:[1,0,0,0]}};
            let left = shHomTerm[E,E,[a,b],1,1];
            let right = shHomTerm[E,E,[b,c],1,1];
            shkuranishi C[E,4] {{shClass(left),shClass(right)}};
            assert shKurA[C,[1,0]] == left;
            assert shKurA[C,[0,1]] == right;
            assert shKurA[C,[1,1]] == shHomTerm[E,E,[a,c],1,1];
            assert shKurR[C,[1,1]] == zero[ShHomC(E,E,2)];
            let arc = shKurArc[C,[-1,1],[1,2],4];
            assert shDefTerm[1](arc) == shScale[-1](left);
            assert shDefTerm[2](arc) == right;
            assert shDefTerm[3](arc) == shHomTerm[E,E,[a,c],1,1,-1];
            assert shDefTerm[4](arc) == zero[ShHomC(E,E,1)];
            assert shDefObstruction(arc) == zero[ShExt(E,E,2)];
            shkuranishiDims[C];
            "
        );
        assert_eq!(
            stdout(interpret(&source)),
            "shkuranishiDims[C] = tangent=2, obstruction=0, order=4\n"
        );
    }
}

#[test]
fn active_degree_scheduling_avoids_inactive_monomials() {
    let input = format!(
        "{ARMS}
        shkuranishi Flat[F,1000000000] {{x0}};
        assert shKurA[Flat,[1000000000]] == zero[ShHomC(F,F,1)];
        assert shKurKappa[Flat,[1000000000]] == zero[ShExt(F,F,2)];
        shkuranishiEq[Flat,0];
        shkuranishi Full[F,1];
        shkuranishiDims[Full];
        "
    );
    assert!(stdout(interpret(&input)).ends_with(concat!(
        "shkuranishiEq[Flat, 0] = 0\n",
        "shkuranishiDims[Full] = tangent=6, obstruction=1, order=1\n",
    )));
}

#[test]
fn errors_distinguish_wrong_directions_arcs_and_truncation() {
    for (extra, expected) in [
        (
            "shkuranishi Bad[F,2] {x0,shScale[2](x0)};",
            "linearly independent",
        ),
        ("shkuranishi Bad[F,2] {shRep(x0)};", "degree-one self-Ext"),
        ("show shKurA[K,[0,1]];", "wrong number of directions"),
        (
            "show shKurR[K,[5,0,0,0,0,0]];",
            "beyond the computed truncation",
        ),
        ("shkuranishiEq[K,1];", "equation coordinate is out of range"),
        (
            "show shKurArc[K,[1,1],[1,1],2];",
            "one amplitude and weight",
        ),
        (
            "show shKurArc[K,[1,0,0,0,0,0],[0,1,1,1,1,1],2];",
            "weights must be positive",
        ),
        (
            "show shKurArcKappa[K,[0,0,0,0,0,0],[1,1,1,1,1,1],5];",
            "beyond the computed truncation",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{ARMS} {extra}")));
        assert!(error.contains(expected), "{extra}: {error}");
    }
}
