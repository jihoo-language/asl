mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn same_three_edges_and_different_residue_cross_ratio() {
    let output = stdout(example("rv-stable-marks"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(&lines[..4], [
        "rvStableSplits = [([0, 1, 2, 3], [4, 5, 6, 7, 8], 19), ([0, 1, 2, 3, 6, 7, 8], [4, 5], 11), ([0, 2, 4, 5, 6, 7, 8], [1, 3], 7)]",
        "rvStableDirections[1] = [[0], [1, 3], [2], [4, 5, 6, 7, 8]]",
        "rvStableResidue[1, 8] = 2",
        "rvStableResidue[1, 8] = 5",
    ]);
    assert_eq!(lines[4], concat!(
        "rvStableComponents = [[([0, 1, 2, 3], infinity), ([4, 5], 0), ([6], 1), ([7], 6), ([8], 2)], ",
        "[([0], infinity), ([1, 3], 0), ([2], 1), ([4, 5, 6, 7, 8], 2)], ",
        "[([0, 2, 4, 5, 6, 7, 8], infinity), ([1], 0), ([3], 1)], ",
        "[([0, 1, 2, 3, 6, 7, 8], infinity), ([4], 0), ([5], 1)]]",
    ));
}

fn rational(points: &[[i128; 2]], prime: u32) -> String {
    let mut value = format!("rvStableLine[{prime}]{{");
    for (i, &[x, y]) in points.iter().enumerate() {
        if i != 0 {
            value.push(',');
        }
        write!(value, "[{x},{y}]").unwrap();
    }
    value.push('}');
    value
}

fn matroid(points: &[[i128; 2]], prime: u32) -> String {
    let mut value = format!("rvValMat[{prime},{}]{{", points.len());
    for (axis, join) in [(0, "["), (1, ",[")] {
        value.push_str(join);
        for (i, column) in points.iter().enumerate() {
            if i != 0 {
                value.push(',');
            }
            write!(value, "{}", column[axis]).unwrap();
        }
        value.push(']');
    }
    value.push('}');
    value
}

#[test]
fn mobius_change_and_independent_signed_column_scaling_preserve_every_component() {
    let a = 7_i128.pow(19);
    let points = [
        [0, 1],
        [a, 1],
        [3 * a, 1],
        [a * (1 + 7_i128.pow(7)), 1],
        [1, 1],
        [1 + 7_i128.pow(11), 1],
        [2, 1],
        [3, 1],
        [1, 0],
    ];
    let moved: Vec<_> = points
        .iter()
        .enumerate()
        .map(|(i, &[x, y])| {
            let scale = 7_i128.pow((i % 5) as u32) * if i % 2 == 0 { 1 } else { -3 };
            [scale * y, scale * (x + y)]
        })
        .collect();
    let program = format!(
        "let C={}; let D={}; assert C==D; assert rvStableTree(C)==rvStableTree({});",
        rational(&points, 7),
        rational(&moved, 7),
        matroid(&moved, 7)
    );
    assert_eq!(stdout(interpret(&program)), "");
}

#[test]
fn quadratic_cancellation_and_cross_field_agreement() {
    let program = r"
        rvqpadic K[5,-1,2];
        let z=zero[RVQpQuad(K)];
        let o=rvQpQuad[K,1,0,1];
        let two=rvQpQuad[K,2,0,1];
        let b=rvQpQuad[K,-1,-2,1];
        let c=add(o,rvQpQuadPow[19](b));
        let A=rvStableLineQuad[K]{[o,z],[z,o],[o,o],[c,o],[two,o]};
        let B=rvStableLineQuad[K]{[b,z],[z,rvQpQuadInv(b)],[o,o],[c,o],[two,o]};
        assert A==B;
        let T=rvStableTree(A);
        rvStableSplits(T);
        rvStableResidue[0,4](A);
        assert T==rvStableTree(rvValMatQuad[K,5]{[o,z,o,c,two],[z,o,o,o,o]});
    ";
    let mut input = program.to_owned();
    writeln!(
        input,
        "assert A==rvStableLine[5]{{[1,0],[0,1],[1,1],[{},1],[2,1]}};",
        1 + 5_i128.pow(19)
    )
    .unwrap();
    assert_eq!(
        stdout(interpret(&input)),
        "rvStableSplits = [([0, 1, 4], [2, 3], 19)]\nrvStableResidue[0, 4] = 2\n"
    );
}

#[test]
fn anchor_relabeling_three_marks_and_large_edge_lengths() {
    let a = 1_i128 << 120;
    let p = format!(
        r"
        let A=rvStableLine[2]{{[1,0],[0,1],[{a},1],[1,1]}};
        let B=rvStableLine[2]{{[0,1],[{a},1],[1,1],[1,0]}};
        rvStableSplits(rvStableTree(A));
        rvStableSplits(rvStableTree(B));
        let C=rvStableLine[2]{{[1,0],[0,1],[1,1]}};
        let D=rvStableLine[2]{{[4,5],[1,9],[2,3]}};
        assert C==D;
        rvStableDirections[0](rvStableTree(C));
        rvStableResidue[0,0](C);
        rvStableResidue[0,1](C);
        rvStableResidue[0,2](C);
        card (\t:RVStableTree(3) => S.hi);
        card (\t:RVStableMarks(2,3) => S.hi);
    "
    );
    assert_eq!(stdout(interpret(&p)),
        "rvStableSplits = [([0, 3], [1, 2], 120)]\nrvStableSplits = [([0, 1], [2, 3], 120)]\nrvStableDirections[0] = [[0], [1], [2]]\nrvStableResidue[0, 0] = infinity\nrvStableResidue[0, 1] = 0\nrvStableResidue[0, 2] = 1\ncard = 1 / 1\ncard = 1 / 1\n");
}

#[test]
fn geometry_and_queries_validate_projective_distinctness_and_types() {
    let setup = "let C=rvStableLine[5]{[1,0],[0,1],[1,1]}; let T=rvStableTree(C);";
    for (program, expected) in [
        ("show rvStableLine[2]{[1,0],[0,1]};", "at least three marks"),
        (
            "show rvStableLine[3]{[0,0],[1,0],[0,1]};",
            "mark 0 has zero homogeneous coordinates",
        ),
        (
            "show rvStableLine[3]{[1,0],[0,1],[-9,0]};",
            "marks 0 and 2 collide exactly",
        ),
        (
            "show rvStableLine[3]{[1,0],[1,1],[2]};",
            "requires two homogeneous coordinates",
        ),
        (
            "show rvStableTree(rvValMat[3,3]{[1,1,0],[0,0,1]});",
            "full-support rank-two",
        ),
        (
            "show rvStableTree(rvValMat[3,4]{[1,0,0,0]});",
            "rank-two RVValMat",
        ),
        ("rvStableDirections[1](T);", "unknown stable vertex 1"),
        ("rvStableResidue[0,3](C);", "outside the ground set"),
        ("rvStableComponents(T);", "requires an RVStableMarks"),
        ("rvStableSplits(C);", "requires an RVStableTree"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{program}")));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
