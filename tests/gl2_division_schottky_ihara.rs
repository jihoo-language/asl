mod common;

use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

fn fixtures(invariant: usize, height: usize) -> String {
    format!(
        "rvunram B[3,7]{{0,1}}; rvcyclotomiclocal K[B,1];
         rvcubicdiv D[K,{invariant},{height}];
         let o = rvCubicTreeRoot[D];
         let unit = rvCubicDivElem[D]([1]);
         let nul = zero[RVCubicDiv(D)];
         let theta = rvCubicDivTheta[D];
         let k0 = rvCubicTreeMatrix(unit,nul,nul,unit);
         let k1 = rvCubicTreeMatrix(unit,unit,unit,theta);"
    )
}

#[test]
fn runnable_marking_certifies_the_actual_axis_and_ihara_polynomial() {
    let output = stdout(example("gl2-division-schottky-ihara"));
    for expected in [
        "rvCubicSchottkyData[G] = division=D, rank=2, lengths=[1, 2], directions=[[27, 0], [3, 9]], ping_pong=certified",
        "rvCubicSchottkyLengths = letters=3, displacement=4, translation=2, prefix=[1], core=[2], power=1",
        "rvCubicSchottkyPath = radius=4, junctions=[0, 1, 3, 4]",
        "rvCubicSchottkyAxisCheck = vertex_radius=1, image_radius=3, displacement=2, verified=true",
        "rvCubicSchottkyPingPong = letter=1, image_radius=2, attracting_direction=27, verified=true",
        "rvCubicSchottkyCore[G] = vertices=2, edges=[(0, 0), (0, 1), (1, 0)], directed=6, rank=2",
        "rvCubicSchottkyNB[G,0,1] = 0",
        "rvCubicSchottkyWeighted[G,1,2] = u^2",
        "rvCubicIhara[G] = reciprocal=[1, -2, -1, 0, 3, 2, -3], bass_verified=true",
        "rvCubicIharaCoeff[G,2] = 5",
        "rvCubicSchottkyPrimitive[G,8] = 50",
        "rvCubicSchottkyEnumerate[G,8] = 50, verified=true",
    ] {
        assert!(output.contains(expected), "{expected}\n{output}");
    }
}

#[test]
fn both_invariants_act_on_vertices_and_boundary_with_proven_reduced_lengths() {
    for t in [1, 2] {
        let output = stdout(interpret(&format!(
            "{}
             rvcubicschottky G[D]{{[1,k0],[2,k1]}};
             let a = rvCubicSchottkyWord[G,[1]];
             let b = rvCubicSchottkyWord[G,[2]];
             let id = rvCubicSchottkyWord[G,[]];
             let w = rvCubicSchottkyWord[G,[1,2,-1]];
             let p = rvCubicSchottkyPow[-3](w);
             assert rvCubicSchottkyCyclic(w) == b;
             assert rvCubicSchottkyPrefix(w) == a;
             assert rvCubicSchottkyAxis(w) == rvCubicSchottkyApply(a,o);
             assert rvCubicSchottkyWord[G,[1,2,-2,-1,2]] == b;
             assert rvCubicSchottkyMul(w,rvCubicSchottkyInverse(w)) == id;
             assert rvCubicSchottkyPow[0](w) == id;
             assert rvCubicSchottkyPow[1000000000](id) == id;
             assert p == rvCubicSchottkyWord[G,[1,-2,-2,-2,-1]];
             assert rvCubicSchottkyIsPrimitive(w) == S.hi;
             assert rvCubicSchottkyIsPrimitive(p) == S.lo;
             assert rvCubicSchottkyIsPrimitive(id) == S.lo;
             let ab = rvCubicSchottkyMul(a,b);
             let ba = rvCubicSchottkyMul(b,a);
             assert rvCubicSchottkyConjugate(ab,ba) == S.hi;
             assert rvCubicSchottkyConjugate(ab,rvCubicSchottkyInverse(ab)) == S.lo;
             let v = rvCubicSchottkyApply(ab,o);
             assert v == rvCubicSchottkyApply(a,rvCubicSchottkyApply(b,o));
             assert rvCubicSchottkyApply(rvCubicSchottkyInverse(ab),v) == o;
             assert rvCubicTreeLca(v,rvCubicSchottkyApply(ba,o)) == o;
             let dir = rvCubicSchottkyDirection[G,2];
             let other = rvCubicSchottkyDirection[G,-1];
             let one = rvCubicBoundaryConstant[D,1,1];
             let f = rvCubicBoundaryCylinder(dir);
             let fwd = rvCubicSchottkyBoundary(a,f);
             assert fwd == rvCubicBoundaryCylinder(rvCubicSchottkyApply(a,dir));
             assert rvCubicSchottkyBoundary(rvCubicSchottkyInverse(a),fwd) == f;
             assert rvCubicSchottkyBoundary(a,rvCubicBoundaryCylinder(other)) ==
                 add(one,neg(rvCubicBoundaryCylinder(rvCubicSchottkyDirection[G,1])));
             assert rvCubicSchottkyBoundary(ab,f) ==
                 rvCubicSchottkyBoundary(a,rvCubicSchottkyBoundary(b,f));
             assert rvCubicSchottkyCoreVertex[G,0] == o;
             assert rvCubicSchottkyCoreVertex[G,1] ==
                 rvCubicTreeAncestor[1](rvCubicSchottkyApply(b,o));
             rvCubicSchottkyPath(id);
             rvCubicSchottkyPath(ab);
             rvCubicSchottkyPath(ba);
             rvCubicSchottkyPath(p);
             rvCubicSchottkyLengths(p);
             rvCubicSchottkyPingPong[G,-2](other);",
            fixtures(t, 10)
        )));
        for expected in [
            "rvCubicSchottkyPath = radius=0, junctions=[0], address=[], verified=true",
            "rvCubicSchottkyPath = radius=3, junctions=[0, 1, 3], address=[27",
            "rvCubicSchottkyPath = radius=3, junctions=[0, 2, 3], address=[3",
            "rvCubicSchottkyPath = radius=8, junctions=[0, 1, 3, 5, 7, 8]",
            "rvCubicSchottkyLengths = letters=5, displacement=8, translation=6, prefix=[1], core=[-2, -2, -2], power=3",
            "rvCubicSchottkyPingPong = letter=-2, image_radius=3, attracting_direction=9, verified=true",
        ] {
            assert!(output.contains(expected), "{expected}\n{output}");
        }
    }
}

#[test]
fn primitive_oriented_cycles_match_bass_nonbacktracking_and_direct_counts() {
    let mut program = fixtures(1, 4);
    program.push_str(
        "
       rvcubicschottky G[D]{[1,k0],[1,k1]};
       rvcubicschottky H[D]{[1,k0],[2,k1]};
       rvCubicIhara[G]; rvCubicIhara[H];
       rvCubicSchottkyTrace[G,0]; rvCubicSchottkyTrace[H,2];
       rvCubicSchottkyNB[H,0,2]; rvCubicSchottkyNB[H,2,4];
       rvCubicSchottkyWeighted[H,2,-2];
    ",
    );
    for name in ["G", "H"] {
        for n in 1..=8 {
            writeln!(
                program,
                "rvCubicSchottkyPrimitive[{name},{n}]; rvCubicSchottkyEnumerate[{name},{n}];"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&program));
    for expected in [
        "rvCubicIhara[G] = reciprocal=[1, -4, 2, 4, -3], bass_verified=true",
        "rvCubicIhara[H] = reciprocal=[1, -2, -1, 0, 3, 2, -3], bass_verified=true",
        "rvCubicSchottkyTrace[G,0] = 4",
        "rvCubicSchottkyTrace[H,2] = 6",
        "rvCubicSchottkyNB[H,0,2] = 1",
        "rvCubicSchottkyNB[H,2,4] = 1",
        "rvCubicSchottkyWeighted[H,2,-2] = 0",
    ] {
        assert!(output.contains(expected), "{expected}\n{output}");
    }
    let groups = [
        ("G", [4, 4, 8, 18, 48, 116, 312, 810]),
        ("H", [2, 2, 4, 4, 8, 14, 28, 50]),
    ];
    let expected: Vec<_> = groups
        .into_iter()
        .flat_map(|(name, counts)| {
            counts.into_iter().enumerate().flat_map(move |(n, c)| {
                let n = n + 1;
                [
                    format!("rvCubicSchottkyPrimitive[{name},{n}] = {c}"),
                    format!("rvCubicSchottkyEnumerate[{name},{n}] = {c}, verified=true"),
                ]
            })
        })
        .collect();
    assert_eq!(output.lines().skip(7).collect::<Vec<_>>(), expected);
}

#[test]
fn empty_and_rank_one_cores_and_long_formal_lengths() {
    let output = stdout(interpret(&format!(
        "{}
         rvcubicschottky Empty[D]{{}};
         rvcubicschottky Line[D]{{[2,k0]}};
         rvcubicschottky Long[D]{{[300,k0]}};
         assert rvCubicSchottkyApply(rvCubicSchottkyWord[Empty,[]],o) == o;
         rvCubicIhara[Empty]; rvCubicIharaCoeff[Empty,0];
         rvCubicSchottkyEnumerate[Empty,3];
         rvCubicIhara[Line]; rvCubicSchottkyEnumerate[Line,2];
         rvCubicSchottkyEnumerate[Line,4];
         rvCubicSchottkyLengths(rvCubicSchottkyPow[2](rvCubicSchottkyWord[Long,[1]]));",
        fixtures(2, 3)
    )));
    assert_eq!(output, concat!(
        "rvCubicIhara[Empty] = reciprocal=[1], bass_verified=true\n",
        "rvCubicIharaCoeff[Empty,0] = 1\n",
        "rvCubicSchottkyEnumerate[Empty,3] = 0, verified=true\n",
        "rvCubicIhara[Line] = reciprocal=[1, 0, -2, 0, 1], bass_verified=true\n",
        "rvCubicSchottkyEnumerate[Line,2] = 2, verified=true\n",
        "rvCubicSchottkyEnumerate[Line,4] = 0, verified=true\n",
        "rvCubicSchottkyLengths = letters=2, displacement=600, translation=600, prefix=[], core=[1, 1], power=2\n",
    ));
}

#[test]
fn invalid_markings_and_unsupported_precisions_remain_distinct() {
    let base = format!(
        "{} rvcubicschottky G[D]{{[1,k0],[2,k1]}};
      rvcubicschottky H[D]{{[1,k0],[1,k1]}};
      let a=rvCubicSchottkyWord[G,[1]];
      let b=rvCubicSchottkyWord[H,[1]];",
        fixtures(2, 3)
    );
    for (input, message) in [
        ("rvcubicschottky Bad[D]{[0,k0]};", "lengths must be positive"),
        ("rvcubicschottky Bad[D]{[1,k0],[1,k0]};", "signed first directions to be distinct"),
        ("rvcubicschottky Bad[D]{[1,rvCubicTreePi[D]]};", "integral GL2(O_D)"),
        ("rvcubicschottky G[D]{[1,k0]};", "already defined"),
        ("show rvCubicSchottkyWord[G,[0]];", "outside the signed generator range"),
        ("show rvCubicSchottkyWord[G,[-3]];", "outside the signed generator range"),
        ("show rvCubicSchottkyMul(a,b);", "RVCubicSchottkyWord(G)"),
        ("show rvCubicSchottkyPow[20001](a);", "at most 20000 letters"),
        ("rvCubicSchottkyPingPong[G,1](o);", "needs a nonroot ray"),
        ("rvCubicSchottkyPingPong[G,1](rvCubicSchottkyDirection[G,-1]);", "excluded repelling"),
        ("show rvCubicSchottkyCoreVertex[G,2];", "core vertex index is out of range"),
        ("rvCubicSchottkyNB[G,6,0];", "directed edge index is out of range"),
        ("rvcubicschottky J[D]{[300,k0]}; rvCubicIhara[J];", "at most 256 edges"),
        ("rvcubicschottky J[D]{[300,k0]}; show rvCubicSchottkyApply(rvCubicSchottkyWord[J,[1]],o);", "requires more lattice precision"),
        ("rvCubicSchottkyPath(rvCubicSchottkyPow[4](a));", "root path requires more lattice precision"),
        ("rvcubicdiv E[K,2,3]; show rvCubicSchottkyApply(a,rvCubicTreeRoot[E]);", "RVCubicTreeVertex(D)"),
    ] {
        let error = diagnostic(interpret(&format!("{base} {input}")));
        assert!(error.contains(message), "{input}\n{error}");
    }
}
