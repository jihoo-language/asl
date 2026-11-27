mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn potential_trace_and_actual_tangent_duality() {
    let output = stdout(example("deligne-ginzburg-calabi-yau"));
    for line in [
        "rvPartitionDGTrace[D] = []",
        "rvPartitionDGTrace[D] = [(0, [0, 0, 0], 2)]",
        "rvPartitionDGTangent[D] = complex=[4, 4, 4, 4], differentialRanks=[0, 0, 0], ext=[4, 4, 4, 4], tracePairingRanks=[4, 4, 4, 4]",
        "rvPartitionDGTangent[D] = complex=[4, 4, 4, 4], differentialRanks=[2, 2, 2], ext=[2, 0, 0, 2], tracePairingRanks=[2, 0, 0, 2]",
        "rvPartitionDGCohomology[D,0] = basis=[[1, 0, 0, 1], [0, 0, 1, 0]], dualBasis=[[0, 0, 1, 0], [0, 0, 0, 1]], tracePairing=[[0, 1], [1, 0]]",
        "rvPartitionDGPotentialTrace[C] = 0",
        "rvPartitionDGComplex[C] = complex=[2, 2, 2, 2], differentials=[[[2, 1], [0, 0]], [[0, 0], [0, 1]], [[2, 0], [1, 0]]]",
        "rvPartitionDGTangent[C] = complex=[2, 2, 2, 2], differentialRanks=[1, 1, 1], ext=[1, 0, 0, 1], tracePairingRanks=[1, 0, 0, 1]",
        "characteristic=3, symbols=[(0, 0, 0, 0), (1, 0, 0, -1), (2, 0, 0, -2)], cyclicComponents=[([0], 2, 3)], potential=false",
    ] {
        assert!(output.contains(line), "missing {line}:\n{output}");
    }
}

#[test]
fn signed_leibniz_and_nonzero_cyclic_cancellations_over_extensions() {
    let output = stdout(interpret(concat!(
        "rvunram U[2,2]{1,1,1}; rvunram V[3,1]{1,0,1};\n",
        "rvpartitionginzburg F[1,0,U]; rvpartitionginzburg G[2,0,V];\n",
        "let x=rvPartitionDGArrow[F,0]; let a=rvPartitionDGArrow[G,0];\n",
        "assert rvPartitionDGScale[2](rvPartitionDGScale[2](x))==rvPartitionDGScale[3](x);\n",
        "assert rvPartitionDGScale[3](rvPartitionDGScale[3](a))==neg(a);\n",
        "let aStar=rvPartitionDGStar[G,0]; let bStar=rvPartitionDGStar[G,1];\n",
        "assert rvPartitionDGDifferential(aStar)==rvPartitionDGPath[G,2]{1,0,1};\n",
        "assert rvPartitionDGDifferential(bStar)==rvPartitionDGPath[G,0]{0,1,0};\n",
        "let z=rvPartitionDGProduct(aStar,bStar);\n",
        "let dz=add(rvPartitionDGProduct(rvPartitionDGDifferential(aStar),bStar),neg(rvPartitionDGProduct(aStar,rvPartitionDGDifferential(bStar))));\n",
        "assert rvPartitionDGDifferential(z)==dz;\n",
        "show dz;\n",
        "assert rvPartitionDGDifferential(dz)==zero[RVPartitionDG(G,0)];\n",
        "rvPartitionDGTrace(add(z,rvPartitionDGProduct(bStar,aStar)));\n",
        "assert rvPartitionDGProduct(a,a)==zero[RVPartitionDG(G,0)];\n",
        "rvPartitionDGTrace(rvPartitionDGPath[F,0]{1,1});\n",
        "rvPartitionDGTangent[F](rvPartitionHallModule[1,0,U]{(0,2,1)});\n",
        "rvPartitionDGTangent[G](rvPartitionHallModule[2,0,V]{(0,1,1),(2,1,1)});\n",
    )));
    assert!(output.contains("RVPartitionDG(G,-1)[(2, [1, 0, 1, 3], 1), (2, [2, 0, 1, 0], 2)]"));
    assert!(output.contains("rvPartitionDGTrace[G] = []"));
    assert!(output.contains("rvPartitionDGTrace[F] = [(0, [1, 1], 1)]"));
    assert!(output.contains("rvPartitionDGTangent[F] = complex=[4, 4, 4, 4], differentialRanks=[2, 2, 2], ext=[2, 0, 0, 2]"));
    assert!(output.contains("rvPartitionDGTangent[G] = complex=[2, 2, 2, 2], differentialRanks=[0, 0, 0], ext=[2, 2, 2, 2]"));
}

#[test]
fn the_differential_squares_to_zero_on_360_actual_paths() {
    let mut source = String::from("rvpartitionginzburg L[1,0,5]; rvpartitionginzburg C[2,0,5];\n");
    let loop_graph = [(0, 0, 0), (0, 0, -1), (0, 0, -2)];
    let cyclic_graph = [
        (0, 2, 0),
        (2, 0, 0),
        (2, 0, -1),
        (0, 2, -1),
        (0, 0, -2),
        (1, 1, -2),
        (2, 2, -2),
    ];
    let mut paths = 0;
    for (frame, graph, starts) in [
        ("L", loop_graph.as_slice(), &[0][..]),
        ("C", cyclic_graph.as_slice(), &[0, 2][..]),
    ] {
        for &start in starts {
            let mut words = vec![(start, 0, Vec::new())];
            for _ in 1..=4 {
                words = words
                    .into_iter()
                    .flat_map(|(end, degree, letters)| {
                        graph
                            .iter()
                            .enumerate()
                            .filter(move |(_, (s, _, _))| *s == end)
                            .map(move |(i, &(_, t, d))| {
                                let mut extended = letters.clone();
                                extended.push(i);
                                (t, degree + d, extended)
                            })
                    })
                    .collect();
                for (_, degree, letters) in &words {
                    let terms = letters
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",");
                    writeln!(source,"assert rvPartitionDGDifferential(rvPartitionDGDifferential(rvPartitionDGPath[{frame},{start}]{{{terms}}}))==zero[RVPartitionDG({frame},{})];",degree+2).unwrap();
                    paths += 1;
                }
            }
        }
    }
    assert_eq!(paths, 360);
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn gradient_only_specialization_preserves_defining_relations_and_formal_tangent() {
    let output = stdout(interpret(concat!(
        "rvpartitionginzburggradient L[1,0,3]; rvpartitionginzburggradient C[2,0,2];\n",
        "assert rvPartitionDGDifferential(rvPartitionDGStar[L,0])==rvPartitionDGPath[L,0]{0,0};\n",
        "assert rvPartitionDGDerivative[C,0]==rvPartitionDGPath[C,2]{1,0,1};\n",
        "rvPartitionDGTangent[C](rvPartitionHallModule[2,0,2]{(0,2,1)});\n",
        "let M=rvPartitionHallModule[1,0,3]{(0,1,1),(0,2,1)};\n",
        "rvpartitionextensions E(M,M); rvpartitionformal D[E,1,2];\n",
        "rvPartitionFormalFrame[D]; rvPartitionDGTangent[L](M);\n",
        "rvpartitionginzburg S[2,3,2]; let V=rvPartitionHallModule[2,3,2]{(1,1,2)};\n",
        "assert rvPartitionDGPotential[S]==zero[RVPartitionDG(S,0)];\n",
        "assert rvPartitionDGDifferential(rvPartitionDGLoop[S,1])==zero[RVPartitionDG(S,-1)];\n",
        "rvPartitionDGTangent[S](V);\n",
    )));
    for line in [
        "differentialRanks=[1, 1, 1], ext=[1, 0, 0, 1]",
        "derivativeRank=4, tangent=5, selfExt1=1",
        "complex=[9, 9, 9, 9], differentialRanks=[4, 4, 4], ext=[5, 1, 1, 5]",
        "complex=[4, 0, 0, 4], differentialRanks=[0, 0, 0], ext=[4, 0, 0, 4]",
    ] {
        assert!(output.contains(line), "missing {line}:\n{output}");
    }
}

#[test]
fn incompatible_cycles_and_field_or_grade_mismatches_have_checked_errors() {
    let fixture =
        "rvpartitionginzburg D[1,0,5]; rvpartitionginzburg E[1,0,5]; rvpartitionginzburg P[2,0,5];";
    for (suffix, detail) in [
        ("rvpartitionginzburg Bad[1,0,3];", "denominator 3 vanishes in characteristic 3"),
        ("rvpartitionginzburg Bad[2,0,2];", "denominator 2 vanishes in characteristic 2"),
        ("rvpartitionginzburg Bad[2,1,5];", "uniform zero-relation lengths"),
        ("rvpartitionginzburggradient Bad[2,2,5];", "uniform zero-relation lengths"),
        ("rvpartitionginzburg D[1,0,5];", "already defined"),
        ("rvpartitionginzburggradient G[1,0,3]; show rvPartitionDGPotential[G];", "integral differentiated data only"),
        ("show rvPartitionDGProduct(rvPartitionDGArrow[D,0],rvPartitionDGArrow[E,0]);", "same nominal algebra"),
        ("show rvPartitionDGPath[P,0]{1};", "noncomposable symbols"),
        ("show rvPartitionDGScale[5](rvPartitionDGArrow[D,0]);", "outside the finite field codes"),
        ("rvPartitionDGCohomology[D,4](rvPartitionHallModule[1,0,5]{});", "degree must be in 0..=3"),
        ("rvPartitionDGTangent[D](rvPartitionHallModule[1,0,3]{});", "same rank, parameter, and coefficient field"),
        ("rvPartitionDGTangent[D](rvPartitionHallModule[1,0,5]{(0,1,7)});", "total length 6"),
        ("show rvPartitionDGDifferential(zero[RVPartitionDG(D,170141183460469231731687303715884105727)]);", "degree overflows signed 128 bits"),
    ] {
        let error = diagnostic(interpret(&format!("{fixture} {suffix}")));
        assert!(error.contains(detail), "{suffix}: {error}");
    }
    let long = std::iter::repeat_n("0", 65).collect::<Vec<_>>().join(",");
    let error = diagnostic(interpret(&format!(
        "{fixture} show rvPartitionDGPath[D,0]{{{long}}};"
    )));
    assert!(error.contains("word length exceeds 64"), "{error}");
}
