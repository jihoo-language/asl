mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn example_detects_quadratic_and_independent_cubic_monomials() {
    let output = stdout(example("deligne-formal-deformations"));
    assert_eq!(output,concat!(
        "rvPartitionFormalObstruction[D] = nextOrder=2, residual=[([1, 1], [1, 0, 0, 1])], obstructions=[([1, 1], [1, 0, 0, 0], 1)], lifts=0, liftsModuloNewGauge=0\n",
        "rvPartitionFormalObstruction[L] = nextOrder=2, residual=[], obstructions=[], lifts=729, liftsModuloNewGauge=729\n",
        "rvPartitionFormalObstruction[L] = nextOrder=3, residual=[([1, 2], [0, 1]), ([2, 1], [1, 0])], obstructions=[([1, 2], [0, 1], 1), ([2, 1], [1, 0], 1)], lifts=0, liftsModuloNewGauge=0\n",
        "RVPartitionFormal(C) order=3, corrections=[([1], [1, 0, 0, 2]), ([2], [0, 2, 0, 0])]\n",
    ));
}

#[test]
fn every_small_first_matrix_matches_the_exact_square_zero_equations() {
    let mut input = String::new();
    let mut expected = Vec::new();
    for q in [2, 3, 5] {
        writeln!(input,"let M{q}=rvPartitionHallModule[1,0,{q}]{{(0,1,2)}}; rvpartitionextensions E{q}(M{q},M{q}); rvpartitionformal D{q}[E{q},1,2];").unwrap();
        let mut unobstructed = 0;
        for a in 0..q {
            for b in 0..q {
                for c in 0..q {
                    for d in 0..q {
                        writeln!(input,"rvPartitionFormalObstruction(rvPartitionFormalStart[D{q}](rvPartitionExtension[E{q}]{{{a},{b},{c},{d}}}));").unwrap();
                        let residual = [
                            (a * a + b * c) % q,
                            b * (a + d) % q,
                            c * (a + d) % q,
                            (d * d + b * c) % q,
                        ];
                        let display = if residual.iter().all(|&x| x == 0) {
                            unobstructed += 1;
                            "residual=[], obstructions=[]".into()
                        } else {
                            format!("residual=[([2], {residual:?})], obstructions=[")
                        };
                        expected.push(display);
                    }
                }
            }
        }
        assert_eq!(unobstructed, q * q);
    }
    let output = stdout(interpret(&input));
    assert_eq!(output.lines().count(), expected.len());
    for (line, residual) in output.lines().zip(expected) {
        assert!(line.contains(&residual), "{line}");
    }
}

#[test]
fn lookahead_and_user_corrections_distinguish_second_and_third_order() {
    let mut source = String::new();
    let mut expected = String::new();
    for q in [2, 3, 5] {
        writeln!(source,"let M{q}=rvPartitionHallModule[1,0,{q}]{{(0,1,2)}}; rvpartitionextensions E{q}(M{q},M{q}); rvpartitionformal D{q}[E{q},1,3];").unwrap();
        writeln!(
            source,
            "let z{q}=rvPartitionFormalStart[D{q}](zero[RVPartitionExtension(E{q})]);"
        )
        .unwrap();
        writeln!(
            source,
            "let a{q}=rvPartitionFormalStart[D{q}](rvPartitionExtension[E{q}]{{0,1,0,0}});"
        )
        .unwrap();
        writeln!(
            source,
            "rvPartitionFormalLookahead(z{q});\nrvPartitionFormalLookahead(a{q});"
        )
        .unwrap();
        for count in [q * q * q * q, q * q] {
            writeln!(
                expected,
                "rvPartitionFormalLookahead[D{q}] = order=2, lifts={}, extendable={count}",
                q * q * q * q
            )
            .unwrap();
        }
    }
    source.push_str(concat!(
        "let bad=rvPartitionFormalCorrect(a3){([2],[1,0,0,0])};\n",
        "rvPartitionFormalObstruction(bad);\n",
        "let good=rvPartitionFormalCorrect(a3){([2],[1,0,0,2])};\n",
        "let lifted=rvPartitionFormalLift(good);\n",
        "assert rvPartitionFormalFirst[0](lifted)==rvPartitionExtension[E3]{0,1,0,0};\n",
    ));
    expected.push_str("rvPartitionFormalObstruction[D3] = nextOrder=3, residual=[([3], [0, 1, 0, 0])], obstructions=[([3], [0, 1, 0, 0], 1)], lifts=0, liftsModuloNewGauge=0\n");
    assert_eq!(stdout(interpret(&source)), expected);
}

#[test]
fn nonsemisimple_derivative_solves_a_nonzero_residual_and_gauge_is_exact() {
    let output = stdout(interpret(concat!(
        "let J=rvPartitionHallModule[1,0,3]{(0,2,1)};\n",
        "rvpartitionextensions E(J,J);\n",
        "rvpartitionformal G[E,1,3];\n",
        "rvPartitionFormalFrame[G];\n",
        "let first=rvPartitionFormalCorrect(rvPartitionFormalConstant[G]){([1],[1,0,0,2])};\n",
        "rvPartitionFormalObstruction(first);\n",
        "let lifted=rvPartitionFormalLift(rvPartitionFormalLift(first));\n",
        "let trivial=rvPartitionFormalTrivial[G];\n",
        "let g=rvPartitionFormalGauge(trivial){([1],[0,1,0,0])};\n",
        "assert lifted==g;\n",
        "assert rvPartitionFormalGauge(g){([1],[0,2,0,0])}==trivial;\n",
        "assert rvPartitionFormalSlice(g)==trivial;\n",
        "let diagonal=rvPartitionFormalGauge(trivial){([1],[1,0,0,0])};\n",
        "show diagonal;\n",
        "assert rvPartitionFormalFirst[0](diagonal)==zero[RVPartitionExtension(E)];\n",
        "assert rvPartitionFormalSlice(diagonal)==trivial;\n",
    )));
    assert_eq!(output,concat!(
        "rvPartitionFormalFrame[G] = q=3, axes=1, order=3, arrowEntries=4, gaugeEntries=4, relationEntries=4, derivativeRank=2, tangent=2, selfExt1=0\n",
        "rvPartitionFormalObstruction[G] = nextOrder=2, residual=[([2], [1, 0, 0, 1])], obstructions=[], lifts=9, liftsModuloNewGauge=1\n",
        "RVPartitionFormal(G) order=3, corrections=[([1], [0, 0, 2, 0]), ([2], [0, 0, 1, 0]), ([3], [0, 0, 2, 0])]\n",
    ));
}

#[test]
fn named_field_coordinates_multiply_in_nontrivial_mixed_and_gauge_terms() {
    let output=stdout(interpret(concat!(
        "rvunram F[2,1]{1,1,1};\n",
        "let M=rvPartitionHallModule[1,0,F]{(0,1,2)};\n",
        "rvpartitionextensions E(M,M);\n",
        "rvpartitionformal D[E,2,3];\n",
        "let a=rvPartitionExtension[E]{0,2,0,0};\n",
        "let b=rvPartitionExtension[E]{0,0,3,0};\n",
        "rvPartitionFormalObstruction(rvPartitionFormalStart[D](a,b));\n",
        "let v=rvPartitionFormalLift(rvPartitionFormalLift(rvPartitionFormalStart[D](a,zero[RVPartitionExtension(E)])));\n",
        "let moved=rvPartitionFormalGauge(v){([0,1],[0,0,3,0])};\n",
        "assert rvPartitionFormalFirst[0](moved)==a;\n",
        "assert rvPartitionFormalFirst[1](moved)==zero[RVPartitionExtension(E)];\n",
        "let J=rvPartitionHallModule[1,0,F]{(0,2,1)};\n",
        "rvpartitionextensions G(J,J);\n",
        "rvpartitionformal C[G,1,3];\n",
        "let w=rvPartitionFormalGauge(rvPartitionFormalTrivial[C]){([1],[0,2,0,0])};\n",
        "show w;\n",
        "assert rvPartitionFormalSlice(w)==rvPartitionFormalTrivial[C];\n",
    )));
    assert_eq!(output,concat!(
        "rvPartitionFormalObstruction[D] = nextOrder=2, residual=[([1, 1], [1, 0, 0, 1])], obstructions=[([1, 1], [1, 0, 0, 0], 1)], lifts=0, liftsModuloNewGauge=0\n",
        "RVPartitionFormal(C) order=3, corrections=[([1], [2, 0, 0, 2]), ([2], [0, 3, 0, 0])]\n",
    ));
}

#[test]
fn malformed_degrees_gauge_terms_and_nominal_frames_have_direct_diagnostics() {
    let fixture="let M=rvPartitionHallModule[1,0,3]{(0,1,2)}; rvpartitionextensions E(M,M); rvpartitionformal D[E,2,3]; let a=rvPartitionExtension[E]{0,1,0,0}; let b=rvPartitionExtension[E]{0,0,1,0}; let x=rvPartitionFormalStart[D](a,b);";
    for (suffix,expected) in [
        ("rvpartitionformal B[E,0,3];","require 1..=4 variables"),
        ("rvpartitionformal B[E,1,7];","order 1..=6"),
        ("rvpartitionformal D[E,1,3];","already defined"),
        ("show rvPartitionFormalStart[D](a);","one direction per variable"),
        ("rvpartitionextensions F(M,M); show rvPartitionFormalStart[D](a,rvPartitionExtension[F]{0,0,0,0});","wrong nominal self-extension space"),
        ("show rvPartitionFormalLift(x);","obstructed at [1, 1]"),
        ("show rvPartitionFormalCorrect(x){([1,1],[0,0,0,0])};","relation is nonzero at [1, 1]"),
        ("show rvPartitionFormalCorrect(x){([1,0],[0,0,0,0])};","invalid total degree"),
        ("show rvPartitionFormalGauge(x){([0,0],[0,0,0,0])};","invalid total degree"),
        ("show rvPartitionFormalGauge(x){([1,0],[0,0,0,3])};","outside the finite field codes"),
        ("show rvPartitionFormalGauge(x){([1,0],[1,0,0,0]),([1,0],[0,1,0,0])};","repeats a monomial"),
        ("show rvPartitionFormalFirst[2](x);","variable index is out of range"),
        ("rvPartitionFormalObstruction(rvPartitionFormalTrivial[D]);","already at the frame order"),
        ("let N=rvPartitionHallModule[1,0,3]{(0,1,1)};rvpartitionextensions F(M,N);rvpartitionformal K[F,1,2];","requires a self-extension space"),
    ] {
        let message=diagnostic(interpret(&format!("{fixture} {suffix}")));
        assert!(message.contains(expected),"{suffix}: {message}");
    }
}

#[test]
fn the_two_arrow_cubic_has_exactly_two_q_minus_one_admissible_directions() {
    let mut source = String::new();
    let mut expected = Vec::new();
    for q in [2, 3, 5] {
        writeln!(source,"let M{q}=rvPartitionHallModule[2,0,{q}]{{(0,1,1),(2,1,1)}}; rvpartitionextensions E{q}(M{q},M{q}); rvpartitionformal D{q}[E{q},1,3];").unwrap();
        let mut admitted = 0;
        for a in 0..q {
            for b in 0..q {
                writeln!(source,"rvPartitionFormalObstruction(rvPartitionFormalLift(rvPartitionFormalStart[D{q}](rvPartitionExtension[E{q}]{{{a},{b}}})));" ).unwrap();
                let values = [a * a * b % q, a * b * b % q];
                if values == [0, 0] {
                    admitted += 1;
                    expected.push("residual=[], obstructions=[]".into());
                } else {
                    expected.push(format!("residual=[([3], {values:?})], obstructions=["));
                }
            }
        }
        assert_eq!(admitted, 2 * q - 1);
    }
    let output = stdout(interpret(&source));
    assert_eq!(output.lines().count(), expected.len());
    for (line, residual) in output.lines().zip(expected) {
        assert!(line.contains(&residual), "{line}");
    }
}

#[test]
fn a_nonzero_derivative_image_can_leave_a_detectable_cokernel_component() {
    let output=stdout(interpret(concat!(
        "let M=rvPartitionHallModule[1,0,3]{(0,1,1),(0,2,1)};\n",
        "rvpartitionextensions E(M,M);\n",
        "rvpartitionformal D[E,1,2];\n",
        "rvPartitionFormalFrame[D];\n",
        "let x=rvPartitionFormalCorrect(rvPartitionFormalConstant[D]){([1],[1,0,0,0,1,0,0,0,2])};\n",
        "rvPartitionFormalObstruction(x);\n",
    )));
    assert_eq!(output,concat!(
        "rvPartitionFormalFrame[D] = q=3, axes=1, order=2, arrowEntries=9, gaugeEntries=9, relationEntries=9, derivativeRank=4, tangent=5, selfExt1=1\n",
        "rvPartitionFormalObstruction[D] = nextOrder=2, residual=[([2], [1, 0, 0, 0, 1, 0, 0, 0, 1])], obstructions=[([2], [1, 0, 0, 0, 0, 0, 0, 0, 0], 1)], lifts=0, liftsModuloNewGauge=0\n",
    ));
}
