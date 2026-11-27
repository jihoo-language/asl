mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/hyperelliptic-kummer-cm-cycles.asl");

#[test]
fn gaussian_graph_classes_match_independent_f37_intersections_and_defects() {
    let output = stdout(example("hyperelliptic-kummer-cm-cycles"));
    for expected in [
        "product_gram=[[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 2], [1, 1, 2, 0]], product_discriminant=-4, tate_rank=4, locally_saturated=true",
        "gaussian_columns=[[1, 8], [2, 8]], weil_scale=1",
        "cycle_lifts=[[0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 1, 1, 1, 1, 0], [2, 1, 8, 1, 1, 1]]",
        "rvKummerCMCyclePair = product=1, normalized_kummer=1/2",
        "rvKummerCMCyclePair = product=2, normalized_kummer=1",
        "rvKummerCMCyclePair = product=0, normalized_kummer=0",
        "rvKummerCMCycleHermitian = a=1, d=1, b_in_1_omega=[0, 1]",
        "rvKummerCMDefect[A, 1] = elliptic_frobenius=(1, ±6), prime=3, transcendental_smith_exponents=[1, 1]",
        "rvKummerCMDefect[A, 2] = elliptic_frobenius=(1, ±6), prime=3, transcendental_smith_exponents=[1, 1]",
        "rvKummerCMDefect[A, 3] = elliptic_frobenius=(1, ±6), prime=3, transcendental_smith_exponents=[2, 2]",
        "rvKummerCMDefect[A, 9] = elliptic_frobenius=(1, ±6), prime=3, transcendental_smith_exponents=[3, 3]",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn genuine_cycle_and_projection_follow_a_nonproduct_marking_and_a_new_weil_root() {
    let output = stdout(interpret(&format!(
        "{EXAMPLE}
         rvjacpowerbasis BT[SL,3,2](rvUnramElem[L]{{16}}){{e1,e2,rvJacScale[2](add(f1,e1)),rvJacScale[2](f2)}};
         rvkummerbrauerpairing T[BT,S,1];
         rvkummercmcycles AT[T,EB](rvUnramElem[K]{{6}});
         rvKummerCMCycleData[AT];
         let tg=rvKummerCMCycle[AT]{{0,0,0,1}};
         assert rvKummerBrauerChange[T](rvKummerCMCycleLift[2](graph))==rvKummerCMCycleLift[2](tg);
         let tw=rvKummerBrauerChange[T](w);
         assert rvKummerBrauerChange[T](algebraic)==rvKummerCMAlgebraic[AT](tw);
         assert rvKummerBrauerChange[T](transcendental)==rvKummerCMTranscendental[AT](tw);
         rvkummerbrauerquotient QT[T,1];
         assert rvKummerBrauerClass[QT](rvKummerCMCycleLift[1](tg))==zero[RVKummerBrauerClass(QT)];"
    )));
    let data = output
        .lines()
        .find(|row| row.starts_with("rvKummerCMCycleData[AT]"))
        .unwrap();
    assert!(data.contains("weil_scale=5"), "{data}");
}

#[test]
fn nonmaximal_witness_remembers_its_intersections_and_rejects_bad_prime_projection() {
    let setup = format!(
        "{EXAMPLE}
         rvkummercmcycles U[P,EB,3](rvUnramElem[K]{{6}});
         rvKummerCMCycleData[U];
         let ug=rvKummerCMCycle[U]{{0,0,0,1}};
         let ud=rvKummerCMCycle[U]{{0,0,1,0}};
         rvKummerCMCyclePair(ug,ud);
         rvKummerCMCycleHermitian(add(ug,ud));
         assert rvKummerCMCycleLift[1](ug)==rvKummerCMCycleLift[1](left);
         rvkummercmcycles V[P,EB,2](rvUnramElem[K]{{6}});
         rvKummerCMCycleData[V];
         assert rvKummerCMAlgebraic[V](w)==algebraic;"
    );
    let output = stdout(interpret(&setup));
    assert!(output.contains("product_discriminant=-36, tate_rank=4, locally_saturated=false"));
    assert!(output.contains("product_discriminant=-16, tate_rank=4, locally_saturated=true"));
    assert!(output.contains("rvKummerCMCyclePair = product=10, normalized_kummer=5"));
    assert!(output.contains("rvKummerCMCycleHermitian = a=10, d=2, b_in_1_omega=[-1, 1]"));
    let error = diagnostic(interpret(&format!(
        "{setup}\nshow rvKummerCMAlgebraic[U](w);"
    )));
    assert!(error.contains("Gram is not a local unit"), "{error}");
}

#[test]
fn rejects_bad_geometric_witnesses_and_nominal_operands() {
    for (command, expected) in [
        ("rvkummercmcycles Bad[P,B](rvUnramElem[K]{6});", "require an elliptic basis"),
        ("rvkummercmcycles Bad[P,EB](rvUnramElem[K]{7});", "base-field square root of -1"),
        ("rvkummercmcycles Bad[P,EB](rvUnramElem[L]{6});", "Kummer CM Gaussian unit"),
        ("rvkummercmcycles Bad[P,EB,0](rvUnramElem[K]{6});", "index must be positive"),
        ("rvKummerCMDefect[A,0];", "degree must be positive"),
        ("show rvKummerCMCycleLift[3](graph);", "between one and the marking height"),
        ("show rvKummerCMCycle[A]{1,2,3};", "requires four coefficients"),
        ("rvkummercmcycles W[P,EB](rvUnramElem[K]{31}); rvKummerCMCyclePair(graph,rvKummerCMCycle[W]{0,0,0,1});", "Kummer CM intersection argument"),
        ("rvkummerbrauerpairing R[B,S,2]; show rvKummerCMTranscendental[A](zero[RVKummerBrauerFixed(R,1)]);", "Kummer CM projection argument"),
    ] {
        let error = diagnostic(interpret(&format!("{EXAMPLE}\n{command}")));
        assert!(error.contains(expected), "{error}");
    }
}
