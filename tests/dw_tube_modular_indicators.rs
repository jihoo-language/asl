mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-tube-modular-indicators.asl");

fn type_iii() -> &'static str {
    EXAMPLE.split_once("// The same bulk").unwrap().0
}

fn added_exterior_cocycle() -> &'static str {
    let start = EXAMPLE.find("let exterior=").unwrap();
    EXAMPLE[start..]
        .split_once("rvDoubleTubeRotationData(rvDoubleTubeRotate[V,4,2]);")
        .unwrap()
        .0
}

fn cyclic(twisted: bool) -> String {
    let cocycle = if twisted {
        "rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2}"
    } else {
        "rvPostnikovK[P]"
    };
    format!(
        "rvgroup G[3];rvdwcoeff P[G,3];rvdouble D[P]({cocycle});
        rvdoubleboundary Z[D]{{0}};rvdoubledefects C[Z,Z];rvdoubletube X[C];"
    )
}

#[test]
fn example_recovers_modular_data_and_distinguishes_the_two_plane_categories() {
    let result = stdout(example("twisted-double-tube-modular-indicators"));
    assert_eq!(result.matches("center=22, conductor=4, exponent=4, orientation=conjugated_double_crossing, verified=true").count(),2);
    assert!(result.contains("simple=4, dual=4, power=2, dimension=1, trace=[1]/1, verified=true"));
    assert!(result.contains("simple=4, dual=4, power=2, dimension=1, trace=[-1]/1, verified=true"));
    assert_eq!(
        result
            .matches("simple=4, dual=4, power=4, dimension=4, trace=[2]/1, verified=true")
            .count(),
        2
    );
    assert!(result.ends_with("RVCycloQ(9){[0, 0, 0, 1] / 1}\nRVCycloQ(9){[1, 0, 0, 2] / 3}\n"));
}

#[test]
fn all_first_eight_indicators_agree_with_the_independent_type_iii_formula() {
    let mut source = format!(
        "{}rvdoubletube T[PP];{}",
        type_iii(),
        added_exterior_cocycle()
    );
    source.push_str(
        "let one=rvCycloRat[1](rvCyclo[4,0]);let two=add(one,one);let z=zero[RVCycloQ(4)];",
    );
    for (tube, sign) in [("T", "one"), ("V", "neg(one)")] {
        for i in 0..5 {
            for n in 1..=8 {
                let expected = if i == 0 || (i < 4 && n % 2 == 0) {
                    "one"
                } else if i < 4 || n % 2 == 1 {
                    "z"
                } else if n % 4 == 0 {
                    "two"
                } else {
                    sign
                };
                writeln!(
                    source,
                    "assert rvDoubleTubeIndicator[{tube},{i},{n}]=={expected};"
                )
                .unwrap();
            }
        }
        writeln!(
            source,
            "rvDoubleTubeIndicatorPeriod[{tube},0];rvDoubleTubeIndicatorPeriod[{tube},4];"
        )
        .unwrap();
    }
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvDoubleTubeIndicatorPeriod[0] = period=1, exponent=4\n",
            "rvDoubleTubeIndicatorPeriod[4] = period=4, exponent=4\n",
            "rvDoubleTubeIndicatorPeriod[0] = period=1, exponent=4\n",
            "rvDoubleTubeIndicatorPeriod[4] = period=4, exponent=4\n",
        )
    );
}

#[test]
fn pivotal_rotation_uses_nontrivial_matrix_entries_and_checks_the_odd_empty_space() {
    let mut source = format!(
        "{}rvdoubletube T[PP];{}
        let t=rvDoubleTubeRotate[T,4,4];let v=rvDoubleTubeRotate[V,4,4];
        let one=rvCycloRat[1](rvCyclo[4,0]);let z=zero[RVCycloQ(4)];
        assert rvDoubleTubeRotationTrace(t)==add(one,one);
        rvDoubleTubeRotationData(rvDoubleTubeRotate[T,4,3]);",
        type_iii(),
        added_exterior_cocycle()
    );
    for (c, r0) in [0, 3, 2, 1].into_iter().enumerate() {
        for r in 0..4 {
            let expected = if r == r0 { "one" } else { "z" };
            for m in ["t", "v"] {
                writeln!(
                    source,
                    "assert rvDoubleTubeRotationEntry[{r},{c}]({m})=={expected};"
                )
                .unwrap();
            }
        }
    }
    assert_eq!(stdout(interpret(&source)),
        "rvDoubleTubeRotationData = simple=4, dual=4, power=3, dimension=0, trace=[]/1, verified=true\n");
}

#[test]
fn center_words_act_in_idempotent_coordinates_with_the_signed_complex_orientation() {
    let mut source = cyclic(true);
    source.push_str("
        let one=rvCycloRat[1](rvCyclo[9,0]);let phase=rvCycloRat[1](rvCyclo[9,1]);
        let e=rvDoubleTubeCentral[X,3];let e7=rvDoubleTubeCentral[X,7];
        let s=rvDoubleTubeWord{S}(e);
        assert rvDoubleTubeSEntry[X,3,3]==rvCycloRat[3](rvCyclo[9,7]);
        assert rvDoubleTubeBlockEntry[3,0,0](s)==rvDoubleTubeSEntry[X,3,3];
        assert rvDoubleTubeWord{S,S}(e)==e7;
        assert rvDoubleTubeWord{T[1],S,T[1],S,T[1],S}(e)==e7;
        assert rvDoubleTubeWord{T[1]}(e)==rvDoubleTubeScale(phase,e);
        assert rvDoubleTubeWord{T[1],T[-1]}(e)==e;
        assert rvDoubleTubeWord{}(e)==e;
        assert rvDoubleTubeTwist[X,3]==phase;
        assert rvDoubleTubeLens[X,3]==rvDoubleLens[D,3];
        assert rvDoubleTubeLens[X,-3]==rvCycloQGalois[-1](rvDoubleTubeLens[X,3]);
        assert rvDoubleTubeLens[X,3]==rvDoubleTubeBlockEntry[0,0,0](rvDoubleTubeWord{S,T[3],S}(rvDoubleTubeCentral[X,0]));
        assert rvDoubleTubeIndicator[X,1,3]==rvCycloRat[1](rvCyclo[9,3]);
        let r=rvDoubleTubeRotate[X,1,3];
        assert rvDoubleTubeRotationEntry[0,0](r)==rvCycloRat[1](rvCyclo[9,3]);");
    for a in 0..9 {
        writeln!(source, "assert rvDoubleTubeTwist[X,{a}]==rvDoubleT[D,{a}];").unwrap();
        for b in 0..9 {
            writeln!(
                source,
                "assert rvDoubleTubeSEntry[X,{a},{b}]==rvDoubleS[D,{a},{b}];"
            )
            .unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn cyclic_untwisted_and_coboundary_models_keep_their_exact_conductors() {
    let untwisted = format!(
        "{}let r=rvDoubleTubeRotate[X,1,3];
        assert rvDoubleTubeRotationTrace(r)==rvCycloRat[1](rvCyclo[9,0]);
        rvDoubleTubeIndicatorPeriod[X,1];",
        cyclic(false)
    );
    assert_eq!(
        stdout(interpret(&untwisted)),
        "rvDoubleTubeIndicatorPeriod[1] = period=3, exponent=3\n"
    );
    let b = |i| i32::from(i == 1);
    let mu = (1..4)
        .flat_map(|x| (1..4).map(move |y| format!("({},3)", b(x) + b(y) - b(x ^ y))))
        .collect::<Vec<_>>()
        .join(",");
    let coboundary = format!(
        "rvgroup G{{{{0,1,2,3}},{{1,0,3,2}},{{2,3,0,1}},{{3,2,1,0}}}};
        rvdwcoeff P[G,2];rvdouble D[P](rvPostnikovK[P]);rvdoubleboundary H[D]{{0,1,2,3}}{{{mu}}};
        rvdoubledefects T[H,H];rvdoubletube Y[T];rvDoubleTubeModular[Y];
        let r=rvDoubleTubeRotate[Y,1,2];
        assert rvDoubleTubeRotationTrace(r)==rvCycloRat[1](rvCyclo[12,0]);"
    );
    assert_eq!(stdout(interpret(&coboundary)),
        "rvDoubleTubeModular[Y] = center=16, conductor=12, exponent=2, orientation=conjugated_double_crossing, verified=true\n");
}

#[test]
fn modular_and_direct_rotation_queries_enforce_their_domains() {
    let setup = format!("{}rvdoubletube T[PP];", type_iii());
    for (op, message) in [
        (
            "show rvDoubleTubeWord{S}(rvDoubleTubeMatrixUnit[T,8,0,1]);",
            "requires an element of the tube center",
        ),
        (
            "show rvDoubleTubeRotate[T,4,5];",
            "supports indicator degrees 2 through 4",
        ),
        (
            "show rvDoubleTubeIndicator[T,4,0];",
            "indicator index must be positive",
        ),
        (
            "rvDoubleTubeIndicatorPeriod[T,5];",
            "outside the boundary category",
        ),
        (
            "show rvDoubleTubeRotationEntry[4,0](rvDoubleTubeRotate[T,4,4]);",
            "matrix index is out of range",
        ),
        (
            "rvDoubleTubeRotationData(rvDoubleTubeUnit[T]);",
            "requires a pivotal tube rotation",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{op}")));
        assert!(error.contains(message), "{message}: {error}");
    }
}
