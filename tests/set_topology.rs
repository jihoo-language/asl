mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn reflection_over_non_discrete_support_uses_compatible_sections_and_the_adjunction() {
    assert_eq!(
        stdout(example("set-topologies")),
        concat!(
            "settopologysupport[Recovered] = {a, t, b}\n",
            "setreflectcheck[J, F] = no\n",
            "setreflectcheck[J, R] = yes\n",
            "setreflectcheck[J, G] = yes\n",
            "setsections[F] = 1\n",
            "setsections[R] = 2\n",
            "card = 8 / 8\n",
            "card = 4 / 4\n",
            "card = 4 / 4\n",
            "settopologyjoinsteps[E, U] = 4\n",
        )
    );
}

const TINY: &str = r"
    space X {a<b};
    setomega O[X];
    settopology A[O] {a};
    settopology B[O] {b};
    setfamily T[X] {a:1,b:1} {a<b:[0]};
    setnat half[T,O] {a:[1],b:[1]};
    let ab = setComp(setTopMap[A],setTopMap[B]);
    let ba = setComp(setTopMap[B],setTopMap[A]);
    let empty = setBottom[T,O];
    assert setComp(ab,empty) == half;
    assert setComp(ab,setComp(ab,empty)) == setTop[T,O];
    assert setComp(ba,empty) == setTop[T,O];
";

#[test]
fn supplied_maps_reconstruct_the_subset_and_compositions_need_checking() {
    let source = format!(
        r"{TINY}
        setnat declared[O,O] {{a:[1,1,2],b:[1,1]}};
        settopologyMap Recovered(declared);
        assert setTopMap[Recovered] == setTopMap[A];
        settopologysupport[Recovered];
        settopologyMap Constant(ba);
        settopologysupport[Constant];
        settopologyjoinsteps[A,B];
        settopologyjoinsteps[B,A];
        "
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "settopologysupport[Recovered] = {a}\n",
            "settopologysupport[Constant] = {}\n",
            "settopologyjoinsteps[A, B] = 2\n",
            "settopologyjoinsteps[B, A] = 1\n",
        )
    );
    for (statement, expected) in [
        ("settopologyMap Bad(ab);", "candidate is not idempotent"),
        (
            "settopologyMap Bad(setBottom[O,O]);",
            "does not preserve true",
        ),
        (
            "settopologyMap Bad(setId[T]);",
            "not a declared set truth classifier",
        ),
        ("settopologyMap Bad(empty);", "requires an endomorphism"),
        ("settopology Bad[O] {missing};", "unknown point"),
        (
            "settopology Bad[T] {a};",
            "not a declared set truth classifier",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{TINY} {statement}")));
        assert!(error.contains(expected), "{statement}: {error}");
    }
}

#[test]
fn extremal_topologies_and_double_negation_agree_with_their_specialized_models() {
    let source = r"
        space X {r<a,r<b};
        setomega O[X];
        setomega Alias[X];
        settopology E[O] {};
        settopology Full[Alias] {r,a,b};
        settopology Max[O] {a,b};
        assert setTopMap[E] == setTop[O,O];
        assert setTopMap[Full] == setId[Alias];
        assert setTopMap[Max] == setNot(setNot(setId[O]));
        settopologymeet M[Max,Full];
        assert setTopMap[M] == setId[O];
        settopologyjoin P[Full,Max];
        assert setTopMap[P] == setNot(setNot(setId[Alias]));
        setfamily F[X] {r:0,a:2,b:3} {r<a:[],r<b:[]};
        setreflect One[E,F];
        setreflect I[Full,F];
        setreflect D[Max,F];
        setj Old[F];
        assert setComp(setReflectInverse[I],setReflectUnit[I]) == setId[F];
        assert setComp(setReflectUnit[I],setReflectInverse[I]) == setId[I];
        let newer = setJExtend[Old](setReflectUnit[D]);
        let older = setReflectExtend[D](setJUnit[Old]);
        assert setComp(newer,older) == setId[D];
        assert setComp(older,newer) == setId[Old];
        setsections[One];
        setsections[D];
        setreflectcheck[Full,F];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "setsections[One] = 1\nsetsections[D] = 6\nsetreflectcheck[Full, F] = yes\n"
    );
}

#[test]
fn the_join_checks_the_actual_number_of_alternating_passes() {
    for (size, expected) in [(30, 16), (89, 45)] {
        let points: Vec<_> = (0..size).map(|i| format!("x{i}")).collect();
        let even = points
            .iter()
            .step_by(2)
            .cloned()
            .collect::<Vec<_>>()
            .join(",");
        let odd = points
            .iter()
            .skip(1)
            .step_by(2)
            .cloned()
            .collect::<Vec<_>>()
            .join(",");
        let chain = points.join("<");
        let source = format!(
            r"space X {{{chain}}};
            setomega O[X];
            settopology E[O] {{{even}}};
            settopology U[O] {{{odd}}};
            settopologyjoin J[E,U];
            assert setTopMap[J] == setTop[O,O];
            settopologyjoinsteps[E,U];
            "
        );
        assert_eq!(
            stdout(interpret(&source)),
            format!("settopologyjoinsteps[E, U] = {expected}\n")
        );
    }
}

#[test]
fn extension_rejects_an_unreflected_target_and_unrelated_site() {
    let program = include_str!("../examples/set-topologies.asl");
    for (extra, expected) in [
        (
            "show setReflectExtend[R](setId[F]);",
            "target must be a sheaf",
        ),
        (
            "show setReflectExtend[R](setId[R]);",
            "map from the original family",
        ),
        (
            "show setReflectInverse[R];",
            "source of reflection is not a sheaf",
        ),
        (
            "show setReflectUnit[F];",
            "not a declared topology reflection",
        ),
        ("setreflect Elsewhere[J,P];", "topology's named site"),
        ("settopologyjoin Bad[J,E];", "same named site"),
    ] {
        let error = diagnostic(interpret(&format!("{program} {extra}")));
        assert!(error.contains(expected), "{extra}: {error}");
    }
}

#[test]
fn the_adjunction_preserves_an_empty_hom_set() {
    let source = r"
        space X {a};
        setomega O[X];
        settopology Full[O] {a};
        setfamily F[X] {a:1} {};
        setfamily G[X] {a:0} {};
        setreflect R[Full,F];
        let extend = \u : SetNat(F,G) => setReflectExtend[R](u);
        card (\u : SetNat(F,G) => S.hi);
        card (\u : SetNat(R,G) => S.hi);
    ";
    assert_eq!(stdout(interpret(source)), "card = 0 / 0\ncard = 0 / 0\n");
}

#[test]
fn reflection_matches_right_kan_extension_of_restriction_along_the_full_subset() {
    let source = r"
        space X {r,a,b,t,r<a<t,r<b<t};
        space Kept {a,b,t,a<t,b<t};
        setomega O[X];
        settopology J[O] {a,b,t};
        setfamily F[X] {r:1,a:2,b:2,t:2}
          {r<a:[0],r<b:[1],a<t:[0,1],b<t:[1,0]};
        setreflect R[J,F];
        let i = tab Kept -> X {a=>X.a,b=>X.b,t=>X.t};
        setpullback RF[F](i);
        setran K[RF](i);
        setpullback RR[R](i);
        setnat forward[R,K] {r:[0,1],a:[0,1],b:[0,1],t:[0,1]};
        setnat backward[K,R] {r:[0,1],a:[0,1],b:[0,1],t:[0,1]};
        assert setReflectExtend[R](setRanUnit[K]) == forward;
        assert setComp(forward,setReflectUnit[R]) == setRanUnit[K];
        assert setComp(forward,backward) == setId[K];
        assert setComp(backward,forward) == setId[R];
        let u = setRanUntranspose[RR](forward);
        assert setRanTranspose[K](u) == forward;
        assert setComp(u,setPullNat[RF,RR](setReflectUnit[R])) == setId[RF];
        setreflectcheck[J,K];
    ";
    assert_eq!(stdout(interpret(source)), "setreflectcheck[J, K] = yes\n");
}
