mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn later_truth_is_dense_and_excluded_middle_fails() {
    assert_eq!(
        stdout(example("set-logic")),
        "setsieve[O, lo, 1] = {hi}\n{lo: [1], hi: [1]}\n{lo: [2], hi: [1]}\n",
    );
}

#[test]
fn comprehension_reindexes_images_and_preserves_transports() {
    let output = stdout(interpret(
        r"
        setomega O[S];
        setfamily F[S] {lo:3,hi:3} {lo<hi:[1,2,0]};
        setfamily A[S] {lo:2,hi:3} {lo<hi:[0,2]};
        setnat mono[A,F] {lo:[2,0],hi:[0,2,1]};
        let p = setClassify[O](mono);
        setsub B(p);
        let incl = setInclusion[B];
        assert setClassify[O](incl) == p;
        setnat pto[A,B] {lo:[1,0],hi:[0,2,1]};
        setnat back[B,A] {lo:[1,0],hi:[0,2,1]};
        assert setComp(incl,pto) == mono;
        assert setComp(back,pto) == setId[A];
        assert setComp(pto,back) == setId[B];
        show p;
        show incl;
    ",
    ));
    assert_eq!(
        output,
        "{lo: [2, 1, 2], hi: [1, 1, 1]}\n{lo: [0, 2], hi: [0, 1, 2]}\n"
    );
}

#[test]
fn implication_and_quantifier_adjunctions_use_future_transport() {
    let output = stdout(interpret(
        r"
        setomega O[S];
        setfamily A[S] {lo:1,hi:3} {lo<hi:[0]};
        setfamily B[S] {lo:2,hi:2} {lo<hi:[1,0]};
        setfamily Good[S] {lo:1,hi:2} {lo<hi:[0]};
        setnat p[A,B] {lo:[0],hi:[1,1,0]};
        setnat include[Good,A] {lo:[0],hi:[0,2]};
        let P = setClassify[O](include);
        let exists = setExists(p,P);
        let forall = setForall(p,P);
        setnat e[B,O] {lo:[2,1],hi:[1,1]};
        setnat a[B,O] {lo:[0,2],hi:[1,0]};
        assert e == exists;
        assert a == forall;
        let pulled_e = setComp(exists,p);
        let pulled_a = setComp(forall,p);
        assert setAnd(P,pulled_e) == P;
        assert setAnd(P,pulled_a) == pulled_a;
        assert setAnd(exists,setExists(p,pulled_e)) == setExists(p,pulled_e);
        assert setAnd(forall,setForall(p,pulled_a)) == forall;
        assert setImplies(exists,forall) == forall;
        assert setForall(p,setTop[A,O]) == setTop[B,O];
        assert setExists(p,setBottom[A,O]) == setBottom[B,O];
        show exists;
        show forall;
    ",
    ));
    assert_eq!(
        output,
        "{lo: [2, 1], hi: [1, 1]}\n{lo: [0, 2], hi: [1, 0]}\n"
    );
}

#[test]
fn empty_fibers_sites_and_universal_vacuity() {
    assert_eq!(
        stdout(interpret(
            r"
            space Void {};
            setomega W[Void];
            setfamily F[Void] {} {};
            assert setTop[F,W] == setBottom[F,W];
            setsub C(setBottom[F,W]);
            assert setClassify[W](setInclusion[C]) == setTop[F,W];

            setomega O[S];
            setfamily Empty[S] {lo:0,hi:0} {lo<hi:[]};
            setfamily T[S] {lo:1,hi:1} {lo<hi:[0]};
            setnat p[Empty,T] {lo:[],hi:[]};
            let v = setTop[Empty,O];
            assert setExists(p,v) == setBottom[T,O];
            assert setForall(p,v) == setTop[T,O];
        "
        )),
        "",
    );
}

#[test]
fn checks_monic_components_classifier_identity_and_finite_bounds() {
    let mismatch = r"
        setomega O[S];
        setomega W[S];
        setfamily A[S] {lo:2,hi:2} {lo<hi:[0,1]};
        setfamily B[S] {lo:1,hi:1} {lo<hi:[0]};
        setnat p[A,B] {lo:[0,0],hi:[0,0]};
    ";
    for (script, expected) in [
        ("show setClassify[O](p);", "injective component"),
        ("show setAnd(setTop[A,O],setTop[A,W]);", "same classifier"),
        ("show setForall(p,setTop[B,O]);", "source family"),
        ("setsieve[O,lo,3];", "outside its fiber"),
        ("setsub B2(p);", "declared set truth classifier"),
    ] {
        let error = diagnostic(interpret(&format!("{mismatch}{script}")));
        assert!(error.contains(expected), "{error}");
    }
    let branches = (0..13)
        .map(|i| format!("x<y{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let error = diagnostic(interpret(&format!(
        "space X {{{branches}}}; setomega O[X];"
    )));
    assert!(error.contains("limit of 4096 cosieves"), "{error}");
}

#[test]
fn cosieves_transport_across_machine_word_boundaries() {
    let points = (0..70)
        .map(|i| format!("x{i}"))
        .collect::<Vec<_>>()
        .join("<");
    let source = format!("space X {{{points}}}; setomega O[X]; setsieve[O,x0,1]; setsections[O];");
    assert_eq!(
        stdout(interpret(&source)),
        "setsieve[O, x0, 1] = {x69}\nsetsections[O] = 71\n",
    );
}
