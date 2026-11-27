mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn empty_current_fiber_still_has_future_choices() {
    assert_eq!(
        stdout(example("dependent-sets")),
        concat!(
            "{SumTotal.hi_0 => S.hi, SumTotal.hi_1 => S.hi}\n",
            "{PiTotal.lo_0 => S.lo, PiTotal.lo_1 => S.lo, PiTotal.hi_0 => S.hi, PiTotal.hi_1 => S.hi}\n",
            "{lo: [1, 0], hi: [0, 1]}\n",
        )
    );
}

const MERGING: &str = r"
    setfamily F[S] {lo:2,hi:1} {lo<hi:[0,0]};
    elements E[F];
    setfamily G[E] {lo_0:2,lo_1:2,hi_0:2}
        {lo_0<hi_0:[0,1],lo_1<hi_0:[1,0]};
    setsum Sum[E,G];
    setproduct Pi[E,G];
    let p = setProjection[E];
";

#[test]
fn merging_occurrences_force_future_agreement() {
    let program = r"
        setfamily P[S] {lo:2,hi:2} {lo<hi:[0,1]};
        setnat view[Pi,P] {lo:[0,1],hi:[0,1]};
        setnat inv[P,Pi] {lo:[0,1],hi:[0,1]};
        assert setComp(inv,view) == setId[Pi];
        assert setComp(view,inv) == setId[P];
        setfamily Sigma[S] {lo:4,hi:2} {lo<hi:[0,1,1,0]};
        setnat sumView[Sum,Sigma] {lo:[0,1,2,3],hi:[0,1]};
        setnat sumInv[Sigma,Sum] {lo:[0,1,2,3],hi:[0,1]};
        assert setComp(sumView,sumInv) == setId[Sigma];

        setpullback PP[Pi](p);
        setnat expected[PP,G] {lo_0:[0,1],lo_1:[1,0],hi_0:[0,1]};
        assert setEval[Pi,PP] == expected;
        assert setTranspose[Pi](expected) == setId[Pi];

        setpullback PS[Sum](p);
        setnat inclusion[G,PS] {lo_0:[0,1],lo_1:[2,3],hi_0:[0,1]};
        assert setSumRestrict[PS](setId[Sum]) == inclusion;
        assert setSumLift[Sum](inclusion) == setId[Sum];

        let count = \f : SetNat(P,Pi) => S.hi;
        let countPull = \f : SetNat(PP,G) => S.hi;
        card count;
        card countPull;
    ";
    let output = stdout(interpret(&format!("{MERGING}{program}")));
    let mut lines = output.lines();
    let first = lines.next().unwrap();
    let second = lines.next().unwrap();
    assert_eq!(first, second);
    assert!(first.contains('4'), "{output}");
}

#[test]
fn empty_future_is_a_singleton_and_an_empty_required_fiber_annihilates_sections() {
    let output = stdout(interpret(
        r"
        setfamily Void[S] {lo:0,hi:0} {lo<hi:[]};
        elements V[Void];
        setfamily Empty[V] {} {};
        setproduct P[V,Empty];
        setsum Z[V,Empty];
        setsections[P];
        setsections[Z];
        setfamily F[S] {lo:0,hi:1} {lo<hi:[]};
        elements E[F];
        setfamily G[E] {hi_0:0} {};
        setproduct Q[E,G];
        setsections[Q];
        setpullback PQ[Q](setProjection[E]);
        assert setTranspose[Q](setEval[Q,PQ]) == setId[Q];
    ",
    ));
    assert_eq!(
        output,
        "setsections[P] = 1\nsetsections[Z] = 0\nsetsections[Q] = 0\n",
    );
}

#[test]
fn provenance_checks_require_the_actual_map_and_matching_families() {
    let fixtures = [
        (
            "setproduct Bad[E,F];",
            "must live on the declared elements space",
        ),
        (
            "setpullback Wrong[F](\\x : One => x);",
            "family's site as codomain",
        ),
        (
            r"
                let collapse = \x : E => S.hi;
                setpullback Wrong[Pi](collapse);
                show setEval[Pi,Wrong];
            ",
            "pulled back along this elements projection",
        ),
        (
            r"
                let collapse = \x : E => S.hi;
                setpullback Wrong[Pi](collapse);
                setpullback Right[Pi](p);
                show setPullNat[Right,Wrong](setId[Pi]);
            ",
            "along the same map",
        ),
        (
            "setpullback Wrong[F](p); show setEval[Pi,Wrong];",
            "pullback of the indicated product",
        ),
    ];
    for (program, expected) in fixtures {
        let error = diagnostic(interpret(&format!("{MERGING}{program}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn large_unconstrained_future_sections_fail_with_a_bounded_diagnostic() {
    let error = diagnostic(interpret(
        r"
        setfamily F[One] {pt:13} {};
        elements E[F];
        setfamily G[E] {pt_0:2,pt_1:2,pt_2:2,pt_3:2,pt_4:2,pt_5:2,pt_6:2,
            pt_7:2,pt_8:2,pt_9:2,pt_10:2,pt_11:2,pt_12:2} {};
        setproduct P[E,G];
    ",
    ));
    assert!(error.contains("limit of 4096 sections"), "{error}");
}

#[test]
fn transitions_restrict_sections_to_two_different_future_branches() {
    let output = stdout(interpret(
        r"
        space Fork {x<y,x<z};
        setfamily F[Fork] {x:0,y:1,z:1} {x<y:[],x<z:[]};
        elements E[F];
        setfamily G[E] {y_0:2,z_0:3} {};
        setproduct P[E,G];
        setfamily R[Fork] {x:6,y:2,z:3}
            {x<y:[0,0,0,1,1,1],x<z:[0,1,2,0,1,2]};
        setnat view[P,R] {x:[0,1,2,3,4,5],y:[0,1],z:[0,1,2]};
        setnat inverse[R,P] {x:[0,1,2,3,4,5],y:[0,1],z:[0,1,2]};
        assert setComp(inverse,view) == setId[P];
        assert setComp(view,inverse) == setId[R];
        setpullback PR[R](setProjection[E]);
        setnat choose[PR,G] {y_0:[0,1],z_0:[0,1,2]};
        assert setTranspose[P](choose) == inverse;
        assert setUntranspose[PR](inverse) == choose;
        setsections[P];
    ",
    ));
    assert_eq!(output, "setsections[P] = 6\n");
}

#[test]
fn constraints_prune_a_large_family_of_candidates() {
    let occurrences = 48;
    let transitions = vec!["0"; occurrences].join(",");
    let sizes = (0..occurrences)
        .map(|i| format!("lo_{i}:2"))
        .collect::<Vec<_>>()
        .join(",");
    let links = (0..occurrences)
        .map(|i| format!("lo_{i}<hi_0:[0,1]"))
        .collect::<Vec<_>>()
        .join(",");
    let script = format!(
        "setfamily F[S] {{lo:{occurrences},hi:1}} {{lo<hi:[{transitions}]}};
         elements E[F];
         setfamily G[E] {{{sizes},hi_0:2}} {{{links}}};
         setproduct P[E,G];
         setsections[P];
         elements T[P];
         show setProjection[T];"
    );
    assert_eq!(
        stdout(interpret(&script)),
        concat!(
            "setsections[P] = 2\n",
            "{T.lo_0 => S.lo, T.lo_1 => S.lo, T.hi_0 => S.hi, T.hi_1 => S.hi}\n",
        )
    );
}
