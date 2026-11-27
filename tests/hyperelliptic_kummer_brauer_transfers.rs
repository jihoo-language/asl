mod common;

use common::{diagnostic, example, interpret, stdout};

const F37: &str = include_str!("../examples/hyperelliptic-kummer-brauer-transfers.asl");

const F17: &str = r"
rvunram K[17,1]{0,1};
rvhyperelliptic E[K]{rvUnramElem[K]{1},rvUnramElem[K]{0},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvellglue C[E,E](rvUnramElem[K]{0},rvUnramElem[K]{1},rvUnramElem[K]{1},rvUnramElem[K]{0});
rvjacchart S[C](rvUnramElem[K]{4});
rvunram H[17,1]{14,0,1}; rvjacbase EL[E,H];
rvellglue CL[EL,EL](rvUnramElem[H]{0},rvUnramElem[H]{1},rvUnramElem[H]{1},rvUnramElem[H]{0});
rvjacchart SL[CL](rvUnramElem[H]{4});
let p=rvJacPoint[EL](rvUnramElem[H]{1},rvUnramElem[H]{6});
let q=rvJacPoint[EL](rvUnramElem[H]{5},rvUnramElem[H]{0,5});
let e1=rvEllGlueMap[SL](p,zero[RVJac(EL)]);
let e2=rvEllGlueMap[SL](zero[RVJac(EL)],p);
let f1=rvEllGlueMap[SL](q,zero[RVJac(EL)]);
let f2=rvEllGlueMap[SL](zero[RVJac(EL)],q);
rvjacpowerbasis B[SL,3,2](rvUnramElem[H]{7,13}){e1,e2,add(f1,e1),f2};
rvkummerbrauerpairing P[B,S,1];
rvkummerbrauerpairing P2[B,S,2];
";

#[test]
fn certified_classes_and_actual_extension_adjunction() {
    let output = stdout(example("hyperelliptic-kummer-brauer-transfers"));
    for summary in [
        "rvKummerBrauerQuotientData[Q] = level=3, divisor_rank=4, reference_degree=1, invariant_factors=[3, 3]",
        "rvKummerBrauerQuotientData[Q3] = level=3, divisor_rank=4, reference_degree=1, invariant_factors=[3, 3]",
    ] {
        assert!(output.contains(summary), "{output}");
    }
    let links: Vec<_> = output
        .lines()
        .filter(|line| line.starts_with("rvKummerBrauerLink = "))
        .collect();
    let nonzero = [
        "rvKummerBrauerLink = numerator=2, order=3, nonalgebraic_witness=true",
        "rvKummerBrauerLink = numerator=1, order=3, nonalgebraic_witness=true",
    ];
    assert_eq!(
        links,
        [
            nonzero[0],
            nonzero[1],
            nonzero[1],
            nonzero[1],
            nonzero[1],
            "rvKummerBrauerLink = numerator=0, order=1, nonalgebraic_witness=false"
        ]
    );
}

#[test]
fn the_canonical_generators_and_lifts_cover_the_entire_finite_module() {
    let output = stdout(interpret(&format!(
        r"{F37}
        let g=rvKummerBrauerGenerator[Q,0]; let h=rvKummerBrauerGenerator[Q,1];
        assert rvKummerBrauerClass[Q](rvKummerBrauerLift(g))==g;
        assert rvKummerBrauerClass[Q](rvKummerBrauerLift(h))==h;
        assert rvKummerBrauerLinkNonzero(rvKummerBrauerLift(g),rvKummerBrauerLift(h))==S.hi;
        assert (\a:RVKummerBrauerClass(Q)=>rvKummerBrauerClass[Q](rvKummerBrauerLift(a)))
            ==(\a:RVKummerBrauerClass(Q)=>a);
        assert (\a:RVKummerBrauerClass(Q)=>rvKummerBrauerCorestrict[Q](rvKummerBrauerRestrict[Q2](a)))
            ==(\a:RVKummerBrauerClass(Q)=>rvKummerBrauerScale[2](a));
        assert (\a:RVKummerBrauerClass(Q2)=>rvKummerBrauerRestrict[Q2](rvKummerBrauerScale[2](rvKummerBrauerCorestrict[Q](a))))
            ==(\a:RVKummerBrauerClass(Q2)=>a);
        assert (\a:RVKummerBrauerClass(Q3)=>rvKummerBrauerCorestrict[Q](a))
            ==(\a:RVKummerBrauerClass(Q3)=>zero[RVKummerBrauerClass(Q)]);"
    )));
    assert!(output.contains("invariant_factors=[3, 3]"), "{output}");
}

#[test]
fn a_noncoordinate_integral_fixed_lattice_gives_the_actual_norm() {
    let output = stdout(interpret(&format!(
        "{F17}
         rvkummerbrauerquotient Q[P,2];
         rvKummerBrauerQuotientData[Q];
         let fixed=rvKummerBrauerFixed[P,2]{{1,0,0,1,0,0}};
         assert rvKummerBrauerAlgebraic[Q](fixed)==S.hi;
         assert rvKummerBrauerClass[Q](fixed)==zero[RVKummerBrauerClass(Q)];
         let upstairs=rvKummerBrauerFixed[P2,2]{{2,0,0,1,0,0}};
         let norm=rvKummerBrauerCorestrict[P](upstairs);
         assert norm==rvKummerBrauerScale[2](fixed);
         assert rvKummerBrauerRestrict[P2](norm)==rvKummerBrauerFixed[P2,2]{{2,0,0,2,0,0}};"
    )));
    assert!(output.contains("level=9, divisor_rank=4, reference_degree=1, invariant_factors=[]"));
    assert!(output.contains("[1, 0, 0, 1, 0, 0]"), "{output}");
}

#[test]
fn precision_and_nominal_failures_do_not_erase_lift_transfers() {
    for (query, expected) in [
        (
            "rvkummerbrauerquotient Bad[P,2];",
            "needs more marking height",
        ),
        (
            "rvkummerbrauerquotient Bad[P3,1];",
            "needs more marking height",
        ),
        (
            "rvkummerbrauerquotient Bad[P2,1,P3];",
            "same base and marking with nested degrees",
        ),
        (
            "show rvKummerBrauerClass[Q](w2);",
            "Kummer Brauer quotient argument",
        ),
        (
            "show rvKummerBrauerRestrict[P3](x);",
            "unknown Kummer Brauer quotient",
        ),
        (
            "show rvKummerBrauerCorestrict[Q2](x);",
            "same base and marking with nested degrees",
        ),
        (
            "show rvKummerBrauerGenerator[Q,2];",
            "generator index is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F37}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
    for (query, expected) in [
        ("rvkummerbrauerquotient Bad[P2,1,P];", "different abelian fixed rank"),
        ("rvkummerbrauerquotient Q[P,1]; rvkummerbrauerquotient R[P2,2]; show rvKummerBrauerRestrict[R](zero[RVKummerBrauerClass(Q)]);", "requires the same torsion exponent"),
        ("show rvKummerBrauerFixed[P,1]{1,0,0,0,0,0};", "not fixed at the specified exponent"),
    ] {
        let error = diagnostic(interpret(&format!("{F17}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}
