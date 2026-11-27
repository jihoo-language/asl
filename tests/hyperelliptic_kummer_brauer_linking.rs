mod common;

use common::{diagnostic, example, interpret, stdout};

const MIXED: &str = r"
rvunram K[613,1]{0,1};
rvhyperelliptic E[K]{rvUnramElem[K]{0},rvUnramElem[K]{3},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvhyperelliptic F[K]{rvUnramElem[K]{0},rvUnramElem[K]{3},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvellglue C[E,F](rvUnramElem[K]{1},rvUnramElem[K]{131},rvUnramElem[K]{131},rvUnramElem[K]{1});
rvjacchart S[C](rvUnramElem[K]{241});
rvunram H[613,1]{611,0,0,1};
rvjacbase EL[E,H];
rvjacbase FL[F,H];
rvellglue CL[EL,FL](rvUnramElem[H]{1},rvUnramElem[H]{131},rvUnramElem[H]{131},rvUnramElem[H]{1});
rvjacchart SL[CL](rvUnramElem[H]{241});
let x=rvUnramElem[H]{64,113,266};
let y=rvUnramElem[H]{436,276,448};
let u=rvUnramElem[H]{573,165,217};
let v=rvUnramElem[H]{32,599,254};
let e1=rvEllGlueMap[SL](rvJacPoint[EL](x,y),zero[RVJac(FL)]);
let e2=rvEllGlueMap[SL](zero[RVJac(EL)],rvJacPoint[FL](x,y));
let f1=rvEllGlueMap[SL](rvJacPoint[EL](u,v),zero[RVJac(FL)]);
let f2=rvEllGlueMap[SL](zero[RVJac(EL)],rvJacPoint[FL](u,v));
rvjacpowerbasis B27[SL,3,3](rvUnramElem[H]{0,0,226}){e1,e2,f1,f2};
rvkummerbrauerpairing P27[B27,S,1];
let a1=rvKummerBrauerFixed[P27,1]{1,0,0,0,0,0};
let b2=rvKummerBrauerFixed[P27,2]{0,0,1,0,0,0};
rvKummerBrauerLink(a1,b2);
rvKummerBrauerLink(b2,a1);
rvKummerBrauerLink(a1,a1);
";

#[test]
fn independent_f37_linking_form_has_the_exact_sign_and_radical() {
    let output = stdout(example("hyperelliptic-kummer-brauer-linking"));
    let links: Vec<_> = output
        .lines()
        .filter(|line| line.starts_with("rvKummerBrauerLink = "))
        .collect();
    assert_eq!(
        links,
        [
            "rvKummerBrauerLink = numerator=2, order=3, nonalgebraic_witness=true",
            "rvKummerBrauerLink = numerator=1, order=3, nonalgebraic_witness=true",
            "rvKummerBrauerLink = numerator=0, order=1, nonalgebraic_witness=false",
            "rvKummerBrauerLink = numerator=1, order=3, nonalgebraic_witness=true",
            "rvKummerBrauerLink = numerator=2, order=3, nonalgebraic_witness=true",
            "rvKummerBrauerLink = numerator=0, order=1, nonalgebraic_witness=false",
        ]
    );
    let form = output
        .lines()
        .find(|line| line.starts_with("rvKummerBrauerLinkForm[P]"))
        .unwrap();
    for part in [
        "prime=3, fixed_dimension=6, rank=2, radical_dimension=4",
        "gram=[[0, 0, 2, 1, 0, 1], [0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1], [2, 0, 0, 0, 0, 2], [0, 0, 0, 0, 0, 0], [2, 0, 2, 1, 0, 0]]",
        "symplectic_pairs=[[[1, 0, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0]]]",
        "radical_basis=[[0, 1, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 1, 0], [2, 0, 1, 0, 0, 1]]",
    ] {
        assert!(form.contains(part), "{form}");
    }
    assert!(output.contains(
        "rvKummerBrauerLinkForm[P3] = prime=3, fixed_dimension=6, rank=0, radical_dimension=6"
    ));
    assert!(output.contains("brauer_invariant_factors=[9, 9], complete_primary=false"));
    assert!(output.contains("complete_with_artin_tate=true"));
}

#[test]
fn mixed_denominators_use_genuine_level_twenty_seven_and_nontrivial_twist() {
    // On both elliptic components, independent Frobenius columns are
    // (1,18),(18,19) mod 27, with multiplier 613=19 mod 27.
    // The glued Weil root is 226*z^2 for z^3=2.
    assert_eq!(
        stdout(interpret(MIXED)),
        concat!(
            "rvKummerBrauerLink = numerator=1, order=3, nonalgebraic_witness=true\n",
            "rvKummerBrauerLink = numerator=2, order=3, nonalgebraic_witness=true\n",
            "rvKummerBrauerLink = numerator=0, order=1, nonalgebraic_witness=false\n",
        )
    );
}

#[test]
fn fixed_basis_is_computed_on_a_proper_subspace() {
    let output = stdout(interpret(
        "rvunram K[17,1]{0,1};
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
         rvjacpowerbasis N[SL,3,2](rvUnramElem[H]{7,13}){e1,e2,f1,f2};
         rvkummerbrauerpairing L[N,S,1];
         rvKummerBrauerLinkForm[L];",
    ));
    // Frobenius on the two elliptic factors has eigenvalues 1,-1 mod 9.
    // The four mixed wedges are the fixed space after the Tate twist.
    assert!(output.contains("prime=3, fixed_dimension=4, rank=0, radical_dimension=4"));
    assert!(output.contains(
        "fixed_basis=[[0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0]]"
    ), "{output}");
}

#[test]
fn primitive_root_change_also_rescales_higher_precision_lifts() {
    let script = include_str!("../examples/hyperelliptic-kummer-brauer-linking.asl");
    let output = stdout(interpret(&format!(
        "{script}
        rvkummerbrauerpairing T3[Btw,S,3];
        let x=rvKummerBrauerFixed[P3,2]{{0,0,0,0,0,1}};
        let y=rvKummerBrauerFixed[T3,2]{{0,0,0,0,0,5}};
        assert rvKummerBrauerChange[T3](x)==y;
        assert rvKummerBrauerChange[P3](y)==x;
        assert rvKummerBrauerFixed[P3,2]{{0,0,0,0,0,-8}}==x;
        assert add(x,rvKummerBrauerScale[-1](x))==zero[RVKummerBrauerFixed(P3,2)];"
    )));
    assert!(output.contains("complete_with_artin_tate=true"));
}

#[test]
fn enforces_fixedness_nominality_and_pairing_precision() {
    let script = include_str!("../examples/hyperelliptic-kummer-brauer-linking.asl");
    for (query, expected) in [
        ("show rvKummerBrauerFixed[P,2]{1,0,0,0,0,0};", "not fixed at the specified exponent"),
        ("show rvKummerBrauerFixed[P,3]{0,0,0,0,0,0};", "between one and the marking height"),
        ("show rvKummerBrauerFixed[P,1]{0,0,0,0,0};", "requires six coordinates"),
        ("rvKummerBrauerLink(v,w3);", "requires the same named pairing"),
        ("show rvKummerBrauerChange[P3](v);", "requires the same named base curve and extension degree"),
        ("rvKummerBrauerLink(rvKummerBrauerFixed[P3,2]{1,0,0,0,0,0},w3);", "height at least the sum"),
        ("rvkummerbrauerpairing P0[B,S,0];", "extension degree must be positive"),
        ("rvjacpowerbasis Bbad[SL,37,1](rvUnramElem[L]{33}){e1,e2,f1,f2};", "differ from the residue characteristic"),
        ("rvjactorsionbasis B3[SL,3](rvUnramElem[L]{10}){rvJacScale[3](e1),rvJacScale[3](e2),rvJacScale[3](f1),rvJacScale[3](f2)};
          rvkummerbrauerpairing Q[B3,S,1]; rvKummerBrauerLinkForm[Q];", "requires marking height at least two"),
    ] {
        let error = diagnostic(interpret(&format!("{script}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn the_two_primary_kummer_lattice_is_excluded() {
    let script = include_str!("../examples/hyperelliptic-jacobian-level.asl");
    let error = diagnostic(interpret(&format!(
        "{script}\nrvkummerbrauerpairing P[B2,C,1];"
    )));
    assert!(error.contains("requires an odd torsion prime"), "{error}");
}
