mod common;

use common::{diagnostic, example, interpret, stdout};

const F37: &str = include_str!("../examples/hyperelliptic-brauer-twisted-mukai.asl");
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
rvkummerbrauerquotient Q[P,1];
rvkummerbrauerquotient Q2[P,2];
";

#[test]
fn actual_f37_mukai_frobenius_and_exact_period_after_cubic_extension() {
    let output = stdout(example("hyperelliptic-brauer-twisted-mukai"));
    let get = |name| {
        output
            .lines()
            .find(|row| row.starts_with(&format!("rvBrauerMukaiData[{name}, ")))
            .unwrap()
    };
    for (name, parts) in [
        (
            "M",
            &[
                "base_degree=1, precision=1, modulus=3, denominator=3",
                "cohomological_index=3, delta=[1, 0, 1, 2, 0, 0]",
                "frobenius_columns=[[1, 1, 0, 1, 2, 0, 0, 2]",
                "defect_power_ranks_mod_prime=[2, 1, 0]",
            ][..],
        ),
        (
            "M3",
            &[
                "base_degree=3",
                "cohomological_index=3, delta=[0, 0, 0, 0, 0, 0]",
                "defect_power_ranks_mod_prime=[0]",
            ][..],
        ),
        (
            "V",
            &[
                "numerator=[3, 1, 8, 1, 1, 1]",
                "cohomological_index=3, delta=[1, 0, 1, 2, 0, 0]",
            ][..],
        ),
        (
            "G",
            &[
                "cohomological_index=1, delta=[0, 0, 0, 0, 0, 0]",
                "defect_power_ranks_mod_prime=[0]",
            ][..],
        ),
    ] {
        let row = get(name);
        for part in parts {
            assert!(row.contains(part), "{row}");
        }
    }
    assert_eq!(
        output
            .lines()
            .filter(|row| *row == "rvBrauerMukaiPair = residue=2, modulus=3")
            .count(),
        2
    );
}

#[test]
fn signed_integral_gauge_intertwines_nonidentity_frobenius_and_keeps_pairing() {
    let output = stdout(interpret(&format!(
        "{F17}
         rvbrauermukai M[Q,1](zero[RVKummerBrauerFixed(P,1)]);
         rvbrauermukaishift N[M]{{1,0,0,0,0,0}};
         let x=rvBrauerMukai[M]{{1,0,0,0,0,0,1,0}};
         let y=rvBrauerMukai[M]{{0,2,0,0,0,0,0,0}};
         let n=rvBrauerMukaiTransport[N](x);
         let ny=rvBrauerMukaiTransport[N](y);
         // The Mukai row includes the mixed term with the Kummer half-form.
         assert n==rvBrauerMukai[N]{{1,2,0,0,0,0,1,2}};
         assert rvBrauerMukaiTransport[M](n)==x;
         assert rvBrauerMukaiFrob[1](n)==rvBrauerMukaiTransport[N](rvBrauerMukaiFrob[1](x));
         assert rvBrauerMukaiFrob[7](n)==rvBrauerMukaiTransport[N](rvBrauerMukaiFrob[7](x));
         assert rvBrauerMukaiFrob[3](add(n,ny))==add(rvBrauerMukaiFrob[3](n),rvBrauerMukaiFrob[3](ny));
         rvBrauerMukaiPair(x,y);
         rvBrauerMukaiPair(n,ny);
         rvBrauerMukaiData[N];
         let witness=rvBrauerMukaiWitness[N];
         assert rvBrauerMukaiFrob[1](witness)==witness;
         show witness;
         rvbrauermukai M9[Q,1](zero[RVKummerBrauerFixed(P,1)]);
         rvBrauerMukaiData[M9,0];"
    )));
    assert!(
        output.contains("cohomological_index=1, delta=[1, 0, 0, 0, 0, 0]"),
        "{output}"
    );
    assert!(
        output.contains("RVBrauerMukai(N){[1, 2, 0, 0, 0, 0, 0, 0]}"),
        "{output}"
    );
    let pairs: Vec<_> = output
        .lines()
        .filter(|line| line.starts_with("rvBrauerMukaiPair"))
        .collect();
    assert_eq!(pairs, ["rvBrauerMukaiPair = residue=2, modulus=3"; 2]);
}

#[test]
fn actual_change_of_primitive_weil_root_transports_twist_at_full_height() {
    let output = stdout(interpret(&format!(
        "{F37}
        rvjacpowerbasis BT[SL,3,2](rvUnramElem[L]{{16}}){{e1,e2,rvJacScale[2](add(f1,e1)),rvJacScale[2](f2)}};
        rvkummerbrauerpairing T[BT,S,1];
        rvkummerbrauerquotient QT[T,1];
        rvbrauermukaimarking MT[M,QT];
        let shifted=rvBrauerMukaiTransport[MT](fu);
        let ut=rvBrauerMukaiTransport[MT](u);
        assert rvBrauerMukaiFrob[1](ut)==shifted;
        assert rvBrauerMukaiTransport[M](shifted)==fu;
        rvBrauerMukaiPair(shifted,rvBrauerMukaiTransport[MT](t));
        rvBrauerMukaiData[MT];
        let v23=rvKummerBrauerFixed[P,1]{{0,0,0,0,0,1}};
        rvbrauermukai Wide[Q,1](v23);
        rvbrauermukaimarking WideTarget[Wide,QT];
        rvBrauerMukaiData[WideTarget];"
    )));
    assert!(output.contains("rvBrauerMukaiData[MT, 1]"));
    assert!(output.contains("numerator=[2, 0, 0, 0, 0, 0], cohomological_index=3"));
    let carried = output
        .lines()
        .find(|row| row.starts_with("rvBrauerMukaiData[WideTarget, 1]"))
        .unwrap();
    assert!(
        carried.contains("numerator=[0, 0, 8, 0, 0, 5], cohomological_index=3"),
        "{carried}"
    );
    assert_eq!(
        output
            .lines()
            .filter(|row| *row == "rvBrauerMukaiPair = residue=2, modulus=3")
            .count(),
        3
    );
}

#[test]
fn genuine_level_twenty_seven_supports_two_choices_of_denominator_and_precision() {
    // E(F_53^2)[27] has independent +1 and -1 Frobenius eigenvectors
    // (36,21) and (19,12*z), with z^2=2.
    let output = stdout(interpret(
        r"
rvunram K[53,1]{0,1};
rvhyperelliptic E[K]{rvUnramElem[K]{1},rvUnramElem[K]{0},rvUnramElem[K]{0},rvUnramElem[K]{1}};
rvellglue C[E,E](rvUnramElem[K]{0},rvUnramElem[K]{1},rvUnramElem[K]{1},rvUnramElem[K]{0});
rvjacchart S[C](rvUnramElem[K]{23});
rvunram H[53,1]{51,0,1}; rvjacbase EL[E,H];
rvellglue CL[EL,EL](rvUnramElem[H]{0},rvUnramElem[H]{1},rvUnramElem[H]{1},rvUnramElem[H]{0});
rvjacchart SL[CL](rvUnramElem[H]{23});
let p=rvJacPoint[EL](rvUnramElem[H]{36},rvUnramElem[H]{21});
let q=rvJacPoint[EL](rvUnramElem[H]{19},rvUnramElem[H]{0,12});
let e1=rvEllGlueMap[SL](p,zero[RVJac(EL)]);
let e2=rvEllGlueMap[SL](zero[RVJac(EL)],p);
let f1=rvEllGlueMap[SL](q,zero[RVJac(EL)]);
let f2=rvEllGlueMap[SL](zero[RVJac(EL)],q);
rvjacpowerbasis B[SL,3,3](rvUnramElem[H]{17,41}){e1,e2,f1,f2};
rvkummerbrauerpairing P[B,S,1];
rvkummerbrauerquotient Q[P,1];
rvbrauermukai M[Q,2](zero[RVKummerBrauerFixed(P,1)]);
rvbrauermukaishift N[M]{1,0,0,0,0,1};
let x=rvBrauerMukai[M]{1,0,0,0,0,0,0,0};
let nx=rvBrauerMukaiTransport[N](x);
assert nx==rvBrauerMukai[N]{1,8,0,0,0,0,8,4};
assert rvBrauerMukaiFrob[1](nx)==nx;
let y=rvBrauerMukai[N]{1,0,0,0,0,0,0,0};
assert rvBrauerMukaiFrob[1](y)==rvBrauerMukai[N]{1,7,0,0,0,0,7,7};
assert rvBrauerMukaiFrob[2](y)==y;
let witness=rvBrauerMukaiWitness[N];
assert rvBrauerMukaiFrob[1](witness)==witness;
rvBrauerMukaiData[N];

rvkummerbrauerquotient Higher[P,2];
rvbrauermukai W[Higher,1](rvKummerBrauerFixed[P,2]{0,4,0,0,0,0});
rvBrauerMukaiData[W];
",
    ));
    assert!(output.contains("precision=2, modulus=9, denominator=3"));
    assert!(output.contains("cohomological_index=1, delta=[7, 0, 0, 0, 0, 7]"));
    assert!(output.contains("precision=1, modulus=3, denominator=9"));
    assert!(output.contains("numerator=[0, 4, 0, 0, 0, 0], cohomological_index=1"));
}

#[test]
fn rejects_unjustified_precision_gauge_and_nominal_operations() {
    for (query, expected) in [
        (
            "rvbrauermukai Bad[Q,2](v);",
            "s+m at most the marking height",
        ),
        ("rvbrauermukai Bad[Q,0](v);", "positive precision"),
        (
            "rvbrauermukai Bad[Q3,1](v);",
            "Brauer Mukai B-field numerator",
        ),
        (
            "rvbrauermukaialgebraic Bad[M](v);",
            "an explicit Kummer CM cycle",
        ),
        (
            "rvbrauermukaiextend Bad[M3,Q];",
            "same base and marking with nested degrees",
        ),
        (
            "rvbrauermukaimarking Bad[M,Q3];",
            "same named base curve and extension degree",
        ),
        (
            "show rvBrauerMukaiTransport[M](u3);",
            "restriction has no inverse transport",
        ),
        (
            "show rvBrauerMukaiTransport[G](u);",
            "direct registered chart relation",
        ),
        ("rvBrauerMukaiPair(u,u3);", "Brauer Mukai pairing argument"),
        (
            "show rvBrauerMukai[M]{1,2,3};",
            "requires eight coefficients",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F37}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
    let error = diagnostic(interpret(&format!(
        "{F17}
         rvbrauermukai M[Q,1](zero[RVKummerBrauerFixed(P,1)]);
         rvbrauermukaiextend Bad[M,Q2];"
    )));
    assert!(
        error.contains("same torsion exponent and marking level"),
        "{error}"
    );
}
