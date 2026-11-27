mod common;

use common::{diagnostic, example, interpret, stdout};

const ARTIN: &str = include_str!("../examples/set-artin-gluing.asl");
const FANO: &str = include_str!("../examples/set-fano-gluing.asl");

#[test]
fn relabeled_reconstruction_and_non_discrete_open_limit() {
    assert_eq!(
        stdout(example("set-artin-gluing")),
        concat!(
            "{dc: [0, 2, 2], cc: [1, 0]}\n",
            "{v: [0, 0, 1], c: [0, 0], d: [0, 1, 1], u: [0, 0]}\n",
        )
    );
    let tail = r"
        setnat bad[B,B2] {vu:[1,0,1],uu:[0,0]};
        show setGlueMap[G,G2](a,bad);
    ";
    let error = diagnostic(interpret(&format!("{ARTIN}{tail}")));
    assert!(
        error.contains("attachment square fails at closed point 'dc' and element 0"),
        "{error}"
    );
}

#[test]
fn fano_same_restrictions_separate_attachments_and_isomorphic_reflections() {
    let output = stdout(example("set-fano-gluing"));
    let mut lines = output.lines();
    let points = lines.next().unwrap();
    let lines_alpha = lines.next().unwrap();
    let closed_x = ", cx: [0, 57, 114, 171, 228, 285, 342]}";
    assert!(
        points.ends_with(closed_x) && lines_alpha.ends_with(closed_x),
        "{output}"
    );
    assert_ne!(points, lines_alpha);
    assert_eq!(
        lines.collect::<Vec<_>>(),
        [
            "card = 343 / 343",
            "card = 343 / 343",
            "setmonostabilizer[M, 0] = [24, 1, 0]",
            "setmonomaps[M] = 0",
            "card = 0 / 0",
            "card = 1 / 1",
        ]
    );

    let error = diagnostic(interpret(&format!("{FANO}show setGlueMap[GP,GL](a,b);")));
    assert!(
        error.contains("attachment square fails at closed point 'cy'"),
        "{error}"
    );
}

#[test]
fn empty_strata_and_terminal_empty_section() {
    let output = stdout(interpret(
        r"
        space Empty {};
        let z = tab Empty -> One {};
        let one = tab One -> One {pt=>One.pt};
        setfamily F[One] {pt:3} {};
        setrecollement AllClosed(one,z);
        setrestrict A[AllClosed,F,closed];
        setrestrict B[AllClosed,F,open];
        setgluelimit T[AllClosed,B];
        show setAttach[A,T];
        setglue G[AllClosed](setAttach[A,T]);
        assert setComp(setGlueRecover[G,F],setGlueForget[G,F]) == setId[G];
        assert setGlueClosed[G,G](setGlueMap[G,G](setId[A],setId[B])) == setId[A];
        setrecollement AllOpen(z,one);
        setrestrict AC[AllOpen,F,closed];
        setrestrict AO[AllOpen,F,open];
        setgluelimit E[AllOpen,AO];
        show setAttach[AC,E];
        setglue GG[AllOpen](setAttach[AC,E]);
        assert setComp(setGlueForget[GG,F],setGlueRecover[GG,F]) == setId[F];
        assert setGlueOpen[GG,GG](setGlueMap[GG,GG](setId[AC],setId[AO])) == setId[AO];
        let e = tab Empty -> Empty {};
        setfamily EmptyF[Empty] {} {};
        setrecollement Both(e,e);
        setrestrict EC[Both,EmptyF,closed];
        setrestrict EO[Both,EmptyF,open];
        setgluelimit ET[Both,EO];
        setglue EB[Both](setAttach[EC,ET]);
        assert setComp(setGlueForget[EB,EmptyF],setGlueRecover[EB,EmptyF]) == setId[EmptyF];
        let low = tab One -> S {pt=>S.lo};
        let high = tab One -> S {pt=>S.hi};
        setrecollement Chain(low,high);
        setfamily NoValues[S] {lo:0,hi:0} {lo<hi:[]};
        setrestrict Lo[Chain,NoValues,closed];
        setrestrict Hi[Chain,NoValues,open];
        setgluelimit NoSections[Chain,Hi];
        card (\m:SetNat(F,NoSections) => S.hi);
        setglue NoValuesAgain[Chain](setAttach[Lo,NoSections]);
        assert setComp(setGlueForget[NoValuesAgain,NoValues],setGlueRecover[NoValuesAgain,NoValues]) == setId[NoValues];

    ",
    ));
    assert_eq!(output, "{pt: [0, 0, 0]}\n{}\ncard = 0 / 0\n");
}

#[test]
fn wrong_embeddings_origins_and_attachments_report_specific_errors() {
    let site = r"
        space X {a<b<c};
        space Discrete {h,k};
        space Empty {};
        let nothing = tab Empty -> X {};
        let a = tab One -> X {pt=>X.a};
        let c = tab One -> X {pt=>X.c};
        let ends = tab Discrete -> X {h=>X.a,k=>X.b};
        let reverse = tab Discrete -> X {h=>X.b,k=>X.c};
        let id = tab One -> One {pt=>One.pt};
        let left = tab S -> X {lo=>X.a,hi=>X.b};
        let right = tab One -> X {pt=>X.c};
    ";
    for (command, expected) in [
        ("setrecollement Bad(ends,c);", "must reflect the order"),
        ("setrecollement Bad(a,a);", "not injective and disjoint"),
        ("setrecollement Bad(a,c);", "omits ambient point 'b'"),
        ("setrecollement Bad(a,reverse);", "must reflect the order"),
        (
            "setrecollement Bad(right,left);",
            "open image is not upward closed",
        ),
        (
            "setrecollement Bad(nothing,id);",
            "must share a named codomain",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{site}{command}")));
        assert!(error.contains(expected), "{error}");
    }
    for (tail,expected) in [
        ("show setAttach[A,T2];", "is not a set pullback"),
        ("show setGlueRecover[G2,F];", "is not a set pullback"),
        ("setglue Bad[R](a);", "is not a set gluing limit"),
        ("show setGlueMap[G,G2](setId[A],b);", "expected SetNat(A, A2)"),
        ("setrecollement R(closed,open);", "already defined"),
        ("let differentC=tab C -> X {cc=>X.c,dc=>X.u}; let differentU=tab U -> X {uu=>X.d,vu=>X.v}; setrecollement Other(differentC,differentU); setglue Bad[Other](alpha);", "must target a limit of this recollement"),
        ("setnat wrong[A,T] {cc:[0,0],dc:[0,2,2]}; setglue Bad[R](wrong); show setGlueRecover[Bad,F];", "requires the original canonical attachment"),
    ] {
        let error = diagnostic(interpret(&format!("{ARTIN}{tail}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn matching_objects_enforce_the_existing_section_storage_limit() {
    let error = diagnostic(interpret(
        r"
        space X {c<a,c<b};
        space U {u,v};
        let c = tab One -> X {pt=>X.c};
        let open = tab U -> X {u=>X.a,v=>X.b};
        setrecollement R(c,open);
        setfamily B[U] {u:65,v:65} {};
        setgluelimit T[R,B];
    ",
    ));
    assert!(error.contains("limit of 4096 sections"), "{error}");
}
