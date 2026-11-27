mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn fork_kan_extensions_and_geometric_base_change_counterexamples() {
    assert_eq!(
        stdout(example("set-kan")),
        concat!(
            "{lo: [0, 1], hi: [1, 0, 0]}\n",
            "{lo: [0, 1], hi: [1, 0]}\n",
            "{pt: []}\n",
            "{pt: [0, 1]}\n",
            "{pt: [0, 0]}\n",
            "{pt: [0]}\n",
        )
    );
}

#[test]
fn left_and_right_adjunctions_agree_on_entire_hom_spaces() {
    let result = interpret(
        r"
        space X {b,c,a,a<b,a<c};
        setfamily F[X] {a:2,b:2,c:3} {a<b:[0,0],a<c:[1,0]};
        let f = tab X -> S {a=>S.lo,b=>S.lo,c=>S.hi};
        setlan L[F](f); setran R[F](f);
        setfamily H[S] {lo:2,hi:2} {lo<hi:[1,1]};
        setpullback PH[H](f);
        let l = \v : SetNat(L,H) => setLanLift[L](setLanRestrict[PH](v));
        let lp = \v : SetNat(F,PH) => setLanRestrict[PH](setLanLift[L](v));
        let r = \v : SetNat(H,R) => setRanTranspose[R](setRanUntranspose[PH](v));
        let rp = \v : SetNat(PH,F) => setRanUntranspose[PH](setRanTranspose[R](v));
        assert l == \v : SetNat(L,H) => v;
        assert lp == \v : SetNat(F,PH) => v;
        assert r == \v : SetNat(H,R) => v;
        assert rp == \v : SetNat(PH,F) => v;
        card \v : SetNat(L,H) => S.hi;
        card \v : SetNat(F,PH) => S.hi;
        card \v : SetNat(H,R) => S.hi;
        card \v : SetNat(PH,F) => S.hi;
    ",
    );
    let output = stdout(result);
    let counts: Vec<_> = output.lines().collect();
    assert_eq!(counts.len(), 4, "{output}");
    assert_eq!(counts[0], counts[1], "{output}");
    assert_eq!(counts[2], counts[3], "{output}");
    assert!(
        counts.iter().all(|line| !line.starts_with("card = 0")),
        "{output}"
    );
}

#[test]
fn family_projection_agrees_with_dependent_sum_and_product() {
    let result = interpret(
        r"
        setfamily G[S] {lo:2,hi:1} {lo<hi:[0,0]};
        elements E[G];
        setfamily F[E] {lo_0:2,lo_1:2,hi_0:2}
            {lo_0<hi_0:[0,1],lo_1<hi_0:[1,0]};
        let p = setProjection[E];
        setsum Sigma[E,F]; setproduct Pi[E,F];
        setlan L[F](p); setran R[F](p);
        setpullback PS[Sigma](p); setpullback PP[Pi](p);
        setpullback PL[L](p); setpullback PR[R](p);
        let l = setLanLift[L](setSumRestrict[PS](setId[Sigma]));
        let lback = setSumLift[Sigma](setLanUnit[L,PL]);
        let r = setRanTranspose[R](setEval[Pi,PP]);
        let rback = setTranspose[Pi](setRanCounit[R,PR]);
        assert setComp(l,lback) == setId[Sigma];
        assert setComp(lback,l) == setId[L];
        assert setComp(r,rback) == setId[R];
        assert setComp(rback,r) == setId[Pi];
    ",
    );
    assert_eq!(stdout(result), "");
}

#[test]
fn base_change_along_a_constant_family_projection_is_an_isomorphism() {
    let result = interpret(
        r"
        setfamily G[S] {lo:2,hi:2} {lo<hi:[0,1]};
        elements E[G];
        setfamily F[E] {lo_0:2,lo_1:3,hi_0:1,hi_1:2}
            {lo_0<hi_0:[0,0],lo_1<hi_1:[1,0,1]};
        let p = setProjection[E];
        setlan L[F](p); setran R[F](p);
        space Slice {u,v};
        let g = tab One -> S {pt=>S.lo};
        let gp = tab Slice -> E {u=>E.lo_1,v=>E.lo_0};
        let pp = tab Slice -> One {u=>One.pt,v=>One.pt};
        setpullback FP[F](gp);
        setlan LP[FP](pp); setran RP[FP](pp);
        setpullback PL[L](g); setpullback PR[R](g);
        setnat lback[PL,LP] {pt:[3,4,0,1,2]};
        setnat rback[RP,PR] {pt:[0,3,1,4,2,5]};
        let l = setLanBC[LP,PL];
        let r = setRanBC[RP,PR];
        assert setComp(l,lback) == setId[PL];
        assert setComp(lback,l) == setId[LP];
        assert setComp(r,rback) == setId[RP];
        assert setComp(rback,r) == setId[PR];
    ",
    );
    assert_eq!(stdout(result), "");
}

#[test]
fn a_commuting_square_can_repeat_a_section_coordinate() {
    let result = interpret(
        r"
        setfamily F[One] {pt:2} {};
        let id = \x : One => x;
        setlan L[F](id); setran R[F](id);
        setpullback LP[L](id); setpullback RP[R](id);
        space Double {x,y};
        let p = \x : Double => One.pt;
        setpullback FP[F](p);
        setlan LL[FP](p); setran RR[FP](p);
        setnat folded[LL,LP] {pt:[0,1,0,1]};
        setnat repeated[RP,RR] {pt:[0,3]};
        assert setLanBC[LL,LP] == folded;
        assert setRanBC[RR,RP] == repeated;
    ",
    );
    assert_eq!(stdout(result), "");
}

#[test]
fn pushout_fibers_and_a_missing_required_choice() {
    let result = interpret(
        r"
        space X {root<left,root<right};
        setfamily F[X] {root:1,left:2,right:3} {root<left:[1],root<right:[1]};
        let p = \x : X => One.pt;
        setlan L[F](p); setran R[F](p);
        setsections[L]; setsections[R];
        setfamily G[X] {root:0,left:2,right:2} {root<left:[],root<right:[]};
        setlan GL[G](p); setran GR[G](p);
        setsections[GL]; setsections[GR];
    ",
    );
    assert_eq!(
        stdout(result),
        "setsections[L] = 4\nsetsections[R] = 1\nsetsections[GL] = 4\nsetsections[GR] = 0\n"
    );
}

#[test]
fn noncommuting_squares_wrong_provenance_and_large_sections_are_diagnosed() {
    let setup = r"
        setfamily F[One] {pt:2} {};
        let lo = tab One -> S {pt=>S.lo};
        let hi = tab One -> S {pt=>S.hi};
        setlan L[F](lo); setran R[F](lo);
        setpullback Bad[L](hi);
        setpullback Y[F](\x : One => x);
        setlan LY[Y](lo); setran RY[Y](lo);
        setpullback B[L](\s : S => S.hi);
        setpullback C[R](\s : S => S.hi);
    ";
    for (tail, expected) in [
        (
            "show setLanUnit[L,Bad];",
            "pulled back along the Kan extension's map",
        ),
        ("show setLanCounit[L];", "not a set pullback"),
        ("show setRanUnit[R];", "not a set pullback"),
        ("show setLanBC[LY,B];", "commuting square of monotone maps"),
        ("show setRanBC[RY,C];", "commuting square of monotone maps"),
        (
            "setlan Nope[F](S.lo);",
            "requires a monotone map of named sites",
        ),
        (
            "setran Nope[F](\\s : S => s);",
            "input family's site as domain",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{tail}")));
        assert!(error.contains(expected), "{tail}\n{error}");
    }
    let overflow = diagnostic(interpret(
        r"
        space Disc {a,b};
        setfamily Many[Disc] {a:65,b:65} {};
        setran Nope[Many](\x : Disc => One.pt);
    ",
    ));
    assert!(overflow.contains("limit of 4096 sections"), "{overflow}");
}
