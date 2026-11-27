mod common;

use common::{diagnostic, example, interpret, stdout};

const SEED: &str = r"
    shdisk P[5,1,1] [1] [2];
    shdiskseq S[P,P,P,P];
    shhurwitz H[S] [1];
";

#[test]
fn fixed_local_system_has_noncommuting_parabolic_braid_generators() {
    let result = stdout(example("sheaf-parabolic-hurwitz"));
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(lines[0], "shdiskprod[Sphere] = [[1]]");
    assert_eq!(lines[1], "shparrank[Sphere] = 2");
    assert_eq!(lines[2], "shdiskmonosame[Sphere, H1] = true");
    assert_eq!(lines[3], "shdiskseqsame[Sphere, H1] = false");
    assert_ne!(lines[4], lines[5]);
    assert_eq!(lines[6], "shdiskseqmapat[fix, 1] = [[1]] / [[4]]");
    assert_eq!(lines[7], "shdiskseqsame[A, B] = true");
    assert_eq!(lines.len(), 8);
}

#[test]
fn rank_two_mixed_quivers_transport_raw_classes_and_coupled_witnesses() {
    let result = stdout(interpret(
        r"
        shdisk P[5,2,1] [0,1] [-1,0];
        shdisk Q[5,2,2] [1,0,0,0] [0,0,1,0];
        shdisk PI[5,2,1] [0,1] [1,0];
        shdisk QI[5,2,2] [1,0,0,0] [0,0,-1,0];
        shdisk N[5,2,2] [1,0,0,1] [2,0,0,2];
        shdiskseq S[P,Q,QI,PI,N,N];
        shdiskseq T[P,Q,QI,PI,N,N];
        shdiskmap dp[P,P] [2,0,0,2] [2];
        shdiskmap dq[Q,Q] [2,0,0,2] [2,0,0,3];
        shdiskmap dpi[PI,PI] [2,0,0,2] [2];
        shdiskmap dqi[QI,QI] [2,0,0,2] [2,0,0,4];
        shdiskmap dn[N,N] [2,0,0,2] [2,0,0,2];
        shdiskseqmap twice[S,T]{dp,dq,dqi,dpi,dn,dn};
        shhurwitz HS[S] [1,-2,4,5,-3,2,1,-5];
        shhurwitz HT[T] [1,-2,4,5,-3,2,1,-5];
        shhurwitz Back[HS] [5,-1,-2,3,-5,-4,2,-1];
        shdiskprod[HS];
        shdiskseqsame[Back,S];
        shparrank[S];
        shparrank[HS];
        shdiskseqmapbraid braided[HS,HT,twice];
        shdiskseqmapat[braided,0];
        let a=shParRaw[S]{[0,0],[0,0],[0,0],[0,0],[1,0],[1,0]};
        let h=shParClass(a);
        let w=shParVec[S][1,-1];
        let b=shParBoundary(w);
        assert shParPrimitive(b)==w;
        assert shParPrimitive(shParBraid[HS](b))==shParBraid[HS](w);
        assert shParBraid[HT](shParMap[twice](a))==shParMap[braided](shParBraid[HS](a));
        assert shParBraid[HT](shParMap[twice](h))==shParMap[braided](shParBraid[HS](h));
        assert shParSum(shParBraid[HS](a))==zero[ShParVec(HS)];
        assert shParIdentify[S](shParBraid[Back](shParBraid[HS](a)))==a;
        assert shParIdentify[S](shParBraid[Back](shParBraid[HS](h)))==h;
        shparrows(shParBraid[HS](a));
        shparvector(shParPrimitive(b));
        ",
    ));
    assert_eq!(
        result,
        concat!(
            "shdiskprod[HS] = [[1, 0], [0, 1]]\n",
            "shdiskseqsame[Back, S] = true\n",
            "shparrank[S] = 4\n",
            "shparrank[HS] = 4\n",
            "shdiskseqmapat[braided, 0] = [[2, 0], [0, 2]] / [[2, 0], [0, 2]]\n",
            "shparrows = [[0, 4], [0, 0], [4, 3], [0, 0], [0, 0], [1, 0]]\n",
            "shparvector = [1, 4]\n",
        )
    );
}

#[test]
fn affine_tuples_and_degenerate_parabolic_spaces() {
    let result = stdout(interpret(
        r"
        shdisk E[2,1,1] [1] [0];
        shdisk Z[2,0,1] [] [];
        shdiskseq S[E,E,E,E];
        shdiskseq EmptyNearby[Z,Z,Z];
        shhurwitz H[S] [-1,2,1];
        shparrank[S];
        shparrank[H];
        shparrank[EmptyNearby];
        assert zero[ShParRaw(S)] == shParRaw[S]{[2],[0],[0],[0]};
        assert zero[ShParVec(EmptyNearby)] == shParVec[EmptyNearby][];
        assert shParPrimitive(shParBoundary(shParVec[S][1])) == zero[ShParVec(S)];
        card (\x:ShParVec(S) => S.hi);
        card (\x:ShPar(S) => S.hi);
        shdisk E5[5,1,1] [1] [2];
        shdisk Q5[5,1,0] [] [];
        shdiskseq Affine[E5,Q5];
        shhurwitz Moved[Affine] [1];
        shdiskprod[Moved];
        ",
    ));
    assert_eq!(
        result,
        concat!(
            "shparrank[S] = 0\n",
            "shparrank[H] = 0\n",
            "shparrank[EmptyNearby] = 0\n",
            "card = 2 / 2\n",
            "card = 1 / 1\n",
            "shdiskprod[Moved] = [[4]]\n",
        )
    );
}

#[test]
fn rejects_invalid_braids_raw_cochains_and_a_false_nearby_map() {
    let cases = [
        ("shdiskseq Empty[];", "at least one quiver"),
        ("shhurwitz Bad[S] [4];", "generator 4 is out of range"),
        ("shhurwitz Bad[S] [0];", "generator 0 is out of range"),
        ("shdiskpoint Bad[H,4];", "strand index out of range"),
        ("shdiskseq N[P]; shparrank[N];", "identity total monodromy"),
        ("show shParRaw[S]{[0]};", "requires 4 strands"),
        ("show shParVec[S][0,1];", "requires 1 coordinates"),
        (
            "show shParClass(shParRaw[S]{[1],[0],[0],[0]});",
            "transported sum is nonzero",
        ),
        (
            "show shParPrimitive(shParRaw[S]{[1],[0],[-1],[0]});",
            "not a coupled coboundary",
        ),
        (
            "show shParBraid[S](shParGen[S,0]);",
            "declared Hurwitz image",
        ),
        (
            "show shParBraid[H](shParGen[H,0]);",
            "declared source sequence",
        ),
        (
            "shdisk F[5,1,1] [0] [0]; shdiskseq T[F,F,F,F]; show shParIdentify[T](shParGen[S,0]);",
            "equal ordered based monodromies",
        ),
        (
            "shdiskseqmap f[S,S]{shDiskId[P],zero[ShDiskHom(P,P)],shDiskId[P],shDiskId[P]};",
            "differ on the nearby map",
        ),
        (
            "shdisk F[5,1,1] [0] [0]; shdiskseq T[F,F]; show shParRaw[T]{[1],[0]};",
            "outside im(I - T_i)",
        ),
    ];
    for (tail, expected) in cases {
        let error = diagnostic(interpret(&format!("{SEED}{tail}")));
        assert!(error.contains(expected), "{error}");
    }
}
