mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn signs_convolve_to_sl2_with_actual_inverse_convolutions() {
    assert_eq!(
        stdout(example("sheaf-middle-convolution")),
        concat!(
            "shdiskrigidity[Seed] = 2\n",
            "shconvrank[MC] = 2\n",
            "shdiskmat[Left, T] = [[1, 3], [0, 1]]\n",
            "shdiskmat[Right, T] = [[1, 0], [2, 1]]\n",
            "shdiskinfinity[MC] = [[1, 2], [3, 2]]\n",
            "shdiskgroup[MC] = 120\n",
            "shdiskrigidity[MC] = 2\n",
            "shconvcoords = [[3], [4]]\n",
            "shconvcoords = [1, 1]\n",
            "shconvrank[Back] = 1\n",
            "shdiskmonosame[Seed, Back] = true\n",
            "shconvrank[Twist] = 2\n",
            "shconvrank[Untwist] = 1\n",
            "shdiskmonosame[Seed, Untwist] = true\n",
        )
    );
}

#[test]
fn both_relation_spaces_and_functoriality_under_rectangular_maps() {
    let program = r"
        shdisk Sign[5,1,1] [1] [2];
        shdisk A[5,2,2] [1,0,0,1] [0,0,0,2];
        shdisk B[5,2,2] [1,0,0,1] [-2,0,0,2];
        shdiskseq S[Sign,Sign];
        shdiskseq Double[A,B];
        shconv Small[S,2];
        shconv Big[Double,2];
        shconvrank[Small];
        shconvrank[Big];

        let k=shConvRaw[Big]{[1,0],[0,0]};
        let l=shConvRaw[Big]{[3,0],[1,0]};
        assert shConvProject(k) == zero[ShConvVec(Big)];
        assert shConvProject(l) == zero[ShConvVec(Big)];
        assert shConvRawT[1](k) == pscale[2](k);
        assert shConvRawT[2](k) == k;
        assert shConvRawT[1](l) == l;
        assert shConvRawT[2](l) == l;
        let b=shConvRaw[Big]{[4,-2],[1,2]};
        assert shConvRawT[-1](shConvRawT[1](b)) == b;
        assert shConvRawT[2](shConvRawT[-2](b)) == b;

        shdiskmap iA[Sign,A] [0,1] [0,1];
        shdiskmap iB[Sign,B] [0,1] [0,1];
        shdiskmap pA[A,Sign] [0,1] [0,1];
        shdiskmap pB[B,Sign] [0,1] [0,1];
        shdiskseqmap f[S,Double]{iA,iB};
        shdiskseqmap g[Double,S]{pA,pB};
        shdiskseqmap comp[S,S]{shDiskComp(pA,iA),shDiskComp(pB,iB)};
        shdiskseqmap id[S,S]{shDiskId[Sign],shDiskId[Sign]};
        shconvmap cf[Small,Big,f];
        shconvmap cg[Big,Small,g];
        shconvmap cc[Small,Small,comp];
        shconvmap ci[Small,Small,id];
        shdiskseqmapat[cf,0];
        shdiskseqmapat[cg,1];

        let x = add(shConvGen[Small,0],pscale[3](shConvGen[Small,1]));
        let raw = shConvSection(x);
        assert shConvMap[ci](x)==x;
        assert shConvMap[cg](shConvMap[cf](x)) == shConvMap[cc](x);
        assert shConvMap[cc](x)==x;
        assert shConvMap[cg](shConvMap[cf](raw))==raw;
        assert shConvProject(shConvMap[cf](raw))==shConvMap[cf](x);
        assert shConvT[1](shConvMap[cf](x))==shConvMap[cf](shConvT[1](x));
        assert shConvT[-2](shConvMap[cf](x))==shConvMap[cf](shConvT[-2](x));
        assert shConvRawT[-1](shConvMap[cf](raw))==shConvMap[cf](shConvRawT[-1](raw));
        assert shConvRawT[2](shConvMap[cf](raw))==shConvMap[cf](shConvRawT[2](raw));
        card (\q:ShConvVec(Big) => S.hi);

        shdiskseqmap double[S,S]{pscale[2](shDiskId[Sign]),pscale[2](shDiskId[Sign])};
        shdiskseqmap quad[S,S]{pscale[4](shDiskId[Sign]),pscale[4](shDiskId[Sign])};
        shconvmap d[Small,Small,double];
        shconvmap q[Small,Small,quad];
        assert shConvMap[d](shConvMap[d](x)) == shConvMap[q](x);
        assert shConvMap[q](x)==pscale[4](x);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "shconvrank[Small] = 2\n",
            "shconvrank[Big] = 2\n",
            "shdiskseqmapat[cf, 0] = [[1, 0], [0, 1]] / [[1]]\n",
            "shdiskseqmapat[cg, 1] = [[1, 0], [0, 1]] / [[1]]\n",
            "card = 25 / 25\n",
        )
    );
}

#[test]
fn trivial_actions_and_zero_dimensional_quotients() {
    let program = r"
        shdisk Zero[3,0,0] [] [];
        shdisk I[3,1,0] [] [];
        shdiskseq S[Zero,Zero];
        shdiskseq T[I,I,I];
        shconv Empty[S,2];
        shconv Fixed[T,-1];
        let a=shConvRaw[Empty]{[],[]};
        let b=shConvRaw[Fixed]{[1],[2],[-1]};
        assert shConvProject(a) == zero[ShConvVec(Empty)];
        assert shConvProject(b) == zero[ShConvVec(Fixed)];
        assert shConvRawT[-2](shConvRawT[2](b)) == b;
        assert shConvRawT[3](b)==shConvRaw[Fixed]{[1],[2],[1]};
        shconvcoords(a);
        shconvcoords(shConvProject(b));
        shconvrank[Empty];
        shconvrank[Fixed];
        shdiskinfinity[Empty];
        shdiskgroup[Fixed];
        shdiskrigidity[Fixed];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "shconvcoords = [[], []]\n",
            "shconvcoords = []\n",
            "shconvrank[Empty] = 0\n",
            "shconvrank[Fixed] = 0\n",
            "shdiskinfinity[Empty] = []\n",
            "shdiskgroup[Fixed] = 1\n",
            "shdiskrigidity[Fixed] = 0\n",
        )
    );
}

#[test]
fn rejects_wrong_parameters_dimensions_and_morphism_grades() {
    let setup = r"
        shdisk P[5,1,1] [1] [2];
        shdiskseq S[P,P];
        shconv C[S,2];
        shconv D[S,-1];
        shdiskseqmap f[S,S]{shDiskId[P],shDiskId[P]};
        shconvmap cf[C,C,f];
    ";
    let cases = [
        ("shconv Bad[S,0];", "different from zero and one"),
        ("shconv Bad[S,6];", "different from zero and one"),
        ("shconvmap Bad[C,D,f];", "equal field scalars"),
        ("shconvmap Bad[C,C,cf];", "given original endpoints"),
        ("show shConvRaw[C]{[1]};", "requires 2 blocks"),
        ("show shConvRaw[C]{[1,0],[0]};", "requires 1 coordinates"),
        (
            "show shConvRawT[-3](zero[ShConvRaw(C)]);",
            "generator -3 is out of range",
        ),
        (
            "show shConvT[0](zero[ShConvVec(C)]);",
            "generator 0 is out of range",
        ),
        (
            "show shConvProject(zero[ShConvVec(C)]);",
            "requires ShConvRaw",
        ),
        (
            "show shConvMap[cf](zero[ShConvVec(D)]);",
            "declared source quotient",
        ),
        (
            "shconvcoords(S.hi);",
            "expected a middle convolution vector",
        ),
    ];
    for (tail, expected) in cases {
        let error = diagnostic(interpret(&format!("{setup}{tail}")));
        assert!(error.contains(expected), "{error}");
    }
}
