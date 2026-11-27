mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/sheaf-disk.asl");

#[test]
fn distinct_vanishing_cycles_dual_maps_and_actual_relative_lifts() {
    let result = stdout(example("sheaf-disk"));
    let expected = [
        "shdiskmat[P, T] = [[1, 0], [0, 1]]",
        "shdiskmat[Q, T] = [[1, 0], [0, 1]]",
        "shdiskmat[P, S] = [[1, 4], [0, 1]]",
        "shdiskmat[Q, S] = [[1, 0], [0, 1]]",
        "shdiskranks[P] = nearby 2, vanishing 2, stalk [1, 1], costalk [1, 1]",
        "shdiskranks[Q] = nearby 2, vanishing 2, stalk [1, 1], costalk [1, 1]",
        "shdiskhomdim[P, Q] = 3",
        "shdiskmapmat = [[0, 1], [0, 2]] / [[0, 0], [3, 1]]",
        "shdiskmat[Left, T] = [[1, 1], [0, 1]]",
        "shdiskmat[Right, T] = [[1, 0], [4, 1]]",
        "shdiskranks[Mid] = nearby 2, vanishing 1, stalk [1, 0], costalk [0, 1]",
        "shlocaldim[AtTip, pt, relative, 1] = 1",
        "shlocaldim[AtTip, pt, relative, 2] = 1",
        "shlocaldim[AtLeft, pt, relative, 1] = 0",
        "shlocaldim[AtLeft, pt, relative, 2] = 1",
        "shhyperdim[Stalk, 0] = 1",
        "shhyperdim[Stalk, 1] = 1",
        "shhyperdim[Costalk, 1] = 1",
        "shhyperdim[Costalk, 2] = 1",
    ];
    assert_eq!(result.lines().collect::<Vec<_>>(), expected);
}

#[test]
fn characteristic_two_zero_width_and_finite_hom_values() {
    let program = r"
        shdisk A[2,1,1] [1] [0];
        shdisk B[2,0,1] [] [];
        shdisk C[2,1,0] [] [];
        shdiskranks[B];
        shdiskranks[C];
        shdiskdual DA[A];
        shdiskdual DDA[DA];
        assert shDiskComp(shDiskBidualInv[A,DDA],shDiskBidual[A,DDA])==shDiskId[A];
        assert shDiskGen[A,A,0]==shDiskId[A];
        assert pscale[2](shDiskId[A])==zero[ShDiskHom(A,A)];
        let identity = \f:ShDiskHom(A,A) => f;
        let compose = \f:ShDiskHom(A,A) => shDiskComp(f,shDiskId[A]);
        assert identity==compose;
        card (\f:ShDiskHom(A,A) => S.hi);
        shdiskhomdim[B,C];
        shdiskmiddle M[B];
        shdiskranks[M];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "shdiskranks[B] = nearby 0, vanishing 1, stalk [0, 1], costalk [1, 0]
",
            "shdiskranks[C] = nearby 1, vanishing 0, stalk [1, 0], costalk [0, 1]
",
            "card = 2 / 2
",
            "shdiskhomdim[B, C] = 0
",
            "shdiskranks[M] = nearby 0, vanishing 0, stalk [0, 0], costalk [0, 0]
",
        )
    );
}

#[test]
fn actual_constraints_reject_singular_objects_and_wrong_maps() {
    let setup = EXAMPLE.split("// The labels").next().unwrap();
    let errors = [
        ("shdisk Bad[5,1,1] [1] [1];", "nearby monodromy I - v u is singular"),
        ("shdisk Bad[5,1,1] [] [1];", "requires 1 rows and 1 columns"),
        ("shdisk Bad[4,0,0] [] [];", "prime"),
        ("shdiskmap Bad[P,Q] [1,0,0,1] [1,0,0,1];", "intertwine both"),
        ("shdiskmap Bad[P,P] [0,0,0,0] [1,0,0,0];", "intertwine both"),
        ("show shDiskCanonical[P,Q];", "requires a shriek source or a star target"),
        ("show shDiskCanonical[Lower,Q];", "same based nearby monodromy"),
        ("show shDiskBidual[P,DQ];", "requires the declared double dual"),
        ("show shDiskDualMap[DP,DQ](f);", "target dual followed by the source dual"),
        ("show shDiskExtMap[Lower,Upper](twice);", "same kind of extension"),
        ("show shDiskComp(f,f);", "matching middle quivers"),
        ("shdisk Z[2,0,0] [] []; shdiskhomdim[P,Z];", "same prime field"),
        ("space Wrong {a<b<c<d}; let j=tab Wrong -> Wrong {a=>Wrong.a,b=>Wrong.b,c=>Wrong.c,d=>Wrong.d}; shdiskrealize E[P,j,a,a,b,c,d];", "four links"),
    ];
    for (tail, expected) in errors {
        let error = diagnostic(interpret(&format!("{setup}{tail}")));
        assert!(error.contains(expected), "{error}");
    }
}
