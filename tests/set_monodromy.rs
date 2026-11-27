mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn fano_points_and_lines_have_equal_character_and_distinct_modular_linearizations() {
    let result = stdout(example("set-monodromy"));
    for line in [
        "setmonosize[M] = 168",
        "setmonochars[M] = {(0, 0): 48, (1, 1): 98, (3, 3): 21, (7, 7): 1}",
        "setmonostabilizer[M, 0] = [24, 1, 0]",
        "setmonomaps[M] = 0",
        "setmonoorbitals[M] = [28, 21]",
        "shExtDim[P5, L5, 0] = 2",
        "shmaprank[N, x] = 7",
        "shExtDim[P2, L2, 0] = 2",
        "shmaprank[Z, x] = 0",
        "shmaprank[J, x] = 1",
        "shmaprank[Off, x] = 3",
        "shmaprank[Inc, x] = 4",
    ] {
        assert!(result.lines().any(|output| output == line), "{result}");
    }
    assert_eq!(result.matches("[1, 8]").count(), 2, "{result}");
}

const DISCONNECTED: &str = r"
    space X {a<b,c<d};
    setfamily F[X] {a:2,b:2,c:2,d:1} {a<b:[1,0],c<d:[0,0]};
    setfamily G[X] {a:2,b:2,c:0,d:0} {a<b:[0,1],c<d:[]};
    setmonodromy M[F,G,a];
    shlinearize LF[F,5];
    shlinearize LG[G,5];
";

#[test]
fn tree_gauges_transport_signed_orbitals_on_just_the_selected_component() {
    let result = interpret(&format!(
        "{DISCONNECTED}
        setmonosize[M]; setmonoat[M,0]; setmonomaps[M]; setmonoorbitals[M];
        shorbitals A[M,LF,LG] {{0:1,0:-3,3:7}};
        shmap Exact[LF,LG] {{a:[3,0,0,2],b:[0,3,2,0],c:[],d:[]}};
        assert shMapC[A] == shMapC[Exact];
        shmaprank[A,a]; shmaprank[A,b]; shmaprank[A,c];
        shExtDim[LF,LG,0];"
    ));
    assert_eq!(
        stdout(result),
        concat!(
            "setmonosize[M] = 1\nsetmonoat[M, 0] = [0, 1] / [0, 1]\n",
            "setmonomaps[M] = 4\nsetmonoorbitals[M] = [1, 1, 1, 1]\n",
            "shmaprank[A, a] = 2\nshmaprank[A, b] = 2\nshmaprank[A, c] = 0\n",
            "shExtDim[LF, LG, 0] = 4\n",
        )
    );
}

#[test]
fn empty_fibers_have_the_expected_equivariant_maps_and_orbitals() {
    let result = interpret(
        r"
        space X {a<b};
        setfamily E[X] {a:0,b:0} {a<b:[]};
        setfamily F[X] {a:1,b:1} {a<b:[0]};
        setmonodromy A[E,F,a]; setmonodromy B[F,E,b];
        setmonomaps[A]; setmonomaps[B]; setmonoorbitals[A]; setmonochars[B];
        shlinearize LE[E,3]; shlinearize LF[F,3];
        shorbitals M[A,LE,LF] {};
        assert shClass(shMapC[M]) == zero[ShExt(LE,LF,0)];
    ",
    );
    assert_eq!(
        stdout(result),
        concat!(
            "setmonomaps[A] = 1\nsetmonomaps[B] = 0\nsetmonoorbitals[A] = []\n",
            "setmonochars[B] = {(1, 0): 1}\n",
        )
    );
}

#[test]
fn bad_transports_and_invalid_orbital_types_are_rejected() {
    for (tail, expected) in [
        ("setmonodromy Bad[F,G,c];", "is not bijective"),
        ("setmonoat[M,1];", "outside the image"),
        ("setmonostabilizer[M,2];", "outside the source fiber"),
        ("shorbital A[M,LF,LG,4];", "outside the product action"),
        (
            "shorbital A[M,LG,LF,0];",
            "must be a linearization of set family 'F'",
        ),
        (
            "shlinearize L3[G,3]; shorbital A[M,LF,L3,0];",
            "the same prime field",
        ),
        ("setmonodromy M[F,G,a];", "already defined"),
    ] {
        let result = interpret(&format!("{DISCONNECTED}{tail}"));
        let error = diagnostic(result);
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn large_product_actions_respect_the_storage_limit() {
    let result = interpret(
        r"
        space X {p};
        setfamily F[X] {p:700} {};
        setmonodromy M[F,F,p];
    ",
    );
    let error = diagnostic(result);
    assert!(
        error.contains("joint monodromy exceeds its storage limit"),
        "{error}"
    );
}
