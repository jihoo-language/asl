mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn diamond_yoneda_products_are_nonzero_and_opposite() {
    assert_eq!(
        stdout(example("sheaf-yoneda")),
        concat!(
            "shExtDim[A, B, 1] = 1\nshExtDim[A, C, 1] = 1\n",
            "shExtDim[A, D, 1] = 0\nshExtDim[A, D, 2] = 1\n",
            "ShExt(A, D, 2){2*g0}\nShExt(A, D, 2){g0}\n",
            "ShHomC(A, D, 2){2*[a<b<d]#0->0}\ncard = 3 / 3\n",
        ),
    );
}

#[test]
fn matrix_pre_and_post_composition_use_the_dg_signs() {
    let program = r"
        space Chain {a<b<c<d};
        sheaf F[Chain,5] {a:2,b:2,c:2,d:2} {a<b:[1,2,0,1],b<c:[1,0,3,1],c<d:[1,2,1,1]};
        sheaf G[Chain,5] {a:2,b:2,c:2,d:2} {a<b:[2,1,4,3],b<c:[1,4,0,1],c<d:[1,0,3,1]};
        let x=shHomTerm[F,G,[a],0,1];
        let y=shHomTerm[F,G,[b],1,0];
        assert shD(x) == add(add(shHomTerm[F,G,[a,b],0,0,-1],shHomTerm[F,G,[a,b],0,1,-3]),add(
            add(shHomTerm[F,G,[a,c],0,0,-13],shHomTerm[F,G,[a,c],0,1,-3]),
            add(shHomTerm[F,G,[a,d],0,0,-13],shHomTerm[F,G,[a,d],0,1,-42])));
        assert shD(shD(y)) == zero[ShHomC(F,G,2)];
        let f=shHomTerm[F,G,[a,b],0,1];
        let g=shHomTerm[G,F,[b,c],1,0];
        let h=shHomTerm[F,G,[c,d],0,1];
        assert shExtCompose(g,f) == shHomTerm[F,F,[a,b,c],0,0,-1];
        assert shD(shExtCompose(g,f)) == add(shExtCompose(shD(g),f),shScale[-1](shExtCompose(g,shD(f))));
        assert shExtCompose(h,shExtCompose(g,f)) == shExtCompose(shExtCompose(h,g),f);
        assert shExtCompose(shHomId[G],f) == f;
        assert shExtCompose(f,shHomId[F]) == f;
        assert shD(shHomId[F]) == zero[ShHomC(F,F,1)];
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn constant_source_recovers_roos_and_natural_maps() {
    let program = r"
        space Circle {a<x,a<z,b<x,b<y,c<y,c<z};
        sheaf One3[Circle,3] {a:1,b:1,c:1,x:1,y:1,z:1} {
            a<x:[1],a<z:[1],b<x:[1],b<y:[1],c<y:[1],c<z:[1]
        };
        sheaf U[Circle,3] {a:2,b:2,c:2,x:2,y:2,z:2} {
            a<x:[1,1,0,1],a<z:[1,0,0,1],b<x:[1,0,0,1],
            b<y:[1,0,0,1],c<y:[1,0,0,1],c<z:[1,0,0,1]
        };
        shmap section[One3,U] {a:[1,0],b:[1,0],c:[1,0],x:[1,0],y:[1,0],z:[1,0]};
        assert shD(shMapC[section]) == zero[ShHomC(One3,U,1)];
        assert shClass(shMapC[section]) == shExtGen[One3,U,0,0];
        let id=\f: ShExt(One3,U,0) => f;
        assert id(shClass(shMapC[section])) == shClass(shMapC[section]);
        shdim[U,0];
        shdim[U,1];
        shExtDim[One3,U,0];
        shExtDim[One3,U,1];
        shExtDim[One3,U,2];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shdim[U, 0] = 1\nshdim[U, 1] = 1\nshExtDim[One3, U, 0] = 1\nshExtDim[One3, U, 1] = 1\nshExtDim[One3, U, 2] = 0\n",
    );
}

#[test]
fn malformed_hom_inputs_are_rejected() {
    let setup = r"
        space Chain {a<b};
        sheaf A[Chain,3] {a:1,b:1} {a<b:[1]};
        sheaf B[Chain,3] {a:1,b:1} {a<b:[1]};
        sheaf A5[Chain,5] {a:1,b:1} {a<b:[1]};
        sheaf Other[One,3] {pt:1} {};
    ";
    for (program, expected) in [
        ("shExtDim[A,A5,0];", "same named base and prime field"),
        ("shExtDim[A,Other,0];", "same named base and prime field"),
        ("show shHomTerm[A,B,[b,a],0,0];", "strict increasing chain"),
        (
            "show shHomTerm[A,B,[a],1,0];",
            "outside the endpoint stalks",
        ),
        (
            "show shClass(shHomTerm[A,B,[a],0,0]);",
            "requires a closed sheaf cochain",
        ),
        (
            "show shExtCompose(shHomId[A],shHomId[B]);",
            "matching middle sheaf",
        ),
        (
            "show shExtCompose(shClass(shHomId[A]),shHomId[A]);",
            "two cochains or two classes",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{program}")));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
