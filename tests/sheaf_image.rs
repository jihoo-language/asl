mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn higher_restrictions_act_on_actual_cocycles_and_compose() {
    assert_eq!(
        stdout(example("sheaf-direct-image")),
        concat!(
            "shdim[R0, 0] = 2\n",
            "shdim[R1, 0] = 2\n",
            "shdim[R2, 0] = 0\n",
            "ShC(R1, 1){4*[low<mid]#0 + 4*[low<top]#0}\n",
            "shdim[Pulled, 1] = 3\n",
        )
    );
}

#[test]
fn identity_pushforward_has_the_original_maps_and_no_higher_stalks() {
    let source = r"
        space X {c,b,a,a<b<c};
        sheaf F[X,5] {a:1,b:1,c:1} {a<b:[2],b<c:[-1]};
        let identity = \x : X => x;
        shpull G[identity,F];
        shpush I[identity,F,0];
        shpush J[identity,F,1];
        shpush L[identity,F,2];
        assert shD(shTerm[G,[a],0]) ==
            add(shTerm[G,[a,b],0,3],shTerm[G,[a,c],0,2]);
        assert shD(shTerm[I,[a],0]) ==
            add(shTerm[I,[a,b],0,3],shTerm[I,[a,c],0,2]);
        assert shD(shTerm[I,[b],0]) ==
            add(shTerm[I,[a,b],0],shTerm[I,[b,c],0]);
        shdim[G,0];
        shdim[I,0];
        shdim[J,0];
        shdim[L,0];
        shdim[I,1];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "shdim[G, 0] = 1\nshdim[I, 0] = 1\nshdim[J, 0] = 0\nshdim[L, 0] = 0\nshdim[I, 1] = 0\n"
    );
    let invalid = diagnostic(interpret(&format!("{source}\nshow shTerm[J,[a],0];")));
    assert!(invalid.contains("dimension 0"), "{invalid}");
}

#[test]
fn point_projection_and_empty_or_equal_opens_have_exactly_the_right_dimensions() {
    let source = r"
        space Circle {a<x,a<z,b<x,b<y,c<y,c<z};
        space Chain {l<m<n};
        let constant = \x : Circle => Chain.m;
        let to_point = \x : Circle => One.pt;
        sheaf F[Circle,3] {a:1,b:1,c:1,x:1,y:1,z:1}
            {a<x:[1],a<z:[1],b<x:[1],b<y:[1],c<y:[1],c<z:[1]};
        sheaf H[Chain,3] {l:1,m:1,n:1} {l<m:[-1],m<n:[0]};
        shpush R[constant,F,1];
        assert shD(shTerm[R,[l],0]) == shTerm[R,[l,m],0,-1];
        shdim[R,0];
        shpush Global0[to_point,F,0];
        shpush Global1[to_point,F,1];
        shpush Global2[to_point,F,2];
        shdim[Global0,0];
        shdim[Global1,0];
        shdim[Global2,0];
        shpull Constant[constant,H];
        shdim[Constant,0];
        shdim[Constant,1];
        space E {};
        let empty = tab E -> Chain {};
        shpush Z[empty,H,0];
    ";
    let wrong_endpoint = diagnostic(interpret(source));
    assert!(
        wrong_endpoint.contains("requires a sheaf on the map's named source"),
        "{wrong_endpoint}"
    );
    let source = source.replace(
        "shpush Z[empty,H,0];",
        "sheaf EF[E,3] {} {}; shpush Z[empty,EF,0]; shpull EZ[empty,H]; shdim[Z,0]; shdim[EZ,0];",
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "shdim[R, 0] = 1\n",
            "shdim[Global0, 0] = 1\n",
            "shdim[Global1, 0] = 1\n",
            "shdim[Global2, 0] = 0\n",
            "shdim[Constant, 0] = 1\n",
            "shdim[Constant, 1] = 1\n",
            "shdim[Z, 0] = 0\n",
            "shdim[EZ, 0] = 0\n",
        )
    );
}

#[test]
fn pullback_uses_the_transitive_restriction_on_a_noncover_image() {
    let source = r"
        space X {x<y};
        space Y {a<b<c};
        let f = tab X -> Y {x=>Y.a,y=>Y.c};
        sheaf F[Y,5] {a:1,b:2,c:1} {a<b:[2,1],b<c:[1,-1]};
        shpull G[f,F];
        assert shD(shTerm[G,[x],0]) == shTerm[G,[x,y],0,-1];
        shdim[G,0];
    ";
    assert_eq!(stdout(interpret(source)), "shdim[G, 0] = 1\n");
    for (statement, reason) in [
        (
            "shpull Bad[f,G];",
            "requires a sheaf on the map's named target",
        ),
        (
            "shpush Bad[f,F,1];",
            "requires a sheaf on the map's named source",
        ),
        ("shpush G[f,G,0];", "already defined"),
        ("shpull Bad[missing,F];", "unknown named map"),
        (
            "let z = Y.a; shpull Bad[z,F];",
            "requires a named monotone map",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{source}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
