mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn natural_maps_act_on_roos_cochains_and_their_classes() {
    let program = r"
        space X {a<b};
        sheaf F[X,3] {a:1,b:1} {a<b:[-1]};
        sheaf G[X,3] {a:1,b:1} {a<b:[1]};
        shmap f[F,G] {a:[-1],b:[1]};
        let x = shTerm[F,[a],0];
        let fx = shApply[f](x);
        assert fx == shTerm[G,[a],0,-1];
        assert shD(fx) == shApply[f](shD(x));
        assert shApply[f](shTerm[F,[a,b],0,2]) == shTerm[G,[a,b],0,2];
        let section = shgen[F,0,0];
        let pushed = shApply[f](section);
        assert pushed == shgen[G,0,0];
        assert pushed == shClass(shApply[f](shRep(section)));
        let push = \v : ShH(F,0) => shApply[f](v);
        assert push(shScale[-1](section)) == shScale[-1](pushed);
        show pushed;
    ";
    assert_eq!(stdout(interpret(program)), "ShH(G, 0){g0}\n");
}

#[test]
fn signed_total_differential_embedding_projection_and_functions() {
    let program = r"
        space X {a<b};
        sheaf F[X,3] {a:1,b:1} {a<b:[1]};
        shmap id[F,F] {a:[1],b:[1]};
        shcomplex K {F,id,F};
        let a = shTerm[F,[a],0];
        let ab = shTerm[F,[a,b],0];
        let x = shTotEmbed[K,0](a);
        let vx = shTotEmbed[K,1](a);
        let hx = shTotEmbed[K,0](shScale[-1](ab));
        assert shD(x) == add(hx,vx);
        assert shTotPart[0](shD(x)) == shScale[-1](ab);
        assert shTotPart[1](shD(x)) == a;
        let vv = shTotEmbed[K,1](shScale[-1](ab));
        assert shD(shTotEmbed[K,0](ab)) == vv;
        assert shD(vx) == vv;
        assert shD(shD(x)) == zero[ShTot(K,2)];
        assert shClass(shD(x)) == zero[ShHyper(K,1)];
        assert shTotPart[0](vx) == zero[ShC(F,1)];
        assert add(shScale[-1](vx),vx) == zero[ShTot(K,1)];
        let select = \v : ShTot(K,1) => shTotPart[1](v);
        assert select(shD(x)) == a;
        card (\v : ShTot(K,0) => S.hi);
        card (\v : ShHyper(K,1) => S.hi);
        shhyperdim[K,0];
        shhyperdim[K,1];
        shhyperdim[K,2];
        show shD(x);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "card = 9 / 9\n",
            "card = 1 / 1\n",
            "shhyperdim[K, 0] = 0\n",
            "shhyperdim[K, 1] = 0\n",
            "shhyperdim[K, 2] = 0\n",
            "ShTot(K, 1){2*q0:[a<b]#0 + q1:[a]#0}\n",
        )
    );
}

#[test]
fn sphere_extension_distinguishes_total_hypercohomology() {
    assert_eq!(
        stdout(example("sheaf-hypercohomology")),
        concat!(
            "shhyperdim[KTwisted, 0] = 1\n",
            "shhyperdim[KTwisted, 1] = 0\n",
            "shhyperdim[KTwisted, 2] = 0\n",
            "shhyperdim[KTwisted, 3] = 1\n",
            "shhyperdim[KSplit, 0] = 1\n",
            "shhyperdim[KSplit, 1] = 1\n",
            "shhyperdim[KSplit, 2] = 1\n",
            "shhyperdim[KSplit, 3] = 1\n",
        )
    );
}

#[test]
fn multiple_stalk_degrees_and_large_empty_components() {
    let program = r"
        space Empty {};
        sheaf E[Empty,2] {} {};
        shmap e[E,E] {};
        shcomplex EmptyK {E,e,E,e,E};
        shhyperdim[EmptyK,0];
        shhyperdim[EmptyK,2];
        shhyperdim[EmptyK,10000000];
        assert shD(zero[ShTot(EmptyK,103)]) == zero[ShTot(EmptyK,104)];
        sheaf F[One,3] {pt:1} {};
        shmap d[F,F] {pt:[0]};
        shcomplex K {F,d,F,d,F};
        let u = shTotEmbed[K,2](shTerm[F,[pt],0]);
        assert shTotPart[2](u) == shTerm[F,[pt],0];
        assert shTotPart[1](u) == zero[ShC(F,1)];
        assert shClass(u) == shhypergen[K,2,0];
        assert shRep(shClass(u)) == u;
        shhyperdim[K,2];
        shhyperdim[K,3];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shhyperdim[EmptyK, 0] = 0\nshhyperdim[EmptyK, 2] = 0\nshhyperdim[EmptyK, 10000000] = 0\nshhyperdim[K, 2] = 1\nshhyperdim[K, 3] = 0\n"
    );
}

#[test]
fn validate_naturality_and_every_stalk_composite() {
    let setup = r"
        space X {a<b};
        space Y {a<b};
        sheaf F[X,3] {a:1,b:1} {a<b:[1]};
        sheaf G[X,3] {a:1,b:1} {a<b:[1]};
        sheaf W[Y,3] {a:1,b:1} {a<b:[1]};
        sheaf B[X,2] {a:1,b:1} {a<b:[1]};
        shmap id[F,F] {a:[1],b:[1]};
        shmap f[F,G] {a:[1],b:[1]};
        shcomplex K {F,f,G};
    ";
    for (statement, reason) in [
        (
            "shmap n[F,F] {a:[1],b:[2]};",
            "is not natural on cover 'a<b'",
        ),
        ("shmap m[F,W] {a:[1],b:[1]};", "same named base and prime"),
        ("shmap m[F,B] {a:[1],b:[1]};", "same named base and prime"),
        (
            "shmap m[F,F] {a:[1],a:[1]};",
            "duplicate sheaf map component",
        ),
        ("shmap m[F,F] {a:[1]};", "missing sheaf map component"),
        (
            "shmap m[F,F] {a:[1,2],b:[1]};",
            "requires 1 rows and 1 columns",
        ),
        ("shcomplex M {};", "alternating maps and sheaves"),
        (
            "shcomplex M {G,id,F};",
            "has endpoints F -> F, expected G -> F",
        ),
        (
            "shcomplex M {F,id,F,f,G};",
            "composite 'id then f' is nonzero at 'a'",
        ),
        (
            "show shApply[f](zero[ShC(G,0)]);",
            "requires source sheaf 'F'",
        ),
        (
            "show shTotEmbed[K,2](shTerm[F,[a],0]);",
            "stalk degree 2 is out of range",
        ),
        (
            "show shTotEmbed[K,0](shTerm[G,[a],0]);",
            "requires sheaf 'F'",
        ),
        (
            "show shTotPart[1](zero[ShTot(K,0)]);",
            "negative Roos degree",
        ),
        (
            "show shClass(shTotEmbed[K,0](shTerm[F,[a],0]));",
            "requires a closed sheaf cochain",
        ),
        ("show shhypergen[K,1,0];", "generator 0 is out of range"),
    ] {
        let result = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(result.contains(reason), "{statement}\n{result}");
    }
}
