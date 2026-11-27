mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn torus_link_connecting_classes_and_delignes_cutoff() {
    let output = stdout(example("sheaf-support"));
    assert_eq!(
        output.lines().take(6).collect::<Vec<_>>(),
        [
            "shlocaldim[Costalk, pt, relative, 0] = 0",
            "shlocaldim[Costalk, pt, relative, 1] = 0",
            "shlocaldim[Costalk, pt, relative, 2] = 2",
            "shlocaldim[Costalk, pt, relative, 3] = 1",
            "shlocaldim[PerverseCostalk, pt, relative, 2] = 0",
            "shlocaldim[PerverseCostalk, pt, relative, 3] = 1",
        ]
    );
    assert_eq!(output.lines().count(), 9);
    let program = std::fs::read_to_string(format!(
        "{}/examples/sheaf-support.asl",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap();
    let variants = program.replace(",3] {", ",5] {");
    for program in [program, variants] {
        let assertions = r"
            assert shLocalConnect(a) == shlocalgen[Costalk,pt,relative,2,0];
            assert shLocalConnect(b) == shlocalgen[Costalk,pt,relative,2,1];
            assert shLocalConnect(add(a,shScale[2](b))) == add(shLocalConnect(a),shScale[2](shLocalConnect(b)));
            assert shLocalConnect(top) == shlocalgen[PerverseCostalk,pt,relative,3,0];
            assert shLocalConnect(shlocalgen[PerverseCostalk,pt,punctured,1,0]) == zero[ShLocal(PerverseCostalk,pt,relative,2)];
            shhyperdim[Costalk,2];
            shhyperdim[PerverseCostalk,3];
        ";
        let output = stdout(interpret(&(program + assertions)));
        assert!(
            output.ends_with("shhyperdim[Costalk, 2] = 2\nshhyperdim[PerverseCostalk, 3] = 1\n")
        );
    }
}

#[test]
fn local_sequence_uses_signed_coordinates_and_preserves_class_operations() {
    let source = r"
        space Fork {v,a,b,v<a,v<b};
        let i = tab One -> Fork {pt=>Fork.v};
        sheaf F[Fork,5] {v:1,a:1,b:1} {v<a:[1],v<b:[1]};
        shextra K[i,F];
        let u = shLocalClass[K,pt,punctured](shTerm[F,[a],0]);
        let v = shLocalClass[K,pt,punctured](shTerm[F,[b],0]);
        let r = shLocalClass[K,pt,relative](shTerm[F,[v,a],0]);
        assert shLocalConnect(u) == r;
        assert shLocalConnect(v) == shScale[-1](r);
        assert shLocalConnect(add(u,v)) == zero[ShLocal(K,pt,relative,1)];
        let diagonal = shLocalClass[K,pt,absolute](add(shTerm[F,[v],0],add(shTerm[F,[a],0],shTerm[F,[b],0])));
        assert shLocalRestrict(diagonal) == add(u,v);
        assert shLocalConnect(shLocalRestrict(diagonal)) == zero[ShLocal(K,pt,relative,1)];
        assert shLocalForget(r) == zero[ShLocal(K,pt,absolute,1)];
        assert shLocalClass[K,pt,relative](shLocalRep(r)) == r;
        let localnegate = \x : ShLocal(K,pt,punctured,0) => shScale[-1](x);
        assert shLocalConnect(localnegate(u)) == shScale[-1](r);
        card (\x : ShLocal(K,pt,punctured,0) => S.hi);
        shlocaldim[K,pt,absolute,0];
        shlocaldim[K,pt,punctured,0];
        shlocaldim[K,pt,relative,1];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 25 / 25\n",
            "shlocaldim[K, pt, absolute, 0] = 1\n",
            "shlocaldim[K, pt, punctured, 0] = 2\n",
            "shlocaldim[K, pt, relative, 1] = 1\n",
        )
    );
}

#[test]
fn multiple_closed_points_restrict_actual_relative_coordinates() {
    let source = r"
        space Fork {v,a,b,v<a,v<b};
        space Z {bottom<end};
        let i = tab Z -> Fork {bottom=>Fork.v,end=>Fork.a};
        sheaf F[Fork,3] {v:1,a:1,b:1} {v<a:[1],v<b:[1]};
        shextra K[i,F];
        let x = shTerm[K__0,[bottom],0];
        let y = shTerm[K__0,[bottom],1];
        assert shD(x) == zero[ShC(K__0,1)];
        assert shD(y) == shTerm[K__0,[bottom,end],0,-1];
        shlocaldim[K,bottom,relative,0];
        shlocaldim[K,end,relative,0];
        shlocaldim[K,end,punctured,0];
        let c = shlocalgen[K,end,relative,0,0];
        assert shLocalForget(c) == shlocalgen[K,end,absolute,0,0];
        assert shLocalRestrict(shLocalForget(c)) == zero[ShLocal(K,end,punctured,0)];
        shcohomology H[K,0];
        shdim[H,0];
        space Empty {};
        let e = tab Empty -> Fork {};
        shextra Vac[e,F];
        shhyperdim[Vac,0];
        let identity = \x : Fork => x;
        shextra All[identity,F];
        shlocaldim[All,v,relative,0];
        shlocaldim[All,v,relative,1];
        shlocaldim[All,a,punctured,0];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "shlocaldim[K, bottom, relative, 0] = 0\n",
            "shlocaldim[K, end, relative, 0] = 1\n",
            "shlocaldim[K, end, punctured, 0] = 0\n",
            "shdim[H, 0] = 0\n",
            "shhyperdim[Vac, 0] = 0\n",
            "shlocaldim[All, v, relative, 0] = 1\n",
            "shlocaldim[All, v, relative, 1] = 0\n",
            "shlocaldim[All, a, punctured, 0] = 0\n",
        )
    );
}

#[test]
fn supports_check_closed_embeddings_and_local_classes() {
    let source = r"
        space Fork {v,a,b,v<a,v<b};
        let i = tab One -> Fork {pt=>Fork.v};
        sheaf F[Fork,5] {v:1,a:1,b:1} {v<a:[1],v<b:[1]};
        shextra K[i,F];
    ";
    for (statement, reason) in [
        (
            "let bad = tab One -> Fork {pt=>Fork.a}; shextra B[bad,F];",
            "downward-closed image",
        ),
        (
            "space Pair {c,d}; let bad = tab Pair -> Fork {c=>Fork.v,d=>Fork.v}; shextra B[bad,F];",
            "injective closed embedding",
        ),
        (
            "space Pair {c,d}; let bad = tab Pair -> Fork {c=>Fork.v,d=>Fork.a}; shextra B[bad,F];",
            "order-reflecting closed embedding",
        ),
        (
            "show shLocalClass[K,pt,relative](shTerm[F,[a],0]);",
            "must vanish on the punctured open",
        ),
        (
            "show shLocalClass[K,pt,absolute](shTerm[F,[a],0]);",
            "requires a closed sheaf cochain",
        ),
        (
            "show shLocalConnect(zero[ShLocal(K,pt,relative,1)]);",
            "requires a punctured local class",
        ),
        ("shlocaldim[K,a,relative,0];", "unknown point 'a'"),
    ] {
        let error = diagnostic(interpret(&format!("{source}{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
