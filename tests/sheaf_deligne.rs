mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn smart_truncation_uses_the_kernel_and_refactors_the_incoming_map() {
    let program = r"
        sheaf F[One,5] {pt:1} {};
        sheaf G[One,5] {pt:2} {};
        shmap m[F,G] {pt:[2,2]};
        shmap n[G,F] {pt:[1,-1]};
        shcomplex K {F,m,G,n,F};
        shtruncate Low[K,0];
        shtruncate Mid[K,1];
        shtruncate All[K,1000];
        shhyperdim[Low,0];
        shhyperdim[Mid,0];
        shhyperdim[Mid,1];
        shhyperdim[All,2];
        let x = shTotEmbed[Mid,0](shTerm[F,[pt],0]);
        assert shD(x) == shTotEmbed[Mid,1](shTerm[Mid__1,[pt],0,2]);
        assert shD(shD(x)) == zero[ShTot(Mid,2)];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shhyperdim[Low, 0] = 0\nshhyperdim[Mid, 0] = 0\nshhyperdim[Mid, 1] = 0\nshhyperdim[All, 2] = 0\n",
    );
}

#[test]
fn a_terminal_stalk_kernel_inherits_the_actual_restriction() {
    let program = r"
        space X {a<b};
        sheaf F[X,5] {a:2,b:2} {a<b:[2,0,0,2]};
        sheaf G[X,5] {a:1,b:1} {a<b:[2]};
        shmap d[F,G] {a:[1,-1],b:[1,-1]};
        shcomplex K {F,d,G};
        shtruncate T[K,0];
        assert shD(shTerm[T__0,[a],0]) == shTerm[T__0,[a,b],0,-2];
        shhyperdim[T,0];
        shhyperdim[T,1];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shhyperdim[T, 0] = 1\nshhyperdim[T, 1] = 0\n",
    );
}

#[test]
fn the_torus_suspension_has_complementary_deligne_dimensions() {
    let expected = concat!(
        "betti(Susp) = [1, 0, 2, 1]\n",
        "shdim[Whole, 0] = 1\nshhyperdim[Lower, 0] = 1\nshhyperdim[Upper, 3] = 1\n",
        "shdim[Whole, 1] = 0\nshhyperdim[Lower, 1] = 0\nshhyperdim[Upper, 2] = 0\n",
        "shdim[Whole, 2] = 2\nshhyperdim[Lower, 2] = 2\nshhyperdim[Upper, 1] = 2\n",
        "shdim[Whole, 3] = 1\nshhyperdim[Lower, 3] = 1\nshhyperdim[Upper, 0] = 1\n",
        "shedim[Upper, 2, 0, 1] = 4\n",
        "shedim[Upper, 2, 2, 0] = 2\n",
        "shedim[Upper, 3, 0, 1] = 2\n",
        "shedim[Upper, 3, 2, 0] = 0\n",
    );
    assert_eq!(stdout(example("sheaf-deligne")), expected);
    let fixture = include_str!("../examples/sheaf-deligne.asl");
    for prime in [2, 5] {
        let fixture = fixture
            .replace("[Link,3]", &format!("[Link,{prime}]"))
            .replace("[Susp,3]", &format!("[Susp,{prime}]"));
        assert_eq!(stdout(interpret(&fixture)), expected);
    }
}

#[test]
fn single_step_deligne_checks_open_inclusions_and_preserves_discarded_names() {
    let program = r"
        space X {a<b};
        space D {a,b};
        let id = \x : X => x;
        let low = tab One -> X {pt=>X.a};
        let swallow = \x : X => X.b;
        let nonreflecting = tab D -> X {a=>X.a,b=>X.b};
        sheaf F[X,3] {a:1,b:1} {a<b:[1]};
        sheaf P[One,3] {pt:1} {};
        sheaf G[D,3] {a:1,b:1} {};
        sheaf K__1[X,3] {a:0,b:0} {a<b:[]};
        shdeligne K[id,F,0];
        shhyperdim[K,0];
        shhyperdim[K,1];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shhyperdim[K, 0] = 1\nshhyperdim[K, 1] = 0\n"
    );
    for (statement, reason) in [
        ("shdeligne B[low,P,0];", "requires an upward-closed image"),
        (
            "shdeligne B[swallow,F,0];",
            "requires an injective open embedding",
        ),
        (
            "shdeligne B[nonreflecting,G,0];",
            "requires an order-reflecting open embedding",
        ),
        (
            "shdeligne B[id,P,1];",
            "requires a sheaf on the open embedding's named source",
        ),
        (
            "shdeligne K[id,F,0];",
            "sheaf complex 'K' is already defined",
        ),
        (
            "shtruncate K[K,1000];",
            "sheaf complex 'K' is already defined",
        ),
        ("shtruncate T[Missing,0];", "unknown sheaf complex"),
    ] {
        let error = diagnostic(interpret(&format!("{program}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
