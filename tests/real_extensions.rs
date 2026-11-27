mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn hidden_products_compare_to_the_stable_page_and_have_exact_null_witnesses() {
    assert_eq!(
        stdout(example("real-rho-extensions")),
        concat!(
            "RMBar(2, 3, 0){rho^1*b0}\n",
            "rhodepth = 0\n",
            "rhodepth = 0\n",
            "rhodepth = 1\n",
            "rhoorder = 1\n",
            "rhohidden = hidden (0 + 0 < 1)\n",
            "RhoBock(4, 2, 3, 0, 1){g0}\n",
            "rhohidden = hidden (0 + 0 < 2)\n",
            "rhohidden = hidden (0 + 0 < 2)\n",
        )
    );
}

#[test]
fn living_bar_coordinates_are_first_class_and_division_omits_dying_bars() {
    let source = r"
        let h0 = extClass(rmcword[0,0]{([],[1])});
        let h1 = extClass(rmcword[0,0]{([1],[])});
        let sh = rhoShift[1](h1);
        let x = add(h0,sh);
        let bars = add(rmbarTerm[1,1,0,0],rmbarTerm[1,1,0,1]);
        assert rmbarCoords(x) == bars;
        assert rmbarClass(bars) == x;
        assert extClass(rmbarLift(bars)) == x;
        assert rmbarCoords(zero[RMExt(1,1,0)]) == zero[RMBar(1,1,0)];
        assert rmbarClass(rmbarCoords(zero[RMExt(1,1,0)])) == zero[RMExt(1,1,0)];
        card(\r : RMBar(1,1,0) => S.hi);
        card(\r : RMBar(1,0,-1) => S.hi);
        let y = rhoShift[1](x);
        assert rhoDivide[1](y) == sh;
        assert rhoShift[1](rhoDivide[1](y)) == rhoShift[1](x);
        assert rhoDivide[2](y) == h1;
        rhodepth(x);
        rhodepth(y);
        rhoorder(y);
        rhodepth(rhoShift[1](h0));
        let rho = extClass(rmcword[1,0]{});
        rhohidden(rho,h1);
        rhohidden(rho,h0);
        assert rhoProductNull[0](rho,h0) == rmcword[0,1]{};
        let u = rhoShift[20](h1);
        assert rhoDivide[20](u) == h1;
        assert rhoEinf[20](u) == rhobockgen[21,1,-18,-19,20,0];
        rhodepth(u);
        assert rhoDivide[100](zero[RMExt(1,-3,-4)]) == zero[RMExt(1,97,96)];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 4 / 4\n",
            "card = 2 / 2\n",
            "rhodepth = 0\n",
            "rhodepth = 2\n",
            "rhoorder = infinite\n",
            "rhodepth = infinite\n",
            "rhohidden = visible (1 + 0 = 1)\n",
            "rhohidden = zero\n",
            "rhodepth = 20\n",
        )
    );
}

#[test]
fn hidden_products_can_be_divided_in_their_first_bar_query() {
    for n in [1, 2, 4] {
        let source = format!(
            "let a=extClass(rmcword[0,0]{{([],[{n}])}});
             let b=extClass(rmcword[0,0]{{([],[{}])}});
             let p=yoneda(a,b);
             let q=rhoDivide[{n}](p);
             assert rhoShift[{n}](q)==p;
             assert extClass(extRep(q))==q;
             rhodepth(p);",
            2 * n,
        );
        assert_eq!(stdout(interpret(&source)), format!("rhodepth = {n}\n"));
    }
}

#[test]
fn compact_real_rows_compute_the_depth_four_product_and_its_null() {
    let source = r"
        let t4 = extClass(rmcword[0,0]{([],[4])});
        let t8 = extClass(rmcword[0,0]{([],[8])});
        rhohidden(t4,t8);
        let p = yoneda(t4,t8);
        rmbars[2,12];
        let q = rhoDivide[4](p);
        assert rhoShift[4](q) == p;
        let witness = rhoProductNull[4](t4,t8);
        let reconstructed = rhoShift[4](rmbarLift(rmbarCoords(q)));
        assert cobarD(witness) == add(concat(extRep(t4),extRep(t8)),reconstructed);
        rhodepth(p);
        rhoorder(p);
        rhoorder(q);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rhohidden = hidden (0 + 0 < 4)\n",
            "rmbars[2, 12] = [(2, 1), (4, 2), (16, 8), (20, 4)]\n",
            "rhodepth = 4\n",
            "rhoorder = 4\n",
            "rhoorder = 8\n",
        )
    );
}

#[test]
fn extensions_validate_depth_stability_and_types() {
    for (source, reason) in [
        ("show rmbarTerm[1,0,-1,0];", "rho bar does not live"),
        (
            "show rhoDivide[1](extClass(rmcword[0,0]{([],[1])}));",
            "rho depth 0",
        ),
        (
            "show rhoEinf[1,0](extClass(rmcword[0,0]{([],[1])}));",
            "requires page 2 or later",
        ),
        (
            "show rhoEinf[1](extClass(rmcword[0,0]{([],[1])}));",
            "requires depth at least 1",
        ),
        (
            "show rhoProductNull[0](rmextgen[0,0,0,0],rmextgen[0,0,0,0]);",
            "positive total cobar length",
        ),
        ("show rmbarClass(rmextgen[0,0,0,0]);", "RMBar coordinates"),
        ("rhodepth(cword{});", "an RMExt class or RMBar coordinates"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
