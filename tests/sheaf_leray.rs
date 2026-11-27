mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn actual_leray_differential_distinguishes_two_sheaves_with_constant_higher_images() {
    let before = concat!(
        "shdim[R0, 0] = 1\nshdim[R0, 1] = 0\nshdim[R0, 2] = 1\n",
        "shdim[R1, 0] = 1\nshdim[R1, 1] = 0\nshdim[R1, 2] = 1\n",
        "shedim[K, 2, 0, 1] = 1\nshedim[K, 2, 2, 0] = 1\n",
    );
    assert_eq!(
        stdout(example("sheaf-leray")),
        format!(
            "{before}{}",
            concat!(
                "ShE(K, 2, 2, 0){g0}\n",
                "shedim[K, 3, 0, 1] = 0\nshedim[K, 3, 2, 0] = 0\n",
                "shdim[F, 0] = 1\nshdim[F, 1] = 0\nshdim[F, 2] = 0\nshdim[F, 3] = 1\n",
                "shhyperdim[K, 0] = 1\nshhyperdim[K, 1] = 0\n",
                "shhyperdim[K, 2] = 0\nshhyperdim[K, 3] = 1\n",
            ),
        ),
    );
    let twist = "x2a<x4a:[1,1,0,1]";
    let example = include_str!("../examples/sheaf-leray.asl");
    assert_eq!(example.matches(twist).count(), 1);
    let split = example.replace(twist, "x2a<x4a:[1,0,0,1]");
    assert_eq!(
        stdout(interpret(&split)),
        format!(
            "{before}{}",
            concat!(
                "ShE(K, 2, 2, 0){0}\n",
                "shedim[K, 3, 0, 1] = 1\nshedim[K, 3, 2, 0] = 1\n",
                "shdim[F, 0] = 1\nshdim[F, 1] = 1\nshdim[F, 2] = 1\nshdim[F, 3] = 1\n",
                "shhyperdim[K, 0] = 1\nshhyperdim[K, 1] = 1\n",
                "shhyperdim[K, 2] = 1\nshhyperdim[K, 3] = 1\n",
            ),
        ),
    );
}

#[test]
fn point_projection_compares_raw_cochains_and_canonical_classes() {
    let program = r"
        space Circle {a<x,a<z,b<x,b<y,c<y,c<z};
        let to_point = \x : Circle => One.pt;
        sheaf F[Circle,3] {a:1,b:1,c:1,x:1,y:1,z:1}
            {a<x:[1],a<z:[1],b<x:[1],b<y:[1],c<y:[1],c<z:[1]};
        shleray K[to_point,F];
        let a = shTerm[F,[a],0];
        let i = shLerayIncl[K](a);
        assert i == shTotEmbed[K,0](shTotPart[0](i));
        let at = \v : ShTot(K,0) => shTotPart[0](v);
        assert at(i) == shTerm[K__0,[pt],0];
        assert shD(i) == shLerayIncl[K](shD(a));
        assert shClass(shLerayIncl[K](shRep(shgen[F,0,0]))) == shhypergen[K,0,0];
        assert shClass(shLerayIncl[K](shRep(shgen[F,1,0]))) == shhypergen[K,1,0];
        assert shLerayIncl[K](zero[ShC(F,1000000)]) == zero[ShTot(K,1000000)];
        shhyperdim[K,0];
        shhyperdim[K,1];
        shhyperdim[K,2];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shhyperdim[K, 0] = 1\nshhyperdim[K, 1] = 1\nshhyperdim[K, 2] = 0\n",
    );
}

#[test]
fn empty_preimages_and_validation_use_the_named_spaces() {
    let program = r"
        space E {};
        let empty = tab E -> One {};
        sheaf EF[E,5] {} {};
        sheaf F[One,5] {pt:1} {};
        shcomplex Ordinary {F};
        shleray K[empty,EF];
        assert shLerayIncl[K](zero[ShC(EF,0)]) == zero[ShTot(K,0)];
        assert shTotPart[0](zero[ShTot(K,0)]) == zero[ShC(K__0,0)];
        shhyperdim[K,0];
        shhyperdim[K,1000];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shhyperdim[K, 0] = 0\nshhyperdim[K, 1000] = 0\n",
    );
    for (statement, reason) in [
        (
            "shleray Z[empty,F];",
            "requires a sheaf on the map's named source",
        ),
        (
            "shleray K[empty,EF];",
            "sheaf complex 'K' is already defined",
        ),
        (
            "sheaf P__0[One,5] {pt:0} {}; shleray P[empty,EF];",
            "sheaf 'P__0' is already defined",
        ),
        (
            "show shLerayIncl[Ordinary](zero[ShC(F,0)]);",
            "is not a declared Leray complex",
        ),
        (
            "show shLerayIncl[K](zero[ShC(F,0)]);",
            "requires source sheaf 'EF'",
        ),
        (
            "show shLerayIncl[K](zero[ShH(EF,0)]);",
            "requires a ShC cochain",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{program}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
