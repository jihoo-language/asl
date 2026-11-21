#[allow(dead_code)]
mod common;

use std::process::Command;

use common::{diagnostic, interpret, stdout};

fn generated(order: usize, prime: usize, split: bool) -> String {
    let result = Command::new("awk")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(["-v", &format!("r={order}"), "-v", &format!("prime={prime}")])
        .args(["-v", if split { "untwisted=1" } else { "untwisted=0" }])
        .args(["-f", "examples/sheaf-transgression.awk"])
        .output()
        .unwrap();
    stdout(result)
}

#[test]
fn e0_is_signed_stalk_differential_and_e1_is_roos() {
    let program = r"
        space X {a<b};
        sheaf F[X,3] {a:1,b:1} {a<b:[1]};
        shmap id[F,F] {a:[1],b:[1]};
        shcomplex Contractible {F,id,F};
        let a = shTerm[F,[a],0];
        let ab = shTerm[F,[a,b],0];
        let base = shTotEmbed[Contractible,0](ab);
        let vertical = shTotEmbed[Contractible,1](ab);
        let c = shEClass[0,1](base);
        let next = shED(c);
        assert next == shScale[-1](shEClass[0,1](vertical));
        assert shERep(c) == base;
        assert shENext(next) == zero[ShE(Contractible,1,1,1)];
        assert shED(next) == zero[ShE(Contractible,0,1,2)];
        shedim[Contractible,0,1,0];
        shedim[Contractible,1,1,0];

        shcomplex Horizontal {F};
        let u = shTotEmbed[Horizontal,0](a);
        let v = shEClass[1,0](u);
        assert v == shENext(shEClass[0,0](u));
        assert shED(v) == shScale[-1](shegen[Horizontal,1,1,0,0]);
        assert shED(shED(v)) == zero[ShE(Horizontal,1,2,0)];
        assert shENext(shED(v)) == zero[ShE(Horizontal,2,1,0)];
        assert shEClass[1,0](shERep(v)) == v;
        assert shELiftTo[1,0](u) == u;
        let s = shegen[Horizontal,2,0,0,0];
        assert shENext(shEPrev(s)) == s;
        assert shED(s) == zero[ShE(Horizontal,2,2,-1)];
        assert shERep(shED(s)) == zero[ShTot(Horizontal,1)];
        let stable = shegen[Horizontal,1000,0,0,0];
        assert shENext(shEPrev(stable)) == stable;
        shedim[Horizontal,1000,0,0];
        let negate = \x : ShE(Horizontal,2,0,0) => shScale[-1](x);
        assert add(negate(s),s) == zero[ShE(Horizontal,2,0,0)];
        card (\x : ShE(Horizontal,2,0,0) => S.hi);
        shedim[Horizontal,2,0,0];
        shedim[Horizontal,2,1,0];
        shedim[Horizontal,2,2,-1];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "shedim[Contractible, 0, 1, 0] = 1\n",
            "shedim[Contractible, 1, 1, 0] = 0\n",
            "shedim[Horizontal, 1000, 0, 0] = 1\n",
            "card = 3 / 3\n",
            "shedim[Horizontal, 2, 0, 0] = 1\n",
            "shedim[Horizontal, 2, 1, 0] = 0\n",
            "shedim[Horizontal, 2, 2, -1] = 0\n",
        ),
    );
}

#[test]
fn twisted_spheres_support_real_differentials_through_page_six() {
    for (r, p) in [(2, 2), (3, 5), (4, 3), (6, 3)] {
        let mut source = generated(r, p, false);
        source.push_str(&format!(
            "let a = shTotEmbed[K,{q}](shTotPart[{q}](shERep(v2)));\n
             let lift = shELiftTo[{r},0](a);\n
             assert shEClass[{r},0](lift) == v{r};\n
             assert shELiftTo[{r},0](lift) == lift;\n
             assert shENext(shEPrev(v{r})) == v{r};\n
             assert hit == shScale[-1](shegen[K,{r},{r},0,0]);\n",
            q = r - 1
        ));
        let output = stdout(interpret(&source));
        let coefficient = if p == 2 {
            String::new()
        } else {
            format!("{}*", p - 1)
        };
        assert_eq!(
            output,
            format!(
                "shedim[K, 2, 0, {q}] = 1\nShE(K, {r}, {r}, 0){{{coefficient}g0}}\n\
                shedim[K, {s}, 0, {q}] = 0\nshedim[K, {s}, {r}, 0] = 0\n",
                q = r - 1,
                s = r + 1
            )
        );
    }
}

#[test]
fn the_split_sphere_survives_the_same_pages() {
    let mut source = generated(4, 3, true);
    source.push_str(
        "let survivor = shENext(v4);\n
         assert shENext(shEPrev(survivor)) == survivor;\n
         assert shEClass[5,0](shERep(survivor)) == survivor;",
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "shedim[K, 2, 0, 3] = 1\n",
            "ShE(K, 4, 4, 0){0}\n",
            "shedim[K, 5, 0, 3] = 1\n",
            "shedim[K, 5, 4, 0] = 1\n",
        ),
    );
}

#[test]
fn page_diagnostics_and_true_obstructions() {
    let source = generated(3, 3, false);
    for (statement, reason) in [
        ("show shENext(v3);", "shENext requires a cycle for d_3"),
        (
            "let a = shTotEmbed[K,2](shTotPart[2](shERep(v2))); show shELiftTo[4,0](a);",
            "cannot lift the prescribed leading part to page 4",
        ),
        ("show shEClass[4,0](shERep(v3));", "page cycle condition"),
        (
            "show shEClass[2,1](shERep(v2));",
            "lies below the requested sheaf filtration",
        ),
        (
            "show shELiftTo[2,1](shERep(v2));",
            "lies below the requested sheaf filtration",
        ),
        (
            "show shEPrev(zero[ShE(K,0,0,0)]);",
            "requires a positive page",
        ),
        ("show shegen[K,1,0,9,0];", "generator 0 is out of range"),
        (
            "show shERep(zero[ShE(K,0,0,-1)]);",
            "no representable nonnegative total degree",
        ),
        ("show shED(zero[ShTot(K,0)]);", "requires a ShE class"),
    ] {
        let error = diagnostic(interpret(&format!("{source}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
    let overflow = diagnostic(interpret(&format!(
        "{source}\nshow shENext(zero[ShE(K,{},0,0)]);",
        usize::MAX
    )));
    assert!(overflow.contains("page is out of range"), "{overflow}");
    let overflow = diagnostic(interpret(&format!(
        "{source}\nshow shED(zero[ShE(K,{},1,0)]);",
        usize::MAX
    )));
    assert!(overflow.contains("bidegree is out of range"), "{overflow}");
}
